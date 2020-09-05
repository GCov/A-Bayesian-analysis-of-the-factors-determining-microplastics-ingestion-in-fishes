#### Load packages ####
library(coda)
library(R2jags)
library(mcmcplots)
library(arm)
library(ggplot2)
library(cowplot)
library(plyr)
library(lattice)
library(colorspace)
library(beepr)
library(ggridges)
library(reshape2)
library(DHARMa)
load.module("glm")

extract.post <- function(x){
  out <- data.frame(x$BUGSoutput$sims.list)
  long <- melt(out)
  long <- long[long$variable != "deviance" &
                 long$variable != "sigma" &
                 long$variable != "alpha_p" &
                 long$variable != "beta_p" &
                 long$variable != "r" &
                 long$variable != "rho" &
                 long$variable != "mu_TL" &
                 long$variable != "phi",]
  long$variable <- as.character(long$variable)
  long$variable <- as.factor(long$variable)
  long
}

#### Set up some aesthetics ####

theme1 <-
  theme_bw() +
  theme(
    text = element_text(size = 7),
    axis.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    panel.grid = element_blank()
  )

pal <- c("#c52f01",  #rust
         "#b88100",  #dark goldenrod 
         "#0076a9",  #celadon blue
         "#5a8400",  #avocado
         "#011e44")  #oxford blue

#### Set up the data ####
trophicfish <- 
  read.csv("trophicfishdata.csv", header = TRUE)

names(trophicfish)
trophicfish <- trophicfish[1:42]

summary(trophicfish)

trophicfish$year <- as.factor(trophicfish$year)

# concentration by weight

trophicfish$conc <- trophicfish$Mpsgut/trophicfish$W
summary(trophicfish)
trophicfish$IR <- trophicfish$ingest.rate/100

## Consolidate averages by study

trophicfish$fork.lengthprod <- trophicfish$fork.length*trophicfish$N
trophicfish$total.lengthprod <- trophicfish$total.length*trophicfish$N
trophicfish$Wprod <- trophicfish$W*trophicfish$N
trophicfish$GWprod <- trophicfish$GW*trophicfish$N
trophicfish$IRprod <- trophicfish$IR*trophicfish$N
trophicfish$Mpsgutprod <- trophicfish$Mpsgut*trophicfish$N
trophicfish$ SDMPsgutprod <- trophicfish$ SDMPsgut*trophicfish$N

trophicfish2 <- ddply(trophicfish, 
                      c('author', 'study.habitat', 'year', 'region', 'species',
                        'life.stage', 'family', 'genus', 'environment', 
                        'climate', 'red.list', 'feeding.type', 'feeding.habit', 
                        'TL', 'min.size', 'float.meth', 'dig.meth', 
                        'count.meth', 'polymer.meth', 'polymer.ID', 'blanks', 
                        'maj.fib', 'maj.under.one.mm', 'maj.polymer', 
                        'maj.col', 'exclude.fib'), 
                      summarise, 
                      fork.lengthprod = sum(fork.lengthprod), 
                      total.lengthprod = sum(total.lengthprod), 
                      Wprod = sum(Wprod), 
                      GWprod = sum(GWprod),
                      N = sum(N), 
                      IRprod = sum(IRprod), 
                      Mpsgutprod = sum(Mpsgutprod), 
                      SDMPsgutprod = sum(SDMPsgutprod)
)

trophicfish2$fork.length <- trophicfish2$fork.lengthprod/trophicfish2$N
trophicfish2$total.length <- trophicfish2$total.lengthprod/trophicfish2$N
trophicfish2$W <- trophicfish2$Wprod/trophicfish2$N
trophicfish2$GW <- trophicfish2$GWprod/trophicfish2$N
trophicfish2$IR <- trophicfish2$IRprod/trophicfish2$N
trophicfish2$Mpsgut <- trophicfish2$Mpsgutprod/trophicfish2$N
trophicfish2$SDMPsgut <- trophicfish2$SDMPsgutprod/trophicfish2$N


summary(trophicfish2)
length(trophicfish2$species) # 873 data points
length(trophicfish$species) # consolidated from 977 data point
length(unique(trophicfish2$species)) # ~647 species
length(unique(trophicfish2$family)) # from 171 families
sum(na.omit(trophicfish2$N))  # 30,543 individuals  

trophicfish2$study <- with(trophicfish2, paste0(author, year))
trophicfish2$study <- as.factor(trophicfish2$study)
head(trophicfish2)
length(unique(trophicfish2$study))  # 124 + 1 studies

summary(trophicfish2$min.size)

hist(trophicfish2$min.size)

trophicfish2$polymer.ID <- 
  mapvalues(trophicfish2$polymer.ID, 
            from = c('yes','no'),
            to = c('Polymer ID used',
                   'Polymer ID not used'))
trophicfish2$blanks <- 
  mapvalues(trophicfish2$blanks,
            from = c('yes', 'no'),
            to = c('Blanks used',
                   'Blanks not used'))

trophicfish2$exclude.fib <- 
  mapvalues(trophicfish2$exclude.fib,
            from = c('yes', 'no'),
            to = c('Fibres excluded',
                   'Fibres not excluded'))

trophicfish2$feeding.habit <- as.factor(trophicfish2$feeding.habit)

levels(trophicfish2$feeding.habit)

trophicfish2$feeding.habit <- mapvalues(trophicfish2$feeding.habit, 
                                        from = levels(trophicfish2$feeding.habit),
                                        to = c("Not listed", 
                                               "Browsing on substrate",
                                               "Detritus",
                                               "Filtering plankton",
                                               "Grazing on aquatic plants", 
                                               "Hunting macrofauna",
                                               "Other",
                                               "Selective plankton feeding", 
                                               "Variable"))

levels(trophicfish2$feeding.habit)

trophicfish2$feeding.habit <- factor(
  trophicfish2$feeding.habit,
  levels = c(
    "Not listed",
    "Detritus",
    "Browsing on substrate",
    "Grazing on aquatic plants",
    "Filtering plankton",
    "Selective plankton feeding",
    "Hunting macrofauna",
    "Variable",
    "Other"
  )
)

levels(trophicfish2$feeding.habit)

trophicfish2$study.habitat <- as.factor(trophicfish2$study.habitat)

trophicfish2$environment <- as.factor(trophicfish2$environment)

trophicfish2$region <- as.factor(trophicfish2$region)


#### TL gut mod - set up the data ####

gutdata <- subset(trophicfish2, Mpsgut != 'NA' & 
                    environment != '')
summary(gutdata)
summary(gutdata$author)
nrow(gutdata) # 735 data points
length(unique(gutdata$species)) # 552 species
length(unique(gutdata$family)) # from 157 families
length(unique(gutdata$study)) # from 106 studies

summary(gutdata)
gutdata$region <- as.character(gutdata$region)
gutdata$region <- as.factor(gutdata$region)
summary(gutdata$region)

gutdata$environment <- as.character(gutdata$environment)
gutdata$environment <- as.factor(gutdata$environment)
summary(gutdata$environment)

summary(gutdata$TL)
gutdata$TL <- as.numeric(gutdata$TL)

summary(gutdata$Mpsgut)
gutdata$Mpsgut <- as.numeric(gutdata$Mpsgut)

gutdata$family <- as.factor(gutdata$family)

gutdata <- subset(gutdata, !is.na(TL) & !is.na(Mpsgut))  # remove NAs

gutdata$polymer.ID <- as.factor(gutdata$polymer.ID)

gutdata$blanks <- as.factor(gutdata$blanks)

gutdata$exclude.fib <- as.factor(gutdata$exclude.fib)

gutdata$region <- as.character(gutdata$region)
gutdata$region <- as.factor(gutdata$region)

## Convert to total MP counts for a species for a study

gutdata$totalcount <- round(with(gutdata, Mpsgut * N), 0)

gutdata$zeros <- ifelse(gutdata$totalcount == 0, 0, 1)

#### TL gut mod - fit model ####

TLgutmod <- function() {
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dnegbin(p[i], r)
    p[i] <- r/(r + mu[i])
    log(mu[i]) <-
      log(sample.size[i]) +
      alpha_region[region[i]] +
      beta_min.size*min.size[i] +
      beta_TL[region[i]]*TL[i] +
      gamma_PID[PID[i]] +
      gamma_exclude.fibs[exclude.fibs[i]] +
      gamma_blanks[blanks[i]]
    
    ## Fitted values
    fitted_counts[i] ~ dnegbin(p[i], r)
    fitted[i] <- fitted_counts[i] / sample.size[i]
  }
  
  ## Priors
  r ~ dexp(1)
  
  for(j in 1:nregion) {
    alpha_region[j] <- B[j,1]
    beta_TL[j] <- B[j,2]
    B[j,1:2] ~ dmnorm(B.hat[j,], Tau.B[,])
    B.hat[j,1] <- mu_region
    B.hat[j,2] <- mu_TL
  }
  mu_region ~ dnorm(0, 1)
  mu_TL ~ dnorm(0, 1)
  
  Tau.B[1:2, 1:2] <- inverse(Sigma.B[,])
  Sigma.B[1,1] <- pow(sigma_region, 2)
  sigma_region ~ dexp(1)
  Sigma.B[2,2] <- pow(sigma_TL, 2)
  sigma_TL ~ dexp(1)
  Sigma.B[1,2] <- rho * sigma_region * sigma_TL
  Sigma.B[2,1] <- Sigma.B[1,2]
  rho ~ dnorm(0, 0.5); T(-1, 1)
  
  beta_min.size ~ dnorm(-1, 1)
  
  for(l in 1:2) {
    gamma_PID[l] ~ dnorm(0, 1)
    gamma_exclude.fibs[l] ~ dnorm(0, 1)
    gamma_blanks[l] ~ dnorm(0, 1)
  }
}
  
## Generate initial values for MCMC

TLgutinit <- function()
{
  list(
    "r" = 83.6,
    "mu_region" = rnorm(1),
    "sigma_region" = 0.4,
    "mu_TL" = rnorm(1),
    "sigma_TL" = 0.1,
    "rho" = -0.35,
    "beta_min.size" = rnorm(1, -.15, 0.1),
    "gamma_PID" = rnorm(2),
    "gamma_exclude.fibs" = rnorm(2),
    "gamma_blanks" = rnorm(2)
  )
}

## Keep track of parameters

TLgutparam <- c("r", "alpha_region", 
                "beta_min.size", "beta_TL", 
                "mu_TL", "gamma_PID", "gamma_exclude.fibs", "gamma_blanks",
                "rho")

## Specify data

TLgutdata <-
  list(
    y = gutdata$totalcount,
    N = nrow(gutdata),
    sample.size = as.numeric(gutdata$N),
    region = as.integer(gutdata$region),
    TL = as.numeric(scale(gutdata$TL, center = TRUE)),
    min.size = as.numeric(scale(gutdata$min.size, center = TRUE)),
    PID = as.integer(gutdata$polymer.ID),
    exclude.fibs = as.integer(gutdata$exclude.fib),
    blanks = as.integer(gutdata$blanks),
    nregion = max(as.integer(gutdata$region))
  )

## Run the model
run1 <- jags.parallel(
  data = TLgutdata,
  inits = TLgutinit,
  parameters.to.save = TLgutparam,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 75000,
  n.burnin = 5000,
  n.thin = 25,
  jags.seed = 3156,
  model = TLgutmod
)

run1
run1mcmc <- as.mcmc(run1)
xyplot(run1mcmc, layout = c(6, ceiling(nvar(run1mcmc)/6)))

#### TL gut mod - diagnostics ####

TLgutparam2 <- c("fitted", "fitted_counts")

run2 <- jags.parallel(
  data = TLgutdata,
  inits = TLgutinit,
  parameters.to.save = TLgutparam2,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 75000,
  n.burnin = 5000,
  n.thin = 25,
  jags.seed = 6546,
  model = TLgutmod
)

TLgutmod.response <- t(run2$BUGSoutput$sims.list$fitted_counts)
TLgutmod.observed <- gutdata$totalcount
TLgutmod.fitted <- apply(t(run2$BUGSoutput$sims.list$fitted_counts),
                         1,
                         median)

check.TLgutmod <- createDHARMa(simulatedResponse = TLgutmod.response,
                               observedResponse = TLgutmod.observed, 
                               fittedPredictedResponse = TLgutmod.fitted,
                               integerResponse = T)

plot(check.TLgutmod)

#### TL gut mod - inference ####

gutmod_paramnames <-
  c(
    "Africa - Inland Waters (FAO area)",
    "Europe - Inland Waters (FAO area)",
    "Indian Ocean, Antarctic (FAO area)",
    "Indian Ocean, Eastern (FAO area)",
    "Indian Ocean, Western (FAO area)",
    "Mediterranean and Black Sea (FAO area)",
    "Pacific, Eastern Central (FAO area)",
    "Pacific, Northeast (FAO area)",
    "Pacific, Northwest (FAO area)",
    "Pacific, Southeast (FAO area)",
    "Pacific, Southwest (FAO area)",
    "America, North - Inland Waters (FAO area)",
    "Pacific, Western Central (FAO area)",
    "America, South - Inland Waters (FAO area)",
    "Asia - Inland Waters (FAO area)",
    "Atlantic, Eastern Central (FAO area)",
    "Atlantic, Northeast (FAO area)",
    "Atlantic: Northwest (FAO area)",
    "Atlantic, Southwest (FAO area)",
    "Atlantic, Western Central (FAO area)",
    "Standardized lowest detectable particle size (microns)",
    "Standardized trophic level:Africa - Inland Waters",
    "Standardized trophic level:Europe - Inland Waters",
    "Standardized trophic level:Indian Ocean, Antarctic",
    "Standardized trophic level:Indian Ocean, Eastern",
    "Standardized trophic level:Indian Ocean, Western",
    "Standardized trophic level:Mediterranean and Black Sea",
    "Standardized trophic level:Pacific, Eastern Central",
    "Standardized trophic level:Pacific, Northeast",
    "Standardized trophic level:Pacific, Northwest",
    "Standardized trophic level:Pacific, Southeast",
    "Standardized trophic level:Pacific, Southwest",
    "Standardized trophic level:America, North - Inland Waters",
    "Standardized trophic level:Pacific, Western Central",
    "Standardized trophic level:America, South - Inland Waters",
    "Standardized trophic level:Asia - Inland Waters",
    "Standardized trophic level:Atlantic, Eastern Central",
    "Standardized trophic level:Atlantic, Northeast",
    "Standardized trophic level:Atlantic, Northwest",
    "Standardized trophic level:Atlantic, Southwest",
    "Standardized trophic level:Atlantic, Western Central",
    "Blanks not used",
    "Blanks used",
    "Fibres excluded",
    "Fibres not excluded",
    "Polymer ID not used",
    "Polymer ID used")

## Posterior density plots

run1long <- extract.post(run1)

run1long$variable <- mapvalues(run1long$variable,
                               from = levels(run1long$variable),
                               to = gutmod_paramnames)
  
run1long$order <- c(nrow(run1long):1)

png(
  'Gut Content Posteriors Plot.png',
  width = 14,
  height = 15,
  units = 'cm',
  res = 500
)

ggplot(run1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[1],
    colour = pal[5],
    alpha = 0.75
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.5,
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-2.5, 3)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()

gutdata$post.predict <-
  apply(TLgutmod.response, 1, median)/gutdata$N
gutdata$lower95 <-
  apply(TLgutmod.response, 1, quantile, prob = 0.025)/gutdata$N
gutdata$upper95 <-
  apply(TLgutmod.response, 1, quantile, prob = 0.975)/gutdata$N

freshplot1 <-
  ggplot(subset(gutdata, study.habitat == 'Freshwater')) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95), 
              alpha = 0.75, fill = pal[3]) +
  geom_line(aes(x = TL, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = TL, y = Mpsgut),
             shape = 1, size = 1) +
  facet_wrap(~ region, ncol = 5,
             labeller = label_wrap_gen(width = 15)) +
  labs(x = 'Trophic Level',
       y = expression(paste(
         '          Microplastic\nConcentration (particles ' ~ ind ^ -1 * ')'
       )),
       size = 'Sample Size') +
  scale_fill_manual(values = pal[c(1, 4)],
                    name = 'Fibres Excluded?') +
  scale_colour_manual(values = pal[c(1, 4)],
                      name = 'Fibres Excluded?') +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 10, 30)) +
  theme1

marineplot1 <-
  ggplot(subset(gutdata, study.habitat == 'Marine')) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95), 
              alpha = 0.75, fill = pal[2]) +
  geom_line(aes(x = TL, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = TL, y = Mpsgut),
             shape = 1, size = 1) +
  facet_wrap(~ region, ncol = 5,
             labeller = label_wrap_gen(width = 15)) +
  labs(x = 'Trophic Level',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       )),
       size = 'Sample Size') +
  scale_fill_manual(values = pal[c(1, 4)],
                    name = 'Fibres Excluded?') +
  scale_colour_manual(values = pal[c(1, 4)],
                      name = 'Fibres Excluded?') +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 10, 30)) +
  theme1


png('Gut Content Bayesian Plot.png', width = 14, height = 16, units = 'cm', 
    res = 500)

plot_grid(freshplot1, marineplot1, labels = c('A', 'B'), rel_heights = c(1,2.9),
          nrow = 2, align = 'v')

dev.off()

## Plot according to min size

freshplot2 <-
  ggplot(subset(gutdata, study.habitat == 'Freshwater')) +
  geom_ribbon(aes(x = min.size, ymin = lower95, ymax = upper95), 
              alpha = 0.75, fill = pal[3]) +
  geom_line(aes(x = min.size, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = min.size, y = Mpsgut),
             shape = 1, size = 1) +
  facet_wrap(~ region, ncol = 5,
             labeller = label_wrap_gen(width = 15)) +
  labs(x = expression(paste('Minimum Detectable Particle Size ('*mu*'m)')),
       y = expression(paste(
         '          Microplastic\nConcentration (particles ' ~ ind ^ -1 * ')'
       )),
       size = 'Sample Size') +
  scale_y_continuous(trans = 'log1p', 
                     breaks = c(0, 1, 10, 30),
                     expand = c(0,0)) +
  scale_x_continuous(trans= 'log1p', breaks = c(0, 1, 10, 100, 500)) +
  theme1

marineplot2 <-
  ggplot(subset(gutdata, study.habitat == 'Marine')) +
  geom_ribbon(aes(x = min.size, ymin = lower95, ymax = upper95), 
              alpha = 0.75, fill = pal[2]) +
  geom_line(aes(x = min.size, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = min.size, y = Mpsgut),
             shape = 1, size = 1) +
  facet_wrap(~ region, ncol = 5,
             labeller = label_wrap_gen(width = 15)) +
  labs(x = expression(paste('Minimum Detectable Particle Size ('*mu*'m)')),
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       )),
       size = 'Sample Size') +
  scale_y_continuous(trans = 'log1p', 
                     breaks = c(0, 1, 10, 30),
                     expand = c(0,0)) +
  scale_x_continuous(trans= 'log1p', breaks = c(0, 1, 10, 100, 500)) +
  theme1

png(
  'Gut Content Bayesian Plot According to Lower Size Limit.png',
  width = 14,
  height = 16,
  units = 'cm',
  res = 500
)

plot_grid(freshplot2, marineplot2, labels = c('A', 'B'), rel_heights = c(1,2.5),
          nrow = 2, align = 'v')

dev.off()

#### TL gut mod - predictions ####

## Simulate results according to trophic level

set.seed(541654)

gutdata.sim1 <-
  data.frame(
    TL = seq(
      from = 1.9,
      to = 5.1,
      length.out = 7000
    ),
    min.size = rep(100, 7000),
    sample.size = rpois(7000, 38.7),
    PID = as.integer(rep(2, 7000)),
    blanks = as.integer(rep(2, 7000)),
    fibres = as.integer(rep(2, 7000)),
    region = as.factor(sample(seq(
      from = 1,
      to = max(as.integer(gutdata$region))
    ),
    7000,
    replace = TRUE))
  )

gutdata.sim1$stand.TL <-
  (gutdata.sim1$TL - mean(gutdata$TL)) /
  sd(gutdata$TL - mean(gutdata$TL))

gutdata.sim1$stand.min.size <-
  (gutdata.sim1$min.size - mean(gutdata$min.size)) /
  sd(gutdata$min.size - mean(gutdata$min.size))

for (i in 1:7000) {
  mu <-
    exp(
      log(gutdata.sim1$sample.size[i]) +
        run1$BUGSoutput$sims.list$alpha_region[, gutdata.sim1$region[i]] +
        run1$BUGSoutput$sims.list$beta_min.size * gutdata.sim1$stand.min.size[i] +
        run1$BUGSoutput$sims.list$beta_TL[, gutdata.sim1$region[i]] *
        gutdata.sim1$stand.TL[i] +
        run1$BUGSoutput$sims.list$gamma_PID[, gutdata.sim1$PID[i]] +
        run1$BUGSoutput$sims.list$gamma_exclude.fibs[, gutdata.sim1$fibres[i]] +
        run1$BUGSoutput$sims.list$gamma_blanks[, gutdata.sim1$blanks[i]]
    )
  r <- run1$BUGSoutput$sims.list$r
  p <- r / (mu + r)
  y <- rnbinom(p,
               prob = p,
               size = r) / gutdata.sim1$sample.size[i]
  gutdata.sim1$median[i] <- median(mu / gutdata.sim1$sample.size[i])
  gutdata.sim1$upper25[i] <- quantile(y, 0.625)
  gutdata.sim1$lower25[i] <- quantile(y, 0.375)
  gutdata.sim1$upper50[i] <- quantile(y, 0.75)
  gutdata.sim1$lower50[i] <- quantile(y, 0.25)
  gutdata.sim1$upper75[i] <- quantile(y, 0.875)
  gutdata.sim1$lower75[i] <- quantile(y, 0.125)
  gutdata.sim1$upper95[i] <- quantile(y, 0.975)
  gutdata.sim1$lower95[i] <- quantile(y, 0.025)
  gutdata.sim1$sample[i] <- sample(y, 1)
}

summary(gutdata.sim1)

gutdata.sim1$region <- mapvalues(gutdata.sim1$region,
                                 from = levels(gutdata.sim1$region),
                                 to = c("Africa - Inland Waters",
                                        "America: North - Inland Waters",
                                        "America: South - Inland Waters",
                                        "Asia - Inland Waters",
                                        "Atlantic: Eastern Central",
                                        "Atlantic: Northeast",
                                        "Atlantic: Northwest",
                                        "Atlantic: Southwest",
                                        "Atlantic: Western Central",
                                        "Europe - Inland Waters",
                                        "Indian Ocean: Antarctic",
                                        "Indian Ocean: Eastern",
                                        "Indian Ocean: Western",
                                        "Mediterranean and Black Sea",
                                        "Pacific: Eastern Central",
                                        "Pacific: Northeast",
                                        "Pacific: Northwest",
                                        "Pacific: Southeast",
                                        "Pacific: Southwest",
                                        "Pacific: Western Central"))

gutdata.sim1$area <-
  with(gutdata.sim1,
       as.factor(
         ifelse(
           region != "Africa - Inland Waters" &
             region != "America: North - Inland Waters" &
             region != "America: South - Inland Waters" &
             region != "Asia - Inland Waters" &
             region != "Europe - Inland Waters",
           "Marine",
           "Freshwater"
         )
       ))

gutdata$area <-
  as.factor(with(
    gutdata,
    ifelse(
      region != "Africa - Inland Waters" &
        region != "America: North - Inland Waters" &
        region != "America: South - Inland Waters" &
        region != "Asia - Inland Waters" &
        region != "Europe - Inland Waters",
      "Marine",
      "Freshwater"
    )
  ))

gutdata.a <-
  ggplot() +
  geom_ribbon(
    data = subset(gutdata.sim1, area == "Freshwater"),
    aes(x = TL,
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(gutdata.sim1, area == "Freshwater"),
    aes(x = TL,
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(gutdata.sim1, area == "Freshwater"),
    aes(x = TL,
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(gutdata.sim1, area == "Freshwater"),
    aes(x = TL,
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(data = subset(gutdata.sim1, area == "Freshwater"),
            aes(x = TL,
                y = median),
            colour = pal[5]) +
  geom_point(
    data = subset(gutdata, area == "Freshwater"),
    aes(x = TL,
        y = Mpsgut),
    shape = 1,
    size = 0.5,
    colour = pal[5]
  ) +
  facet_wrap( ~ region, ncol = 5,
              labeller = label_wrap_gen(width = 15)) +
  labs(x = 'Trophic Level',
       y = "") +
  coord_cartesian(ylim = c(0, 150)) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 10, 100),
                     expand = c(0, 0.05)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(1.9, 5.1)) +
  theme1

gutdata.b <-
  ggplot() +
  geom_ribbon(
    data = subset(gutdata.sim1, area == "Marine"),
    aes(x = TL,
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(gutdata.sim1, area == "Marine"),
    aes(x = TL,
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(gutdata.sim1, area == "Marine"),
    aes(x = TL,
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(gutdata.sim1, area == "Marine"),
    aes(x = TL,
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(data = subset(gutdata.sim1, area == "Marine"),
            aes(x = TL,
                y = median),
            colour = pal[5]) +
  geom_point(
    data = subset(gutdata, area == "Marine"),
    aes(x = TL,
        y = Mpsgut),
    shape = 1,
    size = 0.5,
    colour = pal[5]
  ) +
  facet_wrap( ~ region, ncol = 5,
              labeller = label_wrap_gen(width = 15)) +
  labs(x = 'Trophic Level',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       ))) +
  coord_cartesian(ylim = c(0, 150)) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 10, 100),
                     expand = c(0, 0.05)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(1.9, 5.1)) +
  theme1

png('MPs by Trophic Level Predictions Plot.png', width = 14, height = 13, 
    units = 'cm', res = 500)

plot_grid(
  gutdata.a,
  gutdata.b,
  labels = c('A', 'B'),
  rel_heights = c(1, 2.5),
  nrow = 2,
  align = 'v'
)

dev.off()

## Simulate results according to methodology

set.seed(63189)

gutdata.sim2 <-
  data.frame(TL = rep(3, 7000),
             min.size = seq(from = 0.5, to = 520, length.out = 7000),
             sample.size = rpois(7000, 38.7),
             PID = as.factor(rep(2, 7000)),
             blanks = as.factor(rep(2, 7000)),
             fibres = as.factor(sample(as.integer(gutdata$exclude.fib), 
                                            7000, 
                                            replace = TRUE)),
             region = as.factor(rep(16, 7000)))

gutdata.sim2$stand.TL <-
  (gutdata.sim2$TL - mean(gutdata$TL)) /
  sd(gutdata$TL - mean(gutdata$TL))

gutdata.sim2$stand.min.size <-
  (gutdata.sim2$min.size - mean(gutdata$min.size)) /
  sd(gutdata$min.size - mean(gutdata$min.size))

for (i in 1:7000) {
  mu <-
    exp(
      log(gutdata.sim2$sample.size[i]) +
        run1$BUGSoutput$sims.list$alpha_region[, gutdata.sim2$region[i]] +
        run1$BUGSoutput$sims.list$beta_min.size * gutdata.sim2$stand.min.size[i] +
        run1$BUGSoutput$sims.list$beta_TL[, gutdata.sim2$region[i]] *
        gutdata.sim1$stand.TL[i] +
        run1$BUGSoutput$sims.list$gamma_PID[, gutdata.sim2$PID[i]] +
        run1$BUGSoutput$sims.list$gamma_exclude.fibs[, gutdata.sim2$fibres[i]] +
        run1$BUGSoutput$sims.list$gamma_blanks[, gutdata.sim2$blanks[i]]
    )
  r <- run1$BUGSoutput$sims.list$r
  p <- r / (mu + r)
  y <- rnbinom(p,
               prob = p,
               size = r) / gutdata.sim2$sample.size[i]
  gutdata.sim2$median[i] <- median(mu / gutdata.sim2$sample.size[i])
  gutdata.sim2$upper25[i] <- quantile(y, 0.625)
  gutdata.sim2$lower25[i] <- quantile(y, 0.375)
  gutdata.sim2$upper50[i] <- quantile(y, 0.75)
  gutdata.sim2$lower50[i] <- quantile(y, 0.25)
  gutdata.sim2$upper75[i] <- quantile(y, 0.875)
  gutdata.sim2$lower75[i] <- quantile(y, 0.125)
  gutdata.sim2$upper95[i] <- quantile(y, 0.975)
  gutdata.sim2$lower95[i] <- quantile(y, 0.025)
}

gutdata.sim2$exclude.fib <- mapvalues(gutdata.sim2$fibres,
                                      from = levels(gutdata.sim2$fibres),
                                      to = c("Fibres excluded",
                                             "Fibres not excluded"))

png('MPs Methodology Predictions Plot.png', width = 9, height = 5, 
    units = 'cm', res = 500)

ggplot() +
  geom_ribbon(
    data = gutdata.sim2,
    aes(x = min.size,
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = gutdata.sim2,
    aes(x = min.size,
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = gutdata.sim2,
    aes(x = min.size,
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = gutdata.sim2,
    aes(x = min.size,
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(data = gutdata.sim2,
            aes(x = min.size,
                y = median),
            colour = pal[5]) +
  geom_point(
    data = gutdata,
    aes(x = min.size,
        y = Mpsgut),
    shape = 1,
    size = 0.5,
    colour = pal[5]
  ) +
  facet_grid(. ~ exclude.fib, labeller = label_wrap_gen(width = 15)) +
  coord_cartesian(ylim = c(0, 50), xlim = c(0, 510)) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 10, 50),
                     expand = c(0, 0.05)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = expression(paste('Lowest Detectable Particle Size ('*mu*'m)')),
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       ))) +
  theme1

dev.off()


#### Occurrence mod - set up the data ####

ingestion <- subset(trophicfish2, !is.na(IR))

## Convert to successes/failures

ingestion$successes <- round(with(ingestion, N * IR), digits = 0)
ingestion$failures <- round(with(ingestion,  N * (1 - IR), digits = 0))

length(ingestion$species) # 642 data points
length(unique(ingestion$species)) # 478 species
length(unique(ingestion$family)) # from 157 families
length(unique(ingestion$study)) # 108 studies

ingestion$region <- as.character(ingestion$region)
ingestion$region <- as.factor(ingestion$region)

summary(ingestion)

#### Occurrence mod - fit model ####

ingmod <- function()
{
  # Likelihood
  for(i in 1:N)
  {
    y[i] ~ dbinom(p[i], n[i])
    logit(p[i]) <- alpha_region[region[i]] + 
      beta_TL[region[i]]*TL[i] +
      beta_min.size*min.size[i] +
      gamma_exclude.fibs[exclude.fibs[i]] +
      gamma_study[study[i]]
    
    fitted[i] ~ dbinom(p[i], n[i])
  }
  
  # Priors
  
  for(j in 1:nregion) {
    alpha_region[j] <- B[j,1]
    beta_TL[j] <- B[j,2]
    B[j,1:2] ~ dmnorm(B.hat[j,], Tau.B[,])
    B.hat[j,1] <- mu_region
    B.hat[j,2] <- mu_TL
  }
  mu_region ~ dnorm(0, 1)
  mu_TL ~ dnorm(0, 1)
  
  Tau.B[1:2, 1:2] <- inverse(Sigma.B[,])
  Sigma.B[1,1] <- pow(sigma_region, 2)
  sigma_region ~ dexp(1)
  Sigma.B[2,2] <- pow(sigma_TL, 2)
  sigma_TL ~ dexp(1)
  Sigma.B[1,2] <- rho * sigma_region * sigma_TL
  Sigma.B[2,1] <- Sigma.B[1,2]
  rho ~ dnorm(0, 0.5); T(-1, 1)
  
  beta_min.size ~ dnorm(0, 1)
  
  for(k in 1:2) {
    gamma_exclude.fibs[k] ~ dnorm(0, 1)
  }
  
  for(l in 1:nstudy){
    gamma_study[l] ~ dnorm(0, tau_study)
  }
  tau_study <- pow(sigma_study, -2)
  sigma_study ~ dexp(1)
  
}


## Initial values for MCMC chains
inginit <- function()
{
  list(
    "mu_region" = rnorm(1),
    "sigma_region" = rexp(1),
    "mu_TL" = rnorm(1),
    "sigma_TL" = rexp(1),
    "beta_min.size" = rnorm(1),
    "gamma_exclude.fibs" = rnorm(2),
    "sigma_study" = rexp(1)
  )
}

## Parameters to keep track of
ingparam <- c("alpha_region", "beta_TL", "beta_min.size", "gamma_exclude.fibs",
              "sigma_study")


## Specify data
ingdata <- list(
  N = nrow(ingestion),
  nregion = max(as.integer(ingestion$region)),
  y = as.numeric(ingestion$successes),
  n = as.numeric(ingestion$N),
  TL = as.numeric(scale(ingestion$TL, center = TRUE)),
  region = as.integer(ingestion$region),
  min.size = as.numeric(scale(ingestion$min.size, center = TRUE)),
  exclude.fibs = as.integer(as.factor(ingestion$exclude.fib)),
  study = as.integer(as.factor(ingestion$study)),
  nstudy = max(as.integer(as.factor(ingestion$study)))
)


## Run the model
ingrun1 <- jags.parallel(
  data = ingdata,
  inits = inginit,
  parameters.to.save = ingparam,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 7000,
  n.burnin = 1000,
  n.thin = 4,
  jags.seed = 6174,
  model = ingmod
)

ingrun1
ingrun1mcmc <- as.mcmc(ingrun1)
xyplot(ingrun1mcmc, layout = c(6, ceiling(nvar(ingrun1mcmc)/6)))


#### Occurrence mod - diagnostics ####

## Rerun model to extract fitted values
ingparam2 <- c("fitted")

ingrun2 <- jags.parallel(
  data = ingdata,
  inits = inginit,
  parameters.to.save = ingparam2,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 7000,
  n.burnin = 1000,
  n.thin = 4,
  jags.seed = 6174,
  model = ingmod
)

ing.response <- t(ingrun2$BUGSoutput$sims.list$fitted)
ing.observed <- ingestion$successes
ing.fitted <- apply(t(ingrun2$BUGSoutput$sims.list$fitted),
                         1,
                         median)

check.ingmod <- createDHARMa(simulatedResponse = ing.response,
                             observedResponse = ing.observed, 
                             fittedPredictedResponse = ing.fitted,
                             integerResponse = T,
                             seed = 5151)

plot(check.ingmod)  


#### Occurrence mod - inference ####

ingrunparams <-
  c(
    "Africa - Inland Waters",
    "Indian Ocean, Antarctic",
    "Indian Ocean, Eastern",
    "Indian Ocean, Western",
    "Mediterranean and Black Sea",
    "Pacific, Eastern Central",
    "Pacific, Northeast",
    "Pacific, Northwest",
    "Pacific, Southeast",
    "Pacific, Southwest",
    "Pacific, Western Central",
    "America, North - Inland Waters",
    "Asia - Inland Waters",
    "Atlantic, Eastern Central",
    "Atlantic, Northeast",
    "Atlantic, Northwest",
    "Atlantic, Southwest",
    "Atlantic, Western Central",
    "Europe - Inland Waters",
    "Standardized lowest detectable particle size (microns)",
    "Standardized trophic level:Africa - Inland Waters",
    "Standardized trophic level:Indian Ocean, Antarctic",
    "Standardized trophic level:Indian Ocean, Eastern",
    "Standardized trophic level:Indian Ocean, Western",
    "Standardized trophic level:Mediterranean and Black Sea",
    "Standardized trophic level:Pacific, Eastern Central",
    "Standardized trophic level:Pacific, Northeast",
    "Standardized trophic level:Pacific, Northwest",
    "Standardized trophic level:Pacific, Southeast",
    "Standardized trophic level:Pacific, Southwest",
    "Standardized trophic level:Pacific, Western Central",
    "Standardized trophic level:America, North - Inland Waters",
    "Standardized trophic level:Asia - Inland Waters",
    "Standardized trophic level:Atlantic, Eastern Central",
    "Standardized trophic level:Atlantic, Northeast",
    "Standardized trophic level:Atlantic, Northwest",
    "Standardized trophic level:Atlantic, Southwest",
    "Standardized trophic level:Atlantic, Western Central",
    "Standardized trophic level:Europe - Inland Waters",
    "Fibres excluded",
    "Fibres not excluded",
    "Study standard deviation"
  )

ingrun1long <- extract.post(ingrun1)

ingrun1long$variable <- mapvalues(ingrun1long$variable,
                               from = levels(ingrun1long$variable),
                               to = ingrunparams)

ingrun1long$order <- c(nrow(ingrun1long):1)

png(
  'Ingestion Posteriors Plot.png',
  width = 14,
  height = 13,
  units = 'cm',
  res = 500
)

ggplot(ingrun1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[1],
    colour = pal[5],
    alpha = 0.75
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.5,
    colour = pal[3]
  ) +
  labs(x = "",
       y = "Parameter") +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  theme1

dev.off()



ingestion$post.predict <-
  apply(ing.response, 1, median)/ingestion$N
ingestion$lower95 <-
  apply(ing.response, 1, quantile, prob = 0.025)/ingestion$N
ingestion$upper95 <-
  apply(ing.response, 1, quantile, prob = 0.975)/ingestion$N

freshplot3 <-
  ggplot(subset(ingestion, study.habitat == 'Freshwater')) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95),
              fill = pal[3],
              alpha = 0.75) +
  geom_line(aes(x = TL, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = TL, y = IR, size = N),
             shape = 1) +
  facet_wrap( ~ region, ncol = 4,
              labeller = label_wrap_gen(width = 20)) +
  labs(x = 'Trophic Level',
       y = 'Microplastic Occurrence Rate',
       size = 'Sample Size') +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25)) +
  theme1

marineplot3 <-
  ggplot(subset(ingestion, study.habitat == 'Marine')) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95),
              fill = pal[2],
              alpha = 0.75) +
  geom_line(aes(x = TL, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = TL, y = IR, size = N),
             shape = 1) +
  facet_wrap( ~ region, ncol = 4,
              labeller = label_wrap_gen(width = 20)) +
  labs(x = 'Trophic Level',
       y = 'Microplastic Occurrence Rate',
       size = 'Sample Size') +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25)) +
  theme1


png(
  'Ingestion Rate Bayesian Plot.png',
  width = 14,
  height = 18,
  units = 'cm',
  res = 500
)

plot_grid(
  freshplot3,
  marineplot3,
  labels = c('A', 'B'),
  rel_heights = c(1, 3.2),
  nrow = 2,
  align = 'v'
)

dev.off()

#### Occurrence mod - predictions ####

## Simulate results according to trophic level

set.seed(32156)

ingestion.sim1 <-
  data.frame(
    TL = seq(
      from = 1.9,
      to = 5.1,
      length.out = 10000
    ),
    min.size = rep(100, 10000),
    sample.size = rpois(10000, 38.7),
    fibres = as.integer(rep(2, 10000)),
    region = as.factor(sample(seq(
      from = 1,
      to = max(as.integer(ingestion$region))
    ),
    10000,
    replace = TRUE))
  )

ingestion.sim1$stand.TL <-
  (ingestion.sim1$TL - mean(ingestion$TL)) /
  sd(ingestion$TL - mean(ingestion$TL))

ingestion.sim1$stand.min.size <-
  (ingestion.sim1$min.size - mean(ingestion$min.size)) /
  sd(ingestion$min.size - mean(ingestion$min.size))

for (i in 1:10000) {
  p <-
    plogis(
      ingrun1$BUGSoutput$sims.list$alpha_region[, ingestion.sim1$region[i]] +
        ingrun1$BUGSoutput$sims.list$beta_min.size * ingestion.sim1$stand.min.size[i] +
        ingrun1$BUGSoutput$sims.list$beta_TL[, ingestion.sim1$region[i]] *
        ingestion.sim1$stand.TL[i] +
        ingrun1$BUGSoutput$sims.list$gamma_exclude.fibs[, ingestion.sim1$fibres[i]]
    )
  y <- rbinom(p,
              prob = p,
              size = ingestion.sim1$sample.size[i])/ingestion.sim1$sample.size[i]
  ingestion.sim1$median[i] <- median(p)
  ingestion.sim1$upper25[i] <- quantile(y, 0.625)
  ingestion.sim1$lower25[i] <- quantile(y, 0.375)
  ingestion.sim1$upper50[i] <- quantile(y, 0.75)
  ingestion.sim1$lower50[i] <- quantile(y, 0.25)
  ingestion.sim1$upper75[i] <- quantile(y, 0.875)
  ingestion.sim1$lower75[i] <- quantile(y, 0.125)
  ingestion.sim1$upper95[i] <- quantile(y, 0.975)
  ingestion.sim1$lower95[i] <- quantile(y, 0.025)
  ingestion.sim1$sample[i] <- sample(y, 1)
}

summary(ingestion.sim1)

ingestion.sim1$region <- mapvalues(ingestion.sim1$region,
                                 from = levels(ingestion.sim1$region),
                                 to = c("Africa - Inland Waters",
                                        "America: North - Inland Waters",
                                        "Asia - Inland Waters",
                                        "Atlantic: Eastern Central",
                                        "Atlantic: Northeast",
                                        "Atlantic: Northwest",
                                        "Atlantic: Southwest",
                                        "Atlantic: Western Central",
                                        "Europe - Inland Waters",
                                        "Indian Ocean: Antarctic",
                                        "Indian Ocean: Eastern",
                                        "Indian Ocean: Western",
                                        "Mediterranean and Black Sea",
                                        "Pacific: Eastern Central",
                                        "Pacific: Northeast",
                                        "Pacific: Northwest",
                                        "Pacific: Southeast",
                                        "Pacific: Southwest",
                                        "Pacific: Western Central"))

ingestion.sim1$area <-
  with(ingestion.sim1,
       as.factor(
         ifelse(
           region != "Africa - Inland Waters" &
             region != "America: North - Inland Waters" &
             region != "Asia - Inland Waters" &
             region != "Europe - Inland Waters",
           "Marine",
           "Freshwater"
         )
       ))

ingestion$area <-
  as.factor(with(
    ingestion,
    ifelse(
      region != "Africa - Inland Waters" &
        region != "America: North - Inland Waters" &
        region != "Asia - Inland Waters" &
        region != "Europe - Inland Waters",
      "Marine",
      "Freshwater"
    )
  ))

ing.plot.a <- 
  ggplot() +
  geom_ribbon(
    data = subset(ingestion.sim1, area == "Freshwater"),
    aes(x = TL,
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(ingestion.sim1, area == "Freshwater"),
    aes(x = TL,
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(ingestion.sim1, area == "Freshwater"),
    aes(x = TL,
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(ingestion.sim1, area == "Freshwater"),
    aes(x = TL,
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(data = subset(ingestion.sim1, area == "Freshwater"),
            aes(x = TL,
                y = median),
            colour = pal[5]) +
  geom_point(
    data = subset(ingestion, area == "Freshwater"),
    aes(x = TL,
        y = IR),
    shape = 1,
    size = 0.5,
    colour = pal[5]
  ) +
  facet_wrap(~ region,
             labeller = label_wrap_gen(width = 15),
             ncol = 5) +
  labs(x = '',
       y = "") +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
                     expand = c(0, 0.05),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(1.9, 5.1)) +
  theme1

ing.plot.b <- 
  ggplot() +
  geom_ribbon(
    data = subset(ingestion.sim1, area == "Marine"),
    aes(x = TL,
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(ingestion.sim1, area == "Marine"),
    aes(x = TL,
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(ingestion.sim1, area == "Marine"),
    aes(x = TL,
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(ingestion.sim1, area == "Marine"),
    aes(x = TL,
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(data = subset(ingestion.sim1, area == "Marine"),
            aes(x = TL,
                y = median),
            colour = pal[5]) +
  geom_point(
    data = subset(ingestion, area == "Marine"),
    aes(x = TL,
        y = IR),
    shape = 1,
    size = 0.5,
    colour = pal[5]
  ) +
  facet_wrap(~ region,
             labeller = label_wrap_gen(width = 15),
             ncol = 5) +
  labs(x = 'Trophic Level',
       y = "Occurrence Rate") +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
                     expand = c(0, 0.05),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(1.9, 5.1)) +
  theme1

png('Ingestion Rate Predictions Plot.png', width = 14, height = 12, 
    units = 'cm', res = 500)

plot_grid(
  ing.plot.a,
  ing.plot.b,
  labels = c('A', 'B'),
  rel_heights = c(1, 2.6),
  nrow = 2
)

dev.off()

#### Size mod - set up the data ####

size <- subset(gutdata, total.length != 'NA')
size$life.stage <- as.factor(size$life.stage)
size$region <- as.character(size$region)
size$region <- as.factor(size$region)

length(size$total.length)  # 395 data points remaining
length(unique(size$study))  # 62 studies
length(unique(size$species))  # 327 species

size$totalcounts <- round(with(size, Mpsgut*N))

#### Size mod - fit model ####

sizemod <- function()
{
  # Likelihood
  for(i in 1:N) {
    y[i] ~ dnegbin(p[i], r)
    p[i] <- r/(r + mu[i])
    log(mu[i]) <-
      log(sample.size[i]) +
      alpha_region[region[i]] +
      beta_min.size*min.size[i] +
      beta_length[region[i]]*length[i] +
      gamma_exclude.fibs[exclude.fibs[i]]
    
    ## Fitted values
    fitted_counts[i] ~ dnegbin(p[i], r)
    fitted[i] <- fitted_counts[i] / sample.size[i]
  }
  
  ## Priors
  r ~ dexp(1)
  
  for(j in 1:nregion) {
    alpha_region[j] <- B[j,1]
    beta_length[j] <- B[j,2]
    B[j,1:2] ~ dmnorm(B.hat[j,], Tau.B[,])
    B.hat[j,1] <- mu_region
    B.hat[j,2] <- mu_length
  }
  mu_region ~ dnorm(0, 1)
  mu_length ~ dnorm(0, 1)
  
  Tau.B[1:2, 1:2] <- inverse(Sigma.B[,])
  Sigma.B[1,1] <- pow(sigma_region, 2)
  sigma_region ~ dexp(1)
  Sigma.B[2,2] <- pow(sigma_length, 2)
  sigma_length ~ dexp(1)
  Sigma.B[1,2] <- rho * sigma_region * sigma_length
  Sigma.B[2,1] <- Sigma.B[1,2]
  rho ~ dnorm(0, 0.5); T(-1, 1)
  
  beta_min.size ~ dnorm(-1, 1)
  
  for(l in 1:2) {
    gamma_exclude.fibs[l] ~ dnorm(0, 1)
  }
}

## Initial values
sizemod_init <- function()
{
  list(
    "r" = 46.2,
    "mu_region" = rnorm(1),
    "mu_length" = rnorm(1),
    "sigma_region" = 1.1,
    "sigma_length" = 0.1,
    "beta_min.size" = rnorm(1, -0.4, 0.1),
    "gamma_exclude.fibs" = rnorm(2)
  )
}

## Parameters to keep track of
sizemod_params <- c("r", "alpha_region", "beta_length",
                    "beta_min.size", "gamma_exclude.fibs")

## Specify data

sizedata <- list(
  N = nrow(size),
  y = as.numeric(size$totalcounts),
  sample.size = as.numeric(size$N),
  min.size = as.numeric(scale(size$min.size, center = TRUE)),
  length = as.numeric(scale(size$total.length, center = TRUE)),
  region = as.integer(size$region),
  nregion = max(as.integer(size$region)),
  exclude.fibs = as.integer(size$exclude.fib)
)

## Run the model
sizerun1 <- jags.parallel(
  data = sizedata,
  inits = sizemod_init,
  parameters.to.save = sizemod_params,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 50000,
  n.burnin = 5000,
  n.thin = 10,
  jags.seed = 554,
  model = sizemod
)

sizerun1
sizerun1mcmc <- as.mcmc(sizerun1)
xyplot(sizerun1mcmc, layout = c(6, ceiling(nvar(ingrun1mcmc)/6)))

#### Size mod - diagnostics ####

sizemod_params2 <- c("fitted", "fitted_counts")

sizerun2 <- jags.parallel(
  data = sizedata,
  inits = sizemod_init,
  parameters.to.save = sizemod_params2,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 50000,
  n.burnin = 5000,
  n.thin = 10,
  jags.seed = 554,
  model = sizemod
)

sizemod.response <- t(sizerun2$BUGSoutput$sims.list$fitted_counts)
sizemod.observed <- size$totalcounts
sizemod.fitted <- apply(t(sizerun2$BUGSoutput$sims.list$fitted_counts),
                         1,
                         median)

check.sizemod <- createDHARMa(simulatedResponse = sizemod.response,
                               observedResponse = sizemod.observed, 
                               fittedPredictedResponse = sizemod.fitted,
                               integerResponse = T)

plot(check.sizemod)

#### Size mod - inference ####

## Posterior density plots

sizemod_paramnames <-
  c(
    "Africa - Inland Waters (FAO area)",
    "Europe - Inland Waters (FAO area)",
    "Indian Ocean, Antarctic (FAO area)",
    "Indian Ocean, Eastern (FAO area)",
    "Indian Ocean, Western (FAO area)",
    "Mediterranean and Black Sea (FAO area)",
    "Pacific, Eastern Central (FAO area)",
    "Pacific, Northwest (FAO area)",
    "Pacific, Southeast (FAO area)",
    "Pacific, Southwest (FAO area)",
    "Pacific, Western Central (FAO area)",
    "America, North - Inland Waters (FAO area)",
    "America, South - Inland Waters (FAO area)",
    "Asia - Inland Waters (FAO area)",
    "Atlantic, Eastern Central (FAO area)",
    "Atlantic, Northeast (FAO area)",
    "Atlantic, Northwest (FAO area)",
    "Atlantic, Southwest (FAO area)",
    "Atlantic, Western Central (FAO area)",
    "Standardized total length (cm):Africa - Inland Waters",
    "Standardized total length (cm):Europe - Inland Waters",
    "Standardized total length (cm):Indian Ocean, Antarctic",
    "Standardized total length (cm):Indian Ocean, Eastern",
    "Standardized total length (cm):Indian Ocean, Western",
    "Standardized total length (cm):Mediterranean and Black Sea",
    "Standardized total length (cm):Pacific, Eastern Central",
    "Standardized total length (cm):Pacific, Northwest",
    "Standardized total length (cm):Pacific, Southeast",
    "Standardized total length (cm):Pacific, Southwest",
    "Standardized total length (cm):Pacific, Western Central",
    "Standardized total length (cm):America, North - Inland Waters",
    "Standardized total length (cm):America, South - Inland Waters",
    "Standardized total length (cm):Asia - Inland Waters",
    "Standardized total length (cm):Atlantic, Eastern Central",
    "Standardized total length (cm):Atlantic, Northeast",
    "Standardized total length (cm):Atlantic, Northwest",
    "Standardized total length (cm):Atlantic, Southwest",
    "Standardized total length (cm):Atlantic, Western Central",
    "Standardized lowest detectable particle size (microns)",
    "Fibres excluded",
    "Fibres not excluded")

sizerun1long <- extract.post(sizerun1)

sizerun1long$variable <- mapvalues(sizerun1long$variable,
                                  from = levels(sizerun1long$variable),
                                  to = sizemod_paramnames)

sizerun1long$order <- c(nrow(sizerun1long):1)

png(
  'Body Size Model Posteriors Plot.png',
  width = 14,
  height = 12,
  units = 'cm',
  res = 500
)

ggplot(sizerun1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[1],
    colour = pal[5],
    alpha = 0.75
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.5,
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()

size$post.predict <-
  apply(sizemod.response, 1, median)
size$lower95 <-
  apply(sizemod.response, 1, quantile, prob = 0.025)
size$upper95 <-
  apply(sizemod.response, 1, quantile, prob = 0.975)

## Plot
png(
  'Size Gut MPs Bayesian Plot.png',
  width = 9,
  height = 11,
  units = 'cm',
  res = 500
)

ggplot(size) +
  geom_ribbon(
    aes(x = total.length, ymin = lower95, ymax = upper95),
    fill = pal[1],
    alpha = 0.75
  ) +
  geom_line(aes(x = total.length, y = post.predict),
            size = 0.5,
            alpha = 0.8) +
  geom_point(aes(x = total.length, y = Mpsgut),
             shape = 1,
             size = 0.25) +
  facet_wrap(~ min.size) +
  labs(x = 'Total Length (cm)',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       ))) +
  scale_x_continuous(trans = 'log1p', 
                     breaks = c(0, 1, 10, 100)) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 10, 30)) +
  theme1

dev.off()

#### Size mod - predictions ####

set.seed(6554)

size.sim1 <-
  data.frame(
    length = seq(
      from = 1,
      to = 200,
      length.out = 7000
    ),
    min.size = rep(100, 7000),
    sample.size = rpois(7000, 38.7),
    fibres = as.integer(rep(2, 7000)),
    region = as.factor(sample(seq(
      from = 1,
      to = max(as.integer(size$region))
    ),
    7000,
    replace = TRUE))
  )

size.sim1$stand.length <-
  (size.sim1$length - mean(size$total.length)) /
  sd(size$total.length - mean(size$total.length))

size.sim1$stand.min.size <-
  (size.sim1$min.size - mean(size$min.size)) /
  sd(size$min.size - mean(size$min.size))

for(i in 1:7000) {
  mu <-
    exp(
      log(size.sim1$sample.size[i]) +
        sizerun1$BUGSoutput$sims.list$alpha_region[, size.sim1$region[i]] +
        sizerun1$BUGSoutput$sims.list$beta_min.size * size.sim1$stand.min.size[i] +
        sizerun1$BUGSoutput$sims.list$beta_length[, size.sim1$region[i]] *
        size.sim1$stand.length[i] +
        sizerun1$BUGSoutput$sims.list$gamma_exclude.fibs[, size.sim1$fibres[i]]
    )
  r <- sizerun1$BUGSoutput$sims.list$r
  p <- r / (mu + r)
  y <- rnbinom(p,
               prob = p,
               size = r) / size.sim1$sample.size[i]
  size.sim1$median[i] <- median(mu ) / size.sim1$sample.size[i]
  size.sim1$upper25[i] <- quantile(y, 0.625)
  size.sim1$lower25[i] <- quantile(y, 0.375)
  size.sim1$upper50[i] <- quantile(y, 0.75)
  size.sim1$lower50[i] <- quantile(y, 0.25)
  size.sim1$upper75[i] <- quantile(y, 0.875)
  size.sim1$lower75[i] <- quantile(y, 0.125)
  size.sim1$upper95[i] <- quantile(y, 0.975)
  size.sim1$lower95[i] <- quantile(y, 0.025)
}

summary(size.sim1)

size.sim1$region <- mapvalues(size.sim1$region,
                                 from = levels(size.sim1$region),
                                 to = c("Africa - Inland Waters",
                                        "America: North - Inland Waters",
                                        "America: South - Inland Waters",
                                        "Asia - Inland Waters",
                                        "Atlantic: Eastern Central",
                                        "Atlantic: Northeast",
                                        "Atlantic: Northwest",
                                        "Atlantic: Southwest",
                                        "Atlantic: Western Central",
                                        "Europe - Inland Waters",
                                        "Indian Ocean: Antarctic",
                                        "Indian Ocean: Eastern",
                                        "Indian Ocean: Western",
                                        "Mediterranean and Black Sea",
                                        "Pacific: Eastern Central",
                                        "Pacific: Northwest",
                                        "Pacific: Southeast",
                                        "Pacific: Southwest",
                                        "Pacific: Western Central"))

size.sim1$area <-
  with(size.sim1,
       as.factor(
         ifelse(
           region != "Africa - Inland Waters" &
             region != "America: North - Inland Waters" &
             region != "America: South - Inland Waters" &
             region != "Asia - Inland Waters" &
             region != "Europe - Inland Waters",
           "Marine",
           "Freshwater"
         )
       ))

size$area <-
  as.factor(with(
    size,
    ifelse(
      region != "Africa - Inland Waters" &
        region != "America: North - Inland Waters" &
        region != "America: South - Inland Waters" &
        region != "Asia - Inland Waters" &
        region != "Europe - Inland Waters",
      "Marine",
      "Freshwater"
    )
  ))

size.a <-
  ggplot() +
  geom_ribbon(
    data = subset(size.sim1, area == "Freshwater"),
    aes(x = length,
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(size.sim1, area == "Freshwater"),
    aes(x = length,
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(size.sim1, area == "Freshwater"),
    aes(x = length,
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(size.sim1, area == "Freshwater"),
    aes(x = length,
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(data = subset(size.sim1, area == "Freshwater"),
            aes(x = length,
                y = median),
            colour = pal[5]) +
  geom_point(
    data = subset(size, area == "Freshwater"),
    aes(x = total.length,
        y = Mpsgut),
    shape = 1,
    size = 0.5,
    colour = pal[5]
  ) +
  facet_wrap( ~ region, ncol = 5,
              labeller = label_wrap_gen(width = 15)) +
  labs(x = "",
       y = "") +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 10, 100),
                     expand = c(0, 0.1)) +
  scale_x_continuous(expand = c(0, 0.1),
                     limits = c(1, 200),
                     breaks = c(1, 50, 100, 150, 200)) +
  theme1

size.b <-
  ggplot() +
  geom_ribbon(
    data = subset(size.sim1, area == "Marine"),
    aes(x = length,
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(size.sim1, area == "Marine"),
    aes(x = length,
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(size.sim1, area == "Marine"),
    aes(x = length,
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(size.sim1, area == "Marine"),
    aes(x = length,
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(data = subset(size.sim1, area == "Marine"),
            aes(x = length,
                y = median),
            colour = pal[5]) +
  geom_point(
    data = subset(size, area == "Marine"),
    aes(x = total.length,
        y = Mpsgut),
    shape = 1,
    size = 0.5,
    colour = pal[5]
  ) +
  facet_wrap( ~ region, ncol = 5,
              labeller = label_wrap_gen(width = 15)) +
  labs(x = 'Total Length (cm)',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       ))) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 10, 100),
                     expand = c(0, 0.1)) +
  scale_x_continuous(expand = c(0, 0.1),
                     limits = c(1, 200),
                     breaks = c(1, 50, 100, 150, 200)) +
  theme1

png('MPs by Total Length Predictions Plot.png', width = 14, height = 12, 
    units = 'cm', res = 500)

plot_grid(
  size.a,
  size.b,
  labels = c('A', 'B'),
  rel_heights = c(1, 2.5),
  nrow = 2,
  align = 'v'
)

dev.off()

#### Occurrence x size mod - set up the data ####

sizeing <- subset(size, !is.na(IR))

## Convert to successes/failures

sizeing$successes <- round(with(sizeing, N * IR), digits = 0)
sizeing$failures <- round(with(sizeing,  N * (1 - IR), digits = 0))

length(sizeing$species) # 258 data points
length(unique(sizeing$species)) # 215 species
length(unique(sizeing$family)) # from 93 families
length(unique(sizeing$study)) # 48 studies

sizeing$region <- as.character(sizeing$region)
sizeing$region <- as.factor(sizeing$region)

summary(sizeing)

#### Occurrence x size mod - fit model ####

sizeingmod <- function()
{
  # Likelihood
  for(i in 1:N)
  {
    y[i] ~ dbinom(p[i], n[i])
    logit(p[i]) <- alpha_region[region[i]] + 
      beta_length[region[i]]*length[i] +
      beta_min.size*min.size[i] +
      gamma_exclude.fibs[exclude.fibs[i]] +
      gamma_study[study[i]]
    
    fitted[i] ~ dbinom(p[i], n[i])
  }
  
  # Priors
  
  for(j in 1:nregion) {
    alpha_region[j] <- B[j,1]
    beta_length[j] <- B[j,2]
    B[j,1:2] ~ dmnorm(B.hat[j,], Tau.B[,])
    B.hat[j,1] <- mu_region
    B.hat[j,2] <- mu_length
  }
  mu_region ~ dnorm(0, 1)
  mu_length ~ dnorm(0, 1)
  
  Tau.B[1:2, 1:2] <- inverse(Sigma.B[,])
  Sigma.B[1,1] <- pow(sigma_region, 2)
  sigma_region ~ dexp(1)
  Sigma.B[2,2] <- pow(sigma_length, 2)
  sigma_length ~ dexp(1)
  Sigma.B[1,2] <- rho * sigma_region * sigma_length
  Sigma.B[2,1] <- Sigma.B[1,2]
  rho ~ dnorm(0, 0.5); T(-1, 1)
  
  beta_min.size ~ dnorm(0, 1)
  
  for(k in 1:2) {
    gamma_exclude.fibs[k] ~ dnorm(0, 1)
  }
  
  for(l in 1:nstudy){
    gamma_study[l] ~ dnorm(0, tau_study)
  }
  tau_study <- pow(sigma_study, -2)
  sigma_study ~ dexp(1)
  
}

## Initial values for MCMC chains
sizeinginit <- function()
{
  list(
    "mu_region" = rnorm(1),
    "sigma_region" = rexp(1),
    "mu_length" = rnorm(1),
    "sigma_length" = rexp(1),
    "beta_min.size" = rnorm(1),
    "rho" = runif(1, -1, 1),
    "gamma_exclude.fibs" = rnorm(2),
    "sigma_study" = rexp(1)
  )
}

## Parameters to keep track of
sizeingparam <- c("alpha_region", "beta_length", "beta_min.size", 
                  "gamma_exclude.fibs", "sigma_study")


## Specify data
sizeingdata <- list(
  N = nrow(sizeing),
  nregion = max(as.integer(sizeing$region)),
  y = as.numeric(sizeing$successes),
  n = as.numeric(sizeing$N),
  length = as.numeric(scale(sizeing$total.length, center = TRUE)),
  region = as.integer(sizeing$region),
  min.size = as.numeric(scale(sizeing$min.size, center = TRUE)),
  exclude.fibs = as.integer(as.factor(sizeing$exclude.fib)),
  study = as.integer(as.factor(sizeing$study)),
  nstudy = max(as.integer(as.factor(sizeing$study)))
)

## Run the model
sizeingrun1 <- jags.parallel(
  data = sizeingdata,
  inits = sizeinginit,
  parameters.to.save = sizeingparam,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 100000,
  n.burnin = 2000,
  n.thin = 50,
  jags.seed = 123,
  model = sizeingmod
)

sizeingrun1
sizeingrun1mcmc <- as.mcmc(sizeingrun1)
xyplot(sizeingrun1mcmc, layout = c(6, ceiling(nvar(sizeingrun1mcmc)/6)))

#### Occurrence x size mod - diagnostics ####

## Rerun model to extract p
sizeingparam2 <- c("fitted")

sizeingrun2 <- jags.parallel(
  data = sizeingdata,
  inits = sizeinginit,
  parameters.to.save = sizeingparam2,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 100000,
  n.burnin = 2000,
  n.thin = 50,
  jags.seed = 123,
  model = sizeingmod
)

sizeing.response <- t(sizeingrun2$BUGSoutput$sims.list$fitted)
sizeing.observed <- sizeing$successes
sizeing.fitted <- apply(t(sizeingrun2$BUGSoutput$sims.list$fitted),
                    1,
                    median)

check.sizeingmod <- createDHARMa(simulatedResponse = sizeing.response,
                             observedResponse = sizeing.observed, 
                             fittedPredictedResponse = sizeing.fitted,
                             integerResponse = T,
                             seed = 5151)

plot(check.sizeingmod)  

#### Occurrence x size mod - inference ####

sizeing_paramNames <-
  c(
    "Africa - Inland Waters",
    "Indian Ocean, Antarctic",
    "Indian Ocean, Eastern",
    "Indian Ocean, Western",
    "Mediterranean and Black Sea",
    "Pacific, Eastern Central",
    "Pacific, Northwest",
    "Pacific, Southeast",
    "Pacific, Southwest",
    "America, North - Inland Waters",
    "Asia - Inland Waters",
    "Atlantic, Eastern Central",
    "Atlantic, Northeast",
    "Atlantic, Northwest",
    "Atlantic, Southwest",
    "Atlantic, Western Central",
    "Europe - Inland Waters",
    "Standardized total length:Africa - Inland Waters",
    "Standardized total length:Indian Ocean, Antarctic",
    "Standardized total length:Indian Ocean, Eastern",
    "Standardized total length:Indian Ocean, Western",
    "Standardized total length:Mediterranean and Black Sea",
    "Standardized total length:Pacific, Eastern Central",
    "Standardized total length:Pacific, Northwest",
    "Standardized total length:Pacific, Southeast",
    "Standardized total length:Pacific, Southwest",
    "Standardized total length:America, North - Inland Waters",
    "Standardized total length:Asia - Inland Waters",
    "Standardized total length:Atlantic, Eastern Central",
    "Standardized total length:Atlantic, Northeast",
    "Standardized total length:Atlantic, Northwest",
    "Standardized total length:Atlantic, Southwest",
    "Standardized total length:Atlantic, Western Central",
    "Standardized total length:Europe - Inland Waters",
    "Standardized lowest detectable particle size (microns)",
    "Fibres excluded",
    "Fibres not excluded",
    "Study standard deviation"
  )
## Posterior density plots

sizeingrun1long <- extract.post(sizeingrun1)

sizeingrun1long$variable <- mapvalues(sizeingrun1long$variable,
                                   from = levels(sizeingrun1long$variable),
                                   to = sizeing_paramNames)

sizeingrun1long$order <- c(nrow(sizeingrun1long):1)

png(
  'Body Size Ingestion Model Posteriors Plot.png',
  width = 14,
  height = 13,
  units = 'cm',
  res = 500
)

ggplot(sizeingrun1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[1],
    colour = pal[5],
    alpha = 0.75
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.5,
    colour = pal[3]
  ) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()



sizeing$post.predict <-
  apply(sizeing.response, 1, median)/sizeing$N
sizeing$lower95 <-
  apply(sizeing.response, 1, quantile, prob = 0.025)/sizeing$N
sizeing$upper95 <-
  apply(sizeing.response, 1, quantile, prob = 0.975)/sizeing$N

freshplot4 <-
  ggplot(subset(sizeing, study.habitat == 'Freshwater')) +
  geom_ribbon(aes(x = total.length, ymin = lower95, ymax = upper95),
              fill = pal[3],
              alpha = 0.75) +
  geom_line(aes(x = total.length, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = total.length, y = IR, size = N),
             shape = 1) +
  facet_wrap( ~ region, ncol = 4,
              labeller = label_wrap_gen(width = 20)) +
  labs(x = 'Total Length (cm)',
       y = 'Microplastic Occurrence Rate',
       size = 'Sample Size') +
  scale_x_continuous(trans = 'log',
                     breaks = c(1, 10, 100),
                     limits = c(1, 500)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25)) +
  theme1

marineplot4 <-
  ggplot(subset(sizeing, study.habitat == 'Marine')) +
  geom_ribbon(aes(x = total.length, ymin = lower95, ymax = upper95),
              fill = pal[2],
              alpha = 0.75) +
  geom_line(aes(x = total.length, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = total.length, y = IR, size = N),
             shape = 1) +
  facet_wrap( ~ region, ncol = 4,
              labeller = label_wrap_gen(width = 20)) +
  labs(x = 'Total Length (cm)',
       y = 'Microplastic Occurrence Rate',
       size = 'Sample Size') +
  scale_x_continuous(trans = 'log',
                     breaks = c(1, 10, 100),
                     limits = c(1, 500)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25)) +
  theme1


png(
  'Ingestion Rate Body Size Bayesian Plot.png',
  width = 14,
  height = 15,
  units = 'cm',
  res = 500
)

plot_grid(
  freshplot4,
  marineplot4,
  labels = c('A', 'B'),
  rel_heights = c(1, 2.5),
  nrow = 2,
  align = 'v'
)

dev.off()

#### Occurrence x size mod - predictions ####

## Simulate results according to trophic level

set.seed(6364)

sizeing.sim1 <-
  data.frame(
    length = seq(
      from = 1,
      to = 200,
      length.out = 10000
    ),
    min.size = rep(100, 10000),
    sample.size = rpois(10000, 38.7),
    fibres = as.integer(rep(2, 10000)),
    region = as.factor(sample(seq(
      from = 1,
      to = max(as.integer(sizeing$region))
    ),
    10000,
    replace = TRUE))
  )

sizeing.sim1$stand.length <-
  (sizeing.sim1$length - mean(sizeing$total.length)) /
  sd(sizeing$total.length - mean(sizeing$total.length))

sizeing.sim1$stand.min.size <-
  (sizeing.sim1$min.size - mean(sizeing$min.size)) /
  sd(sizeing$min.size - mean(sizeing$min.size))

for (i in 1:10000) {
  p <-
    plogis(
      sizeingrun1$BUGSoutput$sims.list$alpha_region[, sizeing.sim1$region[i]] +
        sizeingrun1$BUGSoutput$sims.list$beta_min.size * sizeing.sim1$stand.min.size[i] +
        sizeingrun1$BUGSoutput$sims.list$beta_length[, sizeing.sim1$region[i]] *
        sizeing.sim1$stand.length[i] +
        sizeingrun1$BUGSoutput$sims.list$gamma_exclude.fibs[, sizeing.sim1$fibres[i]]
    )
  y <- rbinom(p,
              prob = p,
              size = sizeing.sim1$sample.size[i])/sizeing.sim1$sample.size[i]
  sizeing.sim1$median[i] <- median(p)
  sizeing.sim1$upper25[i] <- quantile(y, 0.625)
  sizeing.sim1$lower25[i] <- quantile(y, 0.375)
  sizeing.sim1$upper50[i] <- quantile(y, 0.75)
  sizeing.sim1$lower50[i] <- quantile(y, 0.25)
  sizeing.sim1$upper75[i] <- quantile(y, 0.875)
  sizeing.sim1$lower75[i] <- quantile(y, 0.125)
  sizeing.sim1$upper95[i] <- quantile(y, 0.975)
  sizeing.sim1$lower95[i] <- quantile(y, 0.025)
  sizeing.sim1$sample[i] <- sample(y, 1)
}

summary(sizeing.sim1)

sizeing.sim1$region <- mapvalues(
  sizeing.sim1$region,
  from = levels(sizeing.sim1$region),
  to = c(
    "Africa - Inland Waters",
    "America: North - Inland Waters",
    "Asia - Inland Waters",
    "Atlantic: Eastern Central",
    "Atlantic: Northeast",
    "Atlantic: Northwest",
    "Atlantic: Southwest",
    "Atlantic: Western Central",
    "Europe - Inland Waters",
    "Indian Ocean: Antarctic",
    "Indian Ocean: Eastern",
    "Indian Ocean: Western",
    "Mediterranean and Black Sea",
    "Pacific: Eastern Central",
    "Pacific: Northwest",
    "Pacific: Southeast",
    "Pacific: Southwest"
  )
)

sizeing.sim1$area <-
  with(sizeing.sim1,
       as.factor(
         ifelse(
           region != "Africa - Inland Waters" &
             region != "America: North - Inland Waters" &
             region != "Asia - Inland Waters" &
             region != "Europe - Inland Waters",
           "Marine",
           "Freshwater"
         )
       ))

sizeing$area <-
  as.factor(with(
    sizeing,
    ifelse(
      region != "Africa - Inland Waters" &
        region != "America: North - Inland Waters" &
        region != "Asia - Inland Waters" &
        region != "Europe - Inland Waters",
      "Marine",
      "Freshwater"
    )
  ))

sizeing.plot.a <- 
  ggplot() +
  geom_ribbon(
    data = subset(sizeing.sim1, area == "Freshwater"),
    aes(x = length,
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(sizeing.sim1, area == "Freshwater"),
    aes(x = length,
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(sizeing.sim1, area == "Freshwater"),
    aes(x = length,
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(sizeing.sim1, area == "Freshwater"),
    aes(x = length,
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(data = subset(sizeing.sim1, area == "Freshwater"),
            aes(x = length,
                y = median),
            colour = pal[5]) +
  geom_point(
    data = subset(sizeing, area == "Freshwater"),
    aes(x = total.length,
        y = IR),
    shape = 1,
    size = 0.5,
    colour = pal[5]
  ) +
  facet_wrap(~ region,
             labeller = label_wrap_gen(width = 15),
             ncol = 4) +
  labs(x = '',
       y = "") +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
                     expand = c(0, 0.05),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0.1),
                     limits = c(1, 200),
                     breaks = c(1, 100, 200)) +
  theme1

sizeing.plot.b <- 
  ggplot() +
  geom_ribbon(
    data = subset(sizeing.sim1, area == "Marine"),
    aes(x = length,
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(sizeing.sim1, area == "Marine"),
    aes(x = length,
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(sizeing.sim1, area == "Marine"),
    aes(x = length,
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = subset(sizeing.sim1, area == "Marine"),
    aes(x = length,
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(data = subset(sizeing.sim1, area == "Marine"),
            aes(x = length,
                y = median),
            colour = pal[5]) +
  geom_point(
    data = subset(sizeing, area == "Marine"),
    aes(x = total.length,
        y = IR),
    shape = 1,
    size = 0.5,
    colour = pal[5]
  ) +
  facet_wrap(~ region,
             labeller = label_wrap_gen(width = 15),
             ncol = 4) +
  labs(x = 'Total Length (cm)',
       y = "Occurrence Rate") +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25),
                     expand = c(0, 0.05),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0.1),
                     limits = c(1, 200),
                     breaks = c(1, 100, 200)) +
  theme1

png('Ingestion Rate by Size Predictions Plot.png', width = 9, height = 13, 
    units = 'cm', res = 500)

plot_grid(
  sizeing.plot.a,
  sizeing.plot.b,
  labels = c('A', 'B'),
  rel_heights = c(1, 3.2),
  nrow = 2,
  align = 'v'
)

dev.off()

#### Family mod - set up data ####

for(i in 1:nrow(size)) {
  size$famcount[i] <- nrow(subset(size, family == family[i]))
}

fam <- subset(size, famcount >= 10)

nrow(fam)  # 180 data points
summary(fam$total.length)  # 1.12-210.83 cm total length
length(unique(fam$family))  # 12 families
length(unique(fam$species))  # 133 species
length(unique(fam$study))  # 46 studies

fam$family <- as.character(fam$family)
fam$family <- as.factor(fam$family)

#### Family mod - fit model ####

fammod <- 
  function() {
    # Likelihood
    for(i in 1:N) {
      y[i] ~ dnegbin(p[i], r)
      p[i] <- r/(r + mu[i])
      log(mu[i]) <-
        log(sample.size[i]) +
        alpha_family[family[i]] +
        beta_min.size*min.size[i] +
        beta_length[family[i]]*length[i] +
        gamma_exclude.fibs[exclude.fibs[i]]
      
      ## Fitted values
      fitted[i] ~ dnegbin(p[i], r)
    }
    ## Priors
    r ~ dexp(0.1)
    
    for(j in 1:nfamily) {
      alpha_family[j] <- B[j,1]
      beta_length[j] <- B[j,2]
      B[j,1:2] ~ dmnorm(B.hat[j,], Tau.B[,])
      B.hat[j,1] <- mu_family
      B.hat[j,2] <- mu_length
    }
    mu_family ~ dnorm(0, 1)
    mu_length ~ dnorm(0, 1)
    
    Tau.B[1:2, 1:2] <- inverse(Sigma.B[,])
    Sigma.B[1,1] <- pow(sigma_family, 2)
    sigma_family ~ dexp(1)
    Sigma.B[2,2] <- pow(sigma_length, 2)
    sigma_length ~ dexp(1)
    Sigma.B[1,2] <- rho * sigma_family * sigma_length
    Sigma.B[2,1] <- Sigma.B[1,2]
    rho ~ dnorm(0, 0.5); T(-1, 1)
    
    beta_min.size ~ dnorm(-1, 1)
    
    for(k in 1:2) {
      gamma_exclude.fibs[k] ~ dnorm(0, 1)
    }
  }

## Generate initial values for MCMC

faminit <- function()
{
  list(
    "r" = dexp(0.1),
    "mu_family" = rnorm(1),
    "mu_length" = rnorm(1),
    "sigma_family" = rexp(1),
    "sigma_length" = rexp(1),
    "beta_min.size" = rnorm(1),
    "gamma_exclude.fibs" = rnorm(2),
    "rho" = 0
  )
}

## Keep track of parameters

famparam <- c("beta_min.size", "alpha_family",
              "beta_length", "gamma_exclude.fibs", "r")

## Specify data

famdata <-
  list(
    y = fam$totalcount,
    length = as.numeric(scale(fam$total.length, center = TRUE)),
    min.size = as.numeric(scale(fam$min.size, center = TRUE)),
    family = as.integer(as.factor(fam$family)),
    N = nrow(fam),
    nfamily = max(as.integer(fam$family)),
    exclude.fibs = as.integer(as.factor(fam$exclude.fib)),
    sample.size = fam$N
  )

## Run the model
famrun1 <- jags.parallel(
  data = famdata,
  inits = faminit,
  parameters.to.save = famparam,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 70000,
  n.burnin = 5000,
  n.thin = 40,
  jags.seed = 5156,
  model = fammod
)

famrun1
famrun1mcmc <- as.mcmc(famrun1)
xyplot(famrun1mcmc, layout = c(6, ceiling(nvar(famrun1mcmc)/6)))

#### Family mod - diagnostics ####

famparam2 <- "fitted"

famrun2 <- jags.parallel(
  data = famdata,
  inits = faminit,
  parameters.to.save = famparam2,
  n.chains = 3,
  n.cluster = 8,
  n.iter = 70000,
  n.burnin = 5000,
  n.thin = 40,
  jags.seed = 5156,
  model = fammod
)

fammod.response <- t(famrun2$BUGSoutput$sims.list$fitted)
fammod.observed <- fam$totalcounts
fammod.fitted <- apply(t(famrun2$BUGSoutput$sims.list$fitted),
                        1,
                        median)

check.fammod <- createDHARMa(simulatedResponse = fammod.response,
                              observedResponse = fammod.observed, 
                              fittedPredictedResponse = fammod.fitted,
                              integerResponse = T)

plot(check.fammod)

#### Family mod - inference ####

famrun_paramNames <-
  c("Acanthuridae",
    "Sciaenidae",
    "Scombridae",
    "Sparidae",
    "Carangidae",
    "Clupeidae",
    "Cyprinidae",
    "Engraulidae",
    "Gerreidae",
    "Gobiidae",
    "Mugilidae",
    "Mullidae",
    "Standarized total length:Acanthuridae",
    "Standarized total length:Sciaenidae",
    "Standarized total length:Scombridae",
    "Standarized total length:Sparidae",
    "Standarized total length:Carangidae",
    "Standarized total length:Clupeidae",
    "Standarized total length:Cyprinidae",
    "Standarized total length:Engraulidae",
    "Standarized total length:Gerreidae",
    "Standarized total length:Gobiidae",
    "Standarized total length:Mugilidae",
    "Standarized total length:Mullidae",
    "Standardized lowest detectable particle size (microns)",
    "Fibres excluded",
    "Fibres not excluded")

## Posterior density plots

famrun1long <- extract.post(famrun1)

famrun1long$variable <- mapvalues(famrun1long$variable,
                                      from = levels(famrun1long$variable),
                                      to = famrun_paramNames)

famrun1long$order <- c(nrow(famrun1long):1)

png(
  'Family Model Posteriors Plot.png',
  width = 14,
  height = 8,
  units = 'cm',
  res = 500
)

ggplot(famrun1long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, order, mean)),
    fill = pal[1],
    colour = pal[5],
    alpha = 0.75
  ) +
  geom_vline(
    aes(xintercept = 0),
    linetype = 'dashed',
    size = 0.5,
    colour = pal[3]
  ) +
  labs(x = "",
       y = "Parameter") +
  theme1

dev.off()

fam$post.predict <-
  apply(fammod.response, 1, median)
fam$lower95 <-
  apply(fammod.response, 1, quantile, prob = 0.025)
fam$upper95 <-
  apply(fammod.response, 1, quantile, prob = 0.975)

png('Gut Content Family Bayesian Plot.png', width = 9, height = 10, units = 'cm', 
    res = 500)

ggplot(fam) +
  geom_ribbon(aes(x = total.length, ymin = lower95/N, ymax = upper95/N), 
              alpha = 0.75, fill = pal[3]) +
  geom_line(aes(x = total.length, y = post.predict/N),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = total.length, y = Mpsgut),
             shape = 1, size = 0.5, alpha = 0.5) +
  facet_wrap(~ reorder(family, Mpsgut, mean), ncol = 3) +
  labs(x = expression(paste('Minimum Detectable Particle Size ('*mu*'m)')),
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       ))) +
  scale_x_continuous(trans = 'log', breaks = c(1, 10, 100, 500)) +
  scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 10, 30)) +
  theme1

dev.off()

#### Family mod - predictions ####

## Simulate predictions for different families if LDPS is held to 100 microns

set.seed(5251)

fam.sim1 <-
  expand.grid(length = mean(fam$total.length),
              min.size = 100,
              fibres = as.integer(2),
              family = unique(as.integer(fam$family)))
fam.sim1$sample.size <- rpois(nrow(fam.sim1), 38.7)

fam.sim1$stand.length <-
  (fam.sim1$length - mean(fam$total.length)) /
  sd(fam$total.length - mean(fam$total.length))

fam.sim1$stand.min.size <-
  (fam.sim1$min.size - mean(fam$min.size)) /
  sd(fam$min.size - mean(fam$min.size))

for(i in 1:nrow(fam.sim1)) {
  mu <-
    exp(
      log(fam.sim1$sample.size[i]) +
        famrun1$BUGSoutput$sims.list$alpha_family[, fam.sim1$family[i]] +
        famrun1$BUGSoutput$sims.list$beta_min.size * fam.sim1$stand.min.size[i] +
        famrun1$BUGSoutput$sims.list$beta_length[, fam.sim1$family[i]] *
        fam.sim1$stand.length[i] +
        famrun1$BUGSoutput$sims.list$gamma_exclude.fibs[, fam.sim1$fibres[i]]
    )
  r <- famrun1$BUGSoutput$sims.list$r
  p <- r / (mu + r)
  y <- rnbinom(p,
               prob = p,
               size = r) / fam.sim1$sample.size[i]
  fam.sim1$median[i] <- median(mu ) / fam.sim1$sample.size[i]
  fam.sim1$upper25[i] <- quantile(y, 0.625)
  fam.sim1$lower25[i] <- quantile(y, 0.375)
  fam.sim1$upper50[i] <- quantile(y, 0.75)
  fam.sim1$lower50[i] <- quantile(y, 0.25)
  fam.sim1$upper75[i] <- quantile(y, 0.875)
  fam.sim1$lower75[i] <- quantile(y, 0.125)
  fam.sim1$upper95[i] <- quantile(y, 0.975)
  fam.sim1$lower95[i] <- quantile(y, 0.025)
}

summary(fam.sim1)

fam.sim1$family <- as.factor(fam.sim1$family)

fam.sim1$family <- mapvalues(fam.sim1$family,
                            from = levels(fam.sim1$family),
                            to = c("Acanthuridae",
                                   "Carangidae",
                                   "Clupeidae",
                                   "Cyprinidae",
                                   "Engraulidae",
                                   "Gerreidae",
                                   "Gobiidae",
                                   "Mugilidae",
                                   "Mullidae",
                                   "Sciaenidae",
                                   "Scombridae",
                                   "Sparidae"))

png(
  'Family Model Predictions Plot.png',
  width = 9,
  height = 9,
  units = 'cm',
  res = 500
)

ggplot() +
  geom_linerange(
    data = fam.sim1,
    aes(x = reorder(family, median, mean),
        ymin = lower25,
        ymax = upper25),
    alpha = 0.75,
    colour = pal[1],
    size = 1
  ) +
  geom_linerange(
    data = fam.sim1,
    aes(x = reorder(family, median, mean),
        ymin = lower50,
        ymax = upper50),
    alpha = 0.5,
    colour = pal[1],
    size = 1
  ) +
  geom_linerange(
    data = fam.sim1,
    aes(x = reorder(family, median, mean),
        ymin = lower75,
        ymax = upper75),
    alpha = 0.25,
    colour = pal[1],
    size = 1
  ) +
  geom_linerange(
    data = fam.sim1,
    aes(x = reorder(family, median, mean),
        ymin = lower95,
        ymax = upper95),
    alpha = 0.05,
    colour = pal[1],
    size = 1
  ) +
  geom_point(data = fam.sim1,
            aes(x = reorder(family, median, mean),
                y = median),
            colour = pal[5],
            fill = pal[3],
            size = 3,
            shape = 21) +
  geom_jitter(
    data = fam,
    aes(x = reorder(family, Mpsgut, mean),
        y = Mpsgut),
    size = 0.75,
    colour = pal[5],
    shape = 1
  ) +
  labs(x = 'Family',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       ))) +
  scale_y_continuous(trans = 'log1p', 
                     expand = c(0, 0), 
                     limits = c(-0.1, 60),
                     breaks = c(0, 1, 10, 50)) +
  theme1 +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

dev.off()

#### Additional plots ####

## Plot according to lower limit of detection

png('Lower Limit Plot.png', width = 19, height = 11, 
    units = 'cm', res = 600)

ggplot(gutdata) +
  geom_point(aes(
    x = reorder(region, Mpsgut, mean),
    y = Mpsgut,
    size = N,
    colour = min.size
  ),
  alpha = 0.5) +
  facet_grid(area ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = 'FAO Area',
    y = expression(paste(
      'Microplastic Concentration (particles ' ~
        ind ^ -1 * ')'
    )),
    colour = expression(paste('Lowest Detectable Particle Size ('*mu*'m)')),
    size = 'Sample Size'
  ) +
  scale_y_continuous(
    breaks = c(0, seq(1, 9, 1), 10, 15, seq(20, 40, 10)),
    expand = c(0, 0.05),
    trans = 'log1p'
  ) +
  coord_flip() +
  scale_colour_gradient2(low = pal[3], mid = pal[2], high = pal[1],
                         midpoint = 500, 
                         breaks = c(0.2, 250, 500, 750, 1000)) +
  scale_size(breaks = c(1, seq(from = 200, to = 1400, by = 200))) +
  theme1 +
  theme(legend.margin = margin(0, 0, 0, 0, unit = 'cm'),
        panel.grid.major.x = element_line(colour = pal[5], 
                                          linetype = "dashed",
                                          size = 0.25),
        legend.position = "bottom",
        legend.spacing = unit(1, units = "cm"),
        legend.justification = "left")

dev.off()

## Plot MP concentrations by region

png('MP Occurrence Rate by Region.png', width = 19, height = 12, 
    units = 'cm', res = 700)

ggplot(ingestion) + 
  geom_density_ridges(aes(x = IR,
                          y = reorder(region, IR, mean)),
                      fill = pal[2],
                      colour = pal[5],
                      alpha = 0.75,
                      scale = 1) +
  geom_point(aes(x = IR,
                 y = region),
             size = 1.5, alpha = 0.3, colour = pal[3]) +
  facet_grid(area ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0.01)) +
  labs(x = "Occurrence Rate",
       y = "FAO Area") +
  theme1 +
  theme(panel.grid.major.y = element_line(colour = pal[5], size = 0.25))

dev.off()

