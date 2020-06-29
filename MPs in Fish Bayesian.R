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
                        'family', 'genus', 'environment', 'climate', 'red.list', 
                        'feeding.type', 'feeding.habit', 'TL', 'min.size', 
                        'float.meth', 'dig.meth', 'count.meth', 'polymer.meth',
                        'polymer.ID', 'blanks', 'maj.fib', 'maj.under.one.mm', 
                        'maj.polymer', 'maj.col', 'exclude.fib'), 
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
            to = c('Polymer ID Method Used',
                   'Polymer ID Method Not Used'))
trophicfish2$blanks <- 
  mapvalues(trophicfish2$blanks,
            from = c('yes', 'no'),
            to = c('Blanks Used',
                   'Blanks Not Used'))

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
trophicfish2$area <- 
  mapvalues(trophicfish2$region,
            from = levels(trophicfish2$region),
            to = c('Freshwater', 'Freshwater',
                   'Freshwater', 'Freshwater',
                   'Atlantic', 'Atlantic',
                   'Atlantic', 'Atlantic',
                   'Freshwater', 'Indian Ocean',
                   'Indian Ocean', 'Indian Ocean',
                   'Mediterranean and Black Sea', 'Pacific',
                   'Pacific', 'Pacific',
                   'Pacific', 'Pacific',
                   'Pacific'))

#### MPs in Guts Model ####

## Set up the data
gutdata <- subset(trophicfish2, Mpsgut != 'NA' & environment != '')
summary(gutdata)
summary(gutdata$author)
length(gutdata$species) # 731 data points
length(unique(gutdata$species)) # 553 species
length(unique(gutdata$family)) # from 158 families
length(unique(gutdata$study)) # from 104 + 1 studies

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

## Specify model

jagsmod1 <- function()
{
  # Likelihood
  for (i in 1:N)
  {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <-
      alpha + beta_TL * TL[i] + beta_PID[PID[i]] + beta_min.size * min.size[i] +
      beta_exclude.fibs[exclude.fibs[i]] + beta_blanks[blanks[i]] +
      beta_sample.size * sample.size[i] + beta_environment[environment[i]] +
      beta_region[region[i]] + beta_interaction[region[i]] * TL[i]
  }
  
  # Prior
  for (j in 1:nregion)
  {
    beta_region[j] ~ dnorm(mu_region, tau_region)
    beta_interaction[j] ~ dnorm(mu_interaction, tau_interaction)
  }
  alpha ~ dexp(1)
  sigma ~ dexp(1)
  tau <- 1 / (sigma * sigma)
  beta_TL ~ dnorm(0, 1)
  beta_min.size ~ dnorm(-1, 1)
  beta_sample.size ~ dnorm(0, 1)
  for (k in 1:2)
  {
    beta_PID[k] ~ dnorm(0, 1)
    beta_exclude.fibs[k] ~ dnorm(0, 1)
    beta_blanks[k] ~ dnorm(0, 1)
  }
  for (l in 1:nenvironment)
  {
    beta_environment[l] ~ dnorm(0, 1)
  }
  mu_region ~ dnorm(0, 1)
  sigma_region ~ dexp(1)
  tau_region <- 1 / (sigma_region * sigma_region)
  mu_interaction ~ dnorm(0, 1)
  sigma_interaction ~ dexp(1)
  tau_interaction <- 1 / (sigma_interaction * sigma_interaction)
}

## Generate initial values for MCMC

init1 <- function()
{
  list(
    "sigma" = 1,
    "alpha" = 1,
    "beta_TL" = rnorm(1),
    "beta_min.size" = rnorm(1),
    "beta_sample.size" = rnorm(1),
    "beta_PID" = rnorm(2),
    "beta_exclude.fibs" = rnorm(2),
    "beta_blanks" = rnorm(2),
    "beta_environment" = rnorm(12),
    "mu_region" = rnorm(1),
    "sigma_region" = 1,
    "mu_interaction" = rnorm(1),
    "sigma_interaction" = 1
  )
}

## Keep track of parameters

param1 <- c("sigma", "alpha", "beta_TL", "beta_min.size", "beta_sample.size", 
            "beta_PID", "beta_exclude.fibs", "beta_blanks", "beta_environment",
            "beta_region", "beta_interaction")

## Specify data

jagsdata1 <-
  list(
    y = log(gutdata$Mpsgut + 1),
    TL = as.numeric(scale(gutdata$TL, center = TRUE)),
    PID = as.integer(gutdata$polymer.ID),
    min.size = as.numeric(scale(gutdata$min.size, center = TRUE)),
    exclude.fibs = as.integer(gutdata$exclude.fib),
    blanks = as.integer(gutdata$blanks),
    sample.size = as.numeric(scale(gutdata$N, center = TRUE)),
    region = as.integer(gutdata$region),
    environment = as.integer(gutdata$environment),
    N = nrow(gutdata),
    nenvironment = max(as.integer(gutdata$environment)),
    nregion = max(as.integer(gutdata$region))
  )

## Run the model
run1 <- jags(
  data = jagsdata1,
  inits = init1,
  parameters.to.save = param1,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = jagsmod1
)
beep(1)

run1mcmc <- as.mcmc(run1)
xyplot(run1mcmc, layout = c(6, ceiling(nvar(run1mcmc)/6)))

## Extend burnin to 5000 and up number of iterations to 10,000

run2 <- jags(
  data = jagsdata1,
  inits = init1,
  parameters.to.save = param1,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 5000,
  n.thin = 5,
  jags.seed = 123,
  model = jagsmod1
)
beep(1)

run2mcmc <- as.mcmc(run2)
xyplot(run2mcmc, layout = c(6, ceiling(nvar(run1mcmc)/6)))

## Increase iteration to 500000 and burnin to 10000

run3 <- jags(
  data = jagsdata1,
  inits = init1,
  parameters.to.save = param1,
  n.chains = 3,
  n.iter = 500000,
  n.burnin = 5000,
  n.thin = 45,
  jags.seed = 123,
  model = jagsmod1
)
beep(1)

run3mcmc <- as.mcmc(run3)
xyplot(run3mcmc, layout = c(6, ceiling(nvar(run1mcmc)/6)))  # looks good

## Inference

HPD1 <- as.data.frame(summary(run3mcmc)$quantiles)
HPD1 <- HPD1[1:60, ]

MAP1 <- as.data.frame(summary(run3mcmc)$statistics)
MAP1 <- MAP1[1:60, ]

## Extract MAP and HPDIs for the parameters
Params1 <- data.frame(
  parameter = as.factor(rownames(MAP1)),
  MAP = MAP1[, 1],
  lower = HPD1[, 1],
  upper = HPD1[, 5]
)

with(gutdata, tapply(as.integer(polymer.ID), polymer.ID, mean))
with(gutdata, tapply(as.integer(blanks), blanks, mean))
with(gutdata, tapply(as.integer(environment), environment, mean))
with(gutdata, tapply(as.integer(exclude.fib), exclude.fib, mean))
with(gutdata, tapply(as.integer(region), region, mean))

Params1$parameter <- as.character(Params1$parameter)
Params1$parameter <- as.factor(Params1$parameter)

Params1$parameter <- mapvalues(
  Params1$parameter,
  from = levels(Params1$parameter),
  to = c(
    "Intercept",
    "Blanks not used",
    "Blanks used",
    "Bathypelagic (environment)",
    "Pelagic-neritic (environment",
    "Pelagic-oceanic (environment)",
    "Reef-associated (environment)",
    "Bathydemersal (environment)",
    "Benthopelagic (environment",
    "Demersal (environment)",
    "Freshwater benthopelagic (environment)",
    "Freshwater demersal (environment)",
    "Freshwater pelagic (environment)",
    "Freshwater pelagic-neritic (environment",
    "Pelagic (environment)",
    "Fibres not excluded",
    "Fibres excluded",
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
    "Standardized trophic level:America, South - Inland Waters",
    "Standardized trophic level:Asia - Inland Waters",
    "Standardized trophic level:Atlantic, Eastern Central",
    "Standardized trophic level:Atlantic, Northeast",
    "Standardized trophic level:Atlantic, Southwest",
    "Standardized trophic level:Atlantic, Western Central",
    "Standardized trophic level:Europe - Inland Waters",
    "Standardized lowest detectable particle size (microns)",
    "Polymer ID not used",
    "Polymer ID used",
    "Africa - Inland Waters (region)",
    "Indian Ocean, Antarctic (region)",
    "Indian Ocean, Eastern (region)",
    "Indian Ocean, Western (region)",
    "Mediterranean and Black Sea (region)",
    "Pacific, Eastern Central (region)",
    "Pacific, Northeast (region)",
    "Pacific, Northwest (region)",
    "Pacific, Southeast (region)",
    "Pacific, Southwest (region)",
    "Pacific, Western Central (region)",
    "America, North - Inland Waters (region)",
    "America, South - Inland Waters (region)",
    "Asia - Inland Waters (region)",
    "Atlantic, Eastern Central (region)",
    "Atlantic, Northeast (region)",
    "Atlantic, Southwest (region)",
    "Atlantic, Western Central (region)",
    "Europe - Inland Waters (region)",
    "Standardized sample size (number of fish)",
    "Standardized trophic level"
  )
)

Params1$order <- c(nrow(Params1):1)

png('Gut Content HPDI Plot.png', 
    width = 14, 
    height = 15, 
    units = 'cm', 
    res = 500)

ggplot(Params1) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed',
             size = 0.5,
             colour = pal[3]) +
  geom_errorbar(aes(x = reorder(parameter, as.numeric(rownames(Params1))),
                    ymin = lower,
                    ymax = upper),
                size = 0.25) +
  geom_point(aes(x = reorder(parameter, as.numeric(rownames(Params1))),
                 y = MAP),
             size = 1,
             shape = 16,
             colour = pal[1]) +
  labs(x = 'Coefficient',
       y = '') +
  coord_flip() +
  theme1

dev.off()

## Run again and estimate mu this time as well

param2 <- c("mu")

run4 <- jags(
  data = jagsdata1,
  inits = init1,
  parameters.to.save = param2,
  n.chains = 3,
  n.iter = 500000,
  n.burnin = 5000,
  n.thin = 24,
  jags.seed = 123,
  model = jagsmod1
)
beep(1)

run4mcmc <- as.mcmc(run4)

# set.seed(1234)
# for(j in 1:nrow(gutdata)) {
#   for(i in 1:nrow(mu1)) {
#   predict[i] <- rnorm(1, mean = mu1[i, j], sd = sigma1[i, ])
#   }
#   gutdata$post.predict[j] <- exp(mean(predict)) - 1
#   gutdata$lower95[j] <- exp(quantile(predict, 0.025)) - 1
#   gutdata$upper95[j] <- exp(quantile(predict, 0.975)) - 1
# }
# 
# set.seed(1234)
# for(j in 1:nrow(gutdata)) {
#   for(i in 1:nrow(mu1)) {
#     mu[i] <- mu1[i, j]
#   }
#   gutdata$post.predict[j] <- exp(mean(mu)) - 1
#   gutdata$lower95[j] <- exp(quantile(mu, 0.025)) - 1
#   gutdata$upper95[j] <- exp(quantile(mu, 0.975)) - 1
# }

gutdata$post.predict <-
  exp(as.data.frame(summary(run4mcmc)$statistics)[2:728, 1]) - 1
gutdata$lower95 <-
  exp(as.data.frame(summary(run4mcmc)$quantiles)[2:728, 1]) - 1
gutdata$upper95 <-
  exp(as.data.frame(summary(run4mcmc)$quantiles)[2:728, 5]) - 1

# gutdata$post.predict[gutdata$post.predict < 0] <- 0
# gutdata$lower95[gutdata$lower95 < 0] <- 0


freshplot1 <-
  ggplot(subset(gutdata, study.habitat == 'Freshwater')) +
  geom_line(aes(x = TL, y = post.predict, colour = exclude.fib),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95, fill = exclude.fib), 
              alpha = 0.3) +
  geom_point(aes(x = TL, y = Mpsgut, colour = exclude.fib),
             shape = 1, size = 0.5) +
  facet_wrap(~ region, scales = 'free_x', ncol = 5,
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
  geom_line(aes(x = TL, y = post.predict, colour = exclude.fib),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95, fill = exclude.fib), 
              alpha = 0.3) +
  geom_point(aes(x = TL, y = Mpsgut, colour = exclude.fib),
             shape = 1, size = 0.5) +
  facet_wrap(~ region, scales = 'free_x', ncol = 5,
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
  geom_line(aes(x = min.size, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = min.size, ymin = lower95, ymax = upper95), 
              alpha = 0.3, fill = pal[1]) +
  geom_point(aes(x = min.size, y = Mpsgut),
             shape = 1, size = 0.5) +
  facet_wrap(~ region, ncol = 5,
             labeller = label_wrap_gen(width = 15)) +
  labs(x = expression(paste('Minimum Detectable Particle Size ('*mu*'m)')),
       y = expression(paste(
         '          Microplastic\nConcentration (particles ' ~ ind ^ -1 * ')'
       )),
       size = 'Sample Size') +
  scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 10, 30)) +
  scale_x_continuous(trans= 'log1p', breaks = c(0, 1, 10, 100, 500)) +
  theme1

marineplot2 <-
  ggplot(subset(gutdata, study.habitat == 'Marine')) +
  geom_line(aes(x = min.size, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = min.size, ymin = lower95, ymax = upper95), 
              alpha = 0.3, fill = pal[1]) +
  geom_point(aes(x = min.size, y = Mpsgut),
             shape = 1, size = 0.5) +
  facet_wrap(~ region, ncol = 5,
             labeller = label_wrap_gen(width = 15)) +
  labs(x = expression(paste('Minimum Detectable Particle Size ('*mu*'m)')),
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       )),
       size = 'Sample Size') +
  scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 10, 30)) +
  scale_x_continuous(trans= 'log1p', breaks = c(0, 1, 10, 100, 500)) +
  theme1

png(
  'Gut Content Bayesian Plot According to Lower Size Limit.png',
  width = 14,
  height = 16,
  units = 'cm',
  res = 500
)

plot_grid(freshplot2, marineplot2, labels = c('A', 'B'), rel_heights = c(1,2.9),
          nrow = 2, align = 'v')

dev.off()

#### Ingestion Rate Model ####

## Set up the data

ingestion <- subset(trophicfish2, !is.na(IR))

## Convert to successes/failures

ingestion$successes <- round(with(ingestion, N * IR), digits = 0)
ingestion$failures <- round(with(ingestion,  N * (1 - IR), digits = 0))

length(ingestion$species) # 638 data points
length(unique(ingestion$species)) # 478 species
length(unique(ingestion$family)) # from 157 families
length(unique(ingestion$study)) # 105 studies

ingestion$region <- as.character(ingestion$region)
ingestion$region <- as.factor(ingestion$region)

summary(ingestion)

## Specify model
jagsmod2 <- function()
{
  # Likelihood
  for(i in 1:N)
  {
    y[i] ~ dbinom(p[i], n[i])
    logit(p[i]) <- alpha + beta_TL * TL[i] + beta_region[region[i]] + 
    beta_interaction[region[i]] * TL[i]
  }
  
  # Prior
  alpha ~ dnorm(0, 1)
  beta_TL ~ dnorm(0, 1)
  for(j in 1:Nregions)
  {
    beta_region[j] ~ dnorm(mu_region, tau_region)
    beta_interaction[j] ~ dnorm(mu_interaction, tau_interaction)
  }
  mu_region ~ dnorm(0, 1)
  tau_region <- 1 / (sigma_region * sigma_region)
  sigma_region ~ dexp(1)
  mu_interaction ~ dnorm(0, 1)
  tau_interaction <- 1 / (sigma_interaction * sigma_interaction)
  sigma_interaction ~ dexp(1)
}


## Initial values for MCMC chains
init2 <- function()
{
  list(
    "alpha" = 0.001,
    "beta_TL" = rnorm(1),
    "mu_region" = rnorm(1),
    "sigma_region" = 1,
    "mu_interaction" = 1,
    "sigma_interaction" = 1
  )
}

## Parameters to keep track of
param2 <- c("alpha", "beta_TL", "beta_region", "beta_interaction")


## Specify data
jagsdata2 <- list(N = nrow(ingestion),
                  Nregions = max(as.integer(ingestion$region)),
                  y = as.numeric(ingestion$successes),
                  n = as.numeric(ingestion$N),
                  TL = as.numeric(scale(ingestion$TL, center = TRUE)),
                  region = as.integer(ingestion$region))


## Run the model
ingrun1 <- jags(
  data = jagsdata2,
  inits = init2,
  parameters.to.save = param2,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = jagsmod2
)

ingrun1mcmc <- as.mcmc(ingrun1)
xyplot(ingrun1mcmc, layout = c(6, ceiling(nvar(ingrun1mcmc)/6)))


# Update burn in to 2000 and iterations to 50000
ingrun2 <- jags(
  data = jagsdata2,
  inits = init2,
  parameters.to.save = param2,
  n.chains = 3,
  n.iter = 50000,
  n.burnin = 2000,
  n.thin = 48,
  jags.seed = 123,
  model = jagsmod2
)

ingrun2mcmc <- as.mcmc(ingrun2)
xyplot(ingrun2mcmc, 
       layout = c(6, ceiling(nvar(ingrun2mcmc)/6)))  # Need more iterations

## Increase iterations to 200000

ingrun3 <- jags(
  data = jagsdata2,
  inits = init2,
  parameters.to.save = param2,
  n.chains = 3,
  n.iter = 200000,
  n.burnin = 5000,
  n.thin = 100,
  jags.seed = 123,
  model = jagsmod2
)
beep(1)

ingrun3mcmc <- as.mcmc(ingrun3)
xyplot(ingrun3mcmc, 
       layout = c(6, ceiling(nvar(ingrun3mcmc)/6)))

## Inference

HPD2 <- as.data.frame(summary(ingrun3mcmc)$quantiles)

MAP2 <- as.data.frame(summary(ingrun3mcmc)$statistics)

## Extract MAP and HPDIs for the parameters
Params2 <- data.frame(
  parameter = as.factor(rownames(MAP2)),
  MAP = MAP2[, 1],
  lower = HPD2[, 1],
  upper = HPD2[, 5]
)

with(ingestion, tapply(as.integer(region), region, mean))

Params2 <- Params2[-39, ]
Params2$parameter <- as.character(Params2$parameter)
Params2$parameter <- as.factor(Params2$parameter)

Params2$parameter <- mapvalues(
  Params2$parameter,
  from = levels(Params2$parameter),
  to = c(
    "Intercept",
    "Standardized trophic level:Africa - Inland Waters",
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
    "Standardized trophic level:Atlantic, Southwest",
    "Standardized trophic level:Atlantic, Western Central",
    "Standardized trophic level:Europe - Inland Waters",
    "Standardized trophic level:Indian Ocean, Antarctic",
    "Africa - Inland Waters (region)",
    "Indian Ocean, Eastern (region)",
    "Indian Ocean, Western (region)",
    "Mediterranean and Black Sea (region)",
    "Pacific, Eastern Central (region)",
    "Pacific, Northeast (region)",
    "Pacific, Northwest (region)",
    "Pacific, Southeast (region)",
    "Pacific, Southwest (region)",
    "Pacific, Western Central (region)",
    "America, North - Inland Waters (region)",
    "Asia - Inland Waters (region)",
    "Atlantic, Eastern Central (region)",
    "Atlantic, Northeast (region)",
    "Atlantic, Southwest (region)",
    "Atlantic, Western Central (region)",
    "Europe - Inland Waters (region)",
    "Indian Ocean, Antarctic (region)",
    "Standardized trophic level"
  )
)

png('Ingestion HPDI Plot.png', 
    width = 14, 
    height = 15, 
    units = 'cm', 
    res = 500)

ggplot(Params2) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed',
             size = 0.5,
             colour = pal[2]) +
  geom_errorbar(aes(x = parameter,
                    ymin = lower,
                    ymax = upper),
                size = 0.25) +
  geom_point(aes(x = parameter,
                 y = MAP),
             size = 1.5,
             shape = 16,
             colour = pal[1]) +
  labs(x = 'Coefficient',
       y = '') +
  coord_flip() +
  theme1

dev.off()

## Rerun model to extract p
param3 <- c("alpha", "beta_TL", "beta_region", "beta_interaction", "p")

ingrun3 <- jags(
  data = jagsdata2,
  inits = init2,
  parameters.to.save = param3,
  n.chains = 3,
  n.iter = 200000,
  n.burnin = 5000,
  n.thin = 100,
  jags.seed = 123,
  model = jagsmod2
)

ingrun3mcmc <- as.mcmc(ingrun3)

ingestion$post.predict <-
  as.data.frame(summary(ingrun3mcmc)$statistics)[40:677, 1]
ingestion$lower95 <-
  as.data.frame(summary(ingrun3mcmc)$quantiles)[40:677, 1]
ingestion$upper95 <-
  as.data.frame(summary(ingrun3mcmc)$quantiles)[40:677, 5]

freshplot3 <-
  ggplot(subset(ingestion, study.habitat == 'Freshwater')) +
  geom_line(aes(x = TL, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95),
              fill = pal[1],
              alpha = 0.3) +
  geom_point(aes(x = TL, y = IR, size = N),
             shape = 1) +
  facet_wrap( ~ region, ncol = 4,
              labeller = label_wrap_gen(width = 20)) +
  labs(x = 'Trophic Level',
       y = 'Ingestion Rate',
       size = 'Sample Size') +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25)) +
  theme1

marineplot3 <-
  ggplot(subset(ingestion, study.habitat == 'Marine')) +
  geom_line(aes(x = TL, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95),
              fill = pal[1],
              alpha = 0.3) +
  geom_point(aes(x = TL, y = IR, size = N),
             shape = 1) +
  facet_wrap( ~ region, ncol = 4,
              labeller = label_wrap_gen(width = 20)) +
  labs(x = 'Trophic Level',
       y = 'Ingestion Rate',
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


#### Allometric model ####

## Set up the data
allo <- subset(gutdata, total.length != 'NA')

length(allo$total.length)  # 342 data point remaining
length(unique(allo$study))  # 52 studies
length(unique(allo$species))  # 300 species

## Specify model
allomod <- function()
{
  # Likelihood
  for (i in 1:N)
  {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta_min.size * min.size[i] + beta_length * length[i]
  }
  # Prior
  alpha ~ dexp(1)
  beta_min.size ~ dnorm(0, 1)
  beta_length ~ dnorm(0, 1)
  tau <- 1 / (sigma * sigma)
  sigma ~ dexp(1)
}

## Initial values
allomod_init <- function()
{
  list(
    "alpha" = 1,
    "beta_min.size" = rnorm(1),
    "beta_length" = rnorm(1),
    "sigma" = 1
  )
}

## Parameters to keep track of
allomod_params <- c("alpha", "beta_min.size", "beta_length", "sigma")

## Specify data

allodata <- list(
  N = nrow(allo),
  y = as.numeric(log(allo$Mpsgut + 1)),
  min.size = as.numeric(scale(allo$min.size, center = TRUE)),
  length = as.numeric(scale(allo$total.length, center = TRUE))
)

## Run the model
allorun1 <- jags(
  data = allodata,
  inits = allomod_init,
  parameters.to.save = allomod_params,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = allomod
)

allorun1mcmc <- as.mcmc(allorun1)
xyplot(allorun1mcmc)

## Increase to 10,000 interation

allorun2 <- jags(
  data = allodata,
  inits = allomod_init,
  parameters.to.save = allomod_params,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = allomod
)

allorun2  # DIC = 753
allorun2mcmc <- as.mcmc(allorun2)
summary(allorun2mcmc)
xyplot(allorun2mcmc)

## Inference

alloHPDI <- as.data.frame(summary(allorun2mcmc)$quantiles)

alloMAP <- as.data.frame(summary(allorun2mcmc)$statistics)

## Extract MAP and HPDIs for the parameters
alloParams <- data.frame(
  parameter = as.factor(rownames(alloMAP)),
  MAP = alloMAP[, 1],
  lower = alloHPDI[, 1],
  upper = alloHPDI[, 5]
)

alloParams <- alloParams[-c(4:5), ]
alloParams$parameter <- as.character(alloParams$parameter)
alloParams$parameter <- as.factor(alloParams$parameter)

alloParams$parameter <- mapvalues(
  alloParams$parameter,
  from = levels(alloParams$parameter),
  to = c(
    "Intercept",
    "Standardized Total Length (cm)",
    "Standardized lowest detectable particle size (microns)"
  )
)

png('Allometry HPDI Plot.png', 
    width = 9, 
    height = 2.25, 
    units = 'cm', 
    res = 500)

ggplot(alloParams) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed',
             size = 0.5,
             colour = pal[2]) +
  geom_errorbar(aes(x = parameter,
                    ymin = lower,
                    ymax = upper),
                size = 0.25) +
  geom_point(aes(x = parameter,
                 y = MAP),
             size = 1,
             shape = 16,
             colour = pal[1]) +
  labs(x = 'Parameter',
       y = '') +
  coord_flip() +
  theme1

dev.off()

## Rerun model and extract mu

allomod_params1 <- c("mu")

allorun3 <- jags(
  data = allodata,
  inits = allomod_init,
  parameters.to.save = allomod_params1,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = allomod
)

allorun3
allorun3mcmc <- as.mcmc(allorun3)

allo$post.predict <- 
  exp(as.data.frame(summary(allorun3mcmc)$statistics)[-1, 1]) - 1
allo$lower95 <- 
  exp(as.data.frame(summary(allorun3mcmc)$quantiles)[-1, 1]) - 1
allo$upper95 <- 
  exp(as.data.frame(summary(allorun3mcmc)$quantiles)[-1, 5]) - 1

## Plot
png(
  'Allometry Bayesian Plot.png',
  width = 9,
  height = 11,
  units = 'cm',
  res = 500
)

ggplot(allo) +
  geom_line(aes(x = total.length, y = post.predict),
            size = 0.5,
            alpha = 0.8) +
  geom_ribbon(
    aes(x = total.length, ymin = lower95, ymax = upper95),
    fill = pal[1],
    alpha = 0.3
  ) +
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

## Simulate results for 1 micron filter size

set.seed(4516)
simdata <-
  data.frame(length = runif(1000, 1, 500),
             min.size = rep(500, 1000))

simdata$stand.length <-
  (simdata$length - mean(allo$total.length)) /
  sd(allo$total.length - mean(allo$total.length))

simdata$stand.min.size <-
  (simdata$min.size - mean(allo$min.size)) /
  sd(allo$min.size - mean(allo$min.size))

simdata$mean <- NA
simdata$upper95 <- NA
simdata$lower95 <- NA

alpha <- rexp(1000, 
              summary(allorun2mcmc)$statistics[1, 1])
beta_length <- rnorm(1000, 
                     summary(allorun2mcmc)$statistics[2,1],
                     summary(allorun2mcmc)$statistics[2,2])
beta_min.size <- rnorm(1000, 
                       summary(allorun2mcmc)$statistics[3,1],
                       summary(allorun2mcmc)$statistics[3,2])
sigma <- dexp(1000,
              summary(allorun2mcmc)$statistics[5,1])

for(i in 1:1000) {
  mu <-
    alpha + beta_length * simdata$stand.length[i] + 
    beta_min.size * simdata$stand.min.size[i]
  y <- exp(rnorm(1000, mu, sigma)) - 1
  simdata$mean[i] <- median(y)
  simdata$upper95[i] <- quantile(y, 0.975)
  simdata$lower95[i] <- quantile(y, 0.025)
}

simdata$predict <-
  exp(
    rnorm(
      1000,
      alpha +
        beta_length * simdata$stand.length +
        beta_min.size * simdata$stand.min.size,
      sigma
    )
  ) - 1

png(
  'Allometry Predictions Plot.png',
  width = 9,
  height = 9,
  units = 'cm',
  res = 500
)

ggplot(simdata) +
  geom_ribbon(aes(x = length,
                  ymin = lower95,
                  ymax = upper95),
              alpha = 0.3, fill = pal[1]) +
  geom_line(aes(x = length,
                 y = mean),
            colour = pal[1]) +
  geom_point(aes(x = length,
                 y = predict),
             size = 0.5) +
  labs(x = 'Total Length (cm)',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'))) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 10, 100, 1000)) +
  theme1

dev.off()
