#### Load packages ####
library(ggplot2)
library(cowplot)
library(plyr)
library(colorspace)
library(beepr)
library(ggridges)
library(reshape2)
library(DHARMa)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(brms)
library(bayesplot)
library(tidybayes)
color_scheme_set("viridis")

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

# MPs in guts model -------------------------------------------------------


#### Set up the data ####

gutdata <- subset(trophicfish2, Mpsgut != 'NA' & environment != '')
summary(gutdata)
summary(gutdata$author)
length(gutdata$species) # 730 data points
length(unique(gutdata$species)) # 550 species
length(unique(gutdata$family)) # from 157 families
length(unique(gutdata$study)) # from 104 studies

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

#### Fit model ####

model1 <-
  brm(bf(totalcount ~
           polymer.ID + blanks + exclude.fib +
           offset(log(N)) +
           scale(TL, center = TRUE) +
           scale(min.size, center = TRUE) +
           (1 + TL | region) + (1 | environment),
         zi ~ scale(min.size, center = TRUE)),
      data = gutdata,
      family = zero_inflated_poisson(link = "log",
                                     link_zi = "logit"),
      prior = c(set_prior("normal(0, 1)", class = "b"),
                set_prior("lkj(2)", class = "cor"),
                set_prior("normal(0, 1)", class = "Intercept"),
                set_prior("exponential(1)", class = "sd"),
                set_prior("normal(0, 1)", class = "b", dpar = "zi")),
    warmup = 500,
    iter = 2000,
    chains = 3,
    inits = "random",
    cores = 8,
    seed = 8913,
    control = list(max_treedepth = 15))

summary(model1)

## Compare without ZI

model2 <-
  brm(bf(totalcount ~
           polymer.ID + blanks + exclude.fib +
           offset(log(N)) +
           scale(TL, center = TRUE) +
           scale(min.size, center = TRUE) +
           (1 + TL | region) + (1 | environment),
         zi ~ scale(min.size, center = TRUE) +
           scale(N, center = TRUE)),
      data = gutdata,
      family = zero_inflated_poisson(link = "log",
                                     link_zi = "logit"),
      prior = c(set_prior("normal(0, 1)", class = "b"),
                set_prior("lkj(2)", class = "cor"),
                set_prior("normal(0, 1)", class = "Intercept"),
                set_prior("exponential(1)", class = "sd"),
                set_prior("normal(0, 1)", class = "b", dpar = "zi")),
      warmup = 500,
      iter = 2000,
      chains = 3,
      inits = "random",
      cores = 8,
      seed = 8913,
      control = list(max_treedepth = 15))

summary(model2)

loo(model1)
loo(model2)

loo_compare(model1, model2)

waic(model1)
waic(model2)

#### Diagnostics ####

mcmc_trace(model1)
mcmc_trace(model2)

set.seed(51685)

response1 <- t(posterior_predict(model1, draws = 5000))
fitted1 <- apply(t(posterior_predict(model1, draws = 5000)), 1, median)

diagnose1 <- createDHARMa(simulatedResponse = response1,
                          observedResponse = gutdata$totalcount,
                          fittedPredictedResponse = fitted1)
plot(diagnose1)

response2 <- t(posterior_predict(model2, draws = 5000))
fitted2 <- apply(t(posterior_predict(model2, draws = 5000)), 1, median)

diagnose2 <- createDHARMa(simulatedResponse = response2,
                          observedResponse = gutdata$totalcount,
                          fittedPredictedResponse = fitted2)
plot(diagnose2)

#### Inference ####

HPD1 <- as.data.frame(posterior_summary(model1))
HPD1$param <- rownames(HPD1)
rownames(HPD1)
HPD1 <- HPD1[-c(1:3,), ]

MAP1 <- as.data.frame(summary(run2mcmc)$statistics)
MAP1 <- MAP1[c(1:19, 21:41, 44:61), ]

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

gutmod_paramnames <-
  c(
    "Africa - Inland Waters (FAO area)",
    "Indian Ocean, Antarctic (FAO area)",
    "Indian Ocean, Eastern (FAO area)",
    "Indian Ocean, Western (FAO area)",
    "Mediterranean and Black Sea (FAO area)",
    "Pacific, Eastern Central (FAO area)",
    "Pacific, Northeast (FAO area)",
    "Pacific, Northwest (FAO area)",
    "Pacific, Southeast (FAO area)",
    "Pacific, Southwest (FAO area)",
    "Pacific, Western Central (FAO area)",
    "America, North - Inland Waters (FAO area)",
    "America, South - Inland Waters (FAO area)",
    "Asia - Inland Waters (FAO area)",
    "Atlantic, Eastern Central (FAO area)",
    "Atlantic, Northeast (FAO area)",
    "Atlantic, Southwest (FAO area)",
    "Atlantic, Western Central (FAO area)",
    "Europe - Inland Waters (FAO area)",
    "Standardized lowest detectable particle size (microns)",
    "Standardized sample size (number of fish)",
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
    "Fibres excluded",
    "Fibres not excluded",
    "Polymer ID not used",
    "Polymer ID used"
  )

Params1$parameter <- mapvalues(Params1$parameter,
                               from = levels(Params1$parameter),
                               to = gutmod_paramnames)

Params1$order <- c(nrow(Params1):1)

## HDPI plots
png('Gut Content HPDI Plot.png', 
    width = 14, 
    height = 15, 
    units = 'cm', 
    res = 500)

ggplot(HPD1) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed',
             size = 0.5,
             colour = pal[2]) +
  geom_errorbar(aes(x = param,
                    ymin = Q2.5,
                    ymax = Q97.5),
                size = 0.25) +
  geom_point(aes(x = param,
                 y = Estimate),
             size = 1,
             shape = 16,
             colour = pal[3]) +
  labs(x = 'Coefficient',
       y = '') +
  coord_flip() +
  theme1

dev.off()

## Posterior density plots

run2long <- extract.post(run2)

run2long$variable <- mapvalues(run2long$variable,
                               from = levels(run2long$variable),
                               to = gutmod_paramnames)

run2long$order <- c(nrow(run2long):1)

png(
  'Gut Content Posteriors Plot.png',
  width = 14,
  height = 19,
  units = 'cm',
  res = 500
)

ggplot(run2long) +
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

## Run again and estimate mu this time as well

gutdata.fitted <- cbind(gutdata, fitted(model1, 
                                        scale = 'response',
                                        re_formula = ~ (TL | region))/gutdata$N)

freshplot1 <-
  ggplot(subset(gutdata.fitted, study.habitat == 'Freshwater')) +
  geom_ribbon(aes(x = TL, ymin = Q2.5, ymax = Q97.5), 
              alpha = 0.75, fill = pal[3]) +
  geom_line(aes(x = TL, y = Estimate),
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
  ggplot(subset(gutdata.fitted, study.habitat == 'Marine')) +
  geom_ribbon(aes(x = TL, ymin = Q2.5, ymax = Q97.5), 
              alpha = 0.75, fill = pal[2]) +
  geom_line(aes(x = TL, y = Estimate),
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
  ggplot(subset(gutdata.fitted, study.habitat == 'Freshwater')) +
  geom_ribbon(aes(x = min.size, ymin = Q2.5, ymax = Q97.5), 
              alpha = 0.75, fill = pal[3]) +
  geom_line(aes(x = min.size, y = Estimate),
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

#### Predictions ####

## Simulate results according to trophic level

set.seed(54166)

gutdata.sim1 <-
  data.frame(TL = seq(from = 1.9, to = 5.1, length.out = 10000),
             min.size = rep(100, 10000),
             N = rep(50, 10000),
             polymer.ID = as.factor(rep(levels(gutdata$polymer.ID)[2], 10000)),
             blanks = as.factor(rep(levels(gutdata$blanks)[2], 10000)),
             environment = as.factor(rep(NA, 10000)),
             exclude.fib = as.factor(rep(levels(gutdata$exclude.fib)[2], 10000)),
             region = as.factor(sample(levels(gutdata$region), 
                                       10000, 
                                       replace = TRUE)))

gutdata.sim1 <- cbind(gutdata.sim1, 
                      predict(model1, 
                              newdata = gutdata.sim1,
                              probs = c(0.025, 0.125, 0.25, 0.375, 0.5,
                                        0.625, 0.75, 0.875, 0.975)))

gutdata.sim1[9:19] <- gutdata.sim1[9:19] / gutdata.sim1$N

png('MPs by Trophic Level Predictions Plot.png', width = 14, height = 16, 
    units = 'cm', res = 500)

ggplot() +
  geom_ribbon(
    data = gutdata.sim1,
    aes(x = TL,
        ymin = Q37.5,
        ymax = Q62.5),
    alpha = 0.75,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = gutdata.sim1,
    aes(x = TL,
        ymin = Q25,
        ymax = Q75),
    alpha = 0.5,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = gutdata.sim1,
    aes(x = TL,
        ymin = Q12.5,
        ymax = Q87.5),
    alpha = 0.25,
    fill = pal[1]
  ) +
  geom_ribbon(
    data = gutdata.sim1,
    aes(x = TL,
        ymin = Q2.5,
          ymax = Q97.5),
    alpha = 0.05,
    fill = pal[1]
  ) +
  geom_line(
    data = gutdata.sim1,
    aes(x = TL,
        y = Q50),
    colour = pal[5]) +
  geom_point(
    data = gutdata,
    aes(x = TL,
        y = Mpsgut),
    size = 0.5,
    colour = pal[5]) +
  facet_wrap( ~ region, ncol = 5,
              labeller = label_wrap_gen(width = 15)) +
  labs(x = 'Trophic Level',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       ))) +
  coord_cartesian(ylim = c(0, 1000)) +
  scale_y_continuous(trans = 'log1p',
                     breaks = c(0, 1, 10, 100, 500),
                     expand = c(0, 0.05)) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(1.9, 5.1)) +
  theme1

dev.off()

## Simulate results according to methodology

set.seed(63189)

gutdata.sim2 <-
  data.frame(TL = rep(3, 12000),
             min.size = seq(from = 0.5, to = 520, length.out = 12000),
             sample.size = rep(50, 12000),
             polymer.ID = as.factor(rep(2, 12000)),
             blanks = as.factor(rep(2, 12000)),
             environment = as.factor(rep(4, 12000)),
             exclude.fib = as.factor(sample(as.integer(gutdata$exclude.fib), 
                                            12000, 
                                            replace = TRUE)),
             region = as.factor(rep(16, 12000)))

gutdata.sim2$stand.TL <-
  (gutdata.sim2$TL - mean(gutdata$TL)) /
  sd(gutdata$TL - mean(gutdata$TL))

gutdata.sim2$stand.min.size <-
  (gutdata.sim2$min.size - mean(gutdata$min.size)) /
  sd(gutdata$min.size - mean(gutdata$min.size))

gutdata.sim2$stand.sample.size <-
  (gutdata.sim2$sample.size - mean(gutdata$N)) /
  sd(gutdata$N - mean(gutdata$N))

for (i in 1:12000) {
  mu <-
    run2$BUGSoutput$sims.list$alpha_region[, gutdata.sim2$region[i]] +
    run2$BUGSoutput$sims.list$beta_sample.size * 
    gutdata.sim2$stand.sample.size[i] +
    run2$BUGSoutput$sims.list$beta_min.size * gutdata.sim2$stand.min.size[i] +
    run2$BUGSoutput$sims.list$beta_TL[, gutdata.sim2$region[i]] *
    gutdata.sim2$stand.TL[i] +
    run2$BUGSoutput$sims.list$gamma_PID[, gutdata.sim2$polymer.ID[i]] +
    run2$BUGSoutput$sims.list$gamma_exclude.fib[, gutdata.sim2$exclude.fib[i]] +
    run2$BUGSoutput$sims.list$gamma_blanks[, gutdata.sim2$blanks[i]] +
    run2$BUGSoutput$sims.list$gamma_environment[, gutdata.sim2$environment[i]]
  sigma <- run2$BUGSoutput$sims.list$sigma
  p <- plogis(run2$BUGSoutput$sims.list$alpha_zero +
                run2$BUGSoutput$sims.list$beta_zero * 
                gutdata.sim2$stand.min.size)
  y <- rlnorm(12000, mu, sigm quantile(y, 0.75)
              gutdata.sim2$lower50[i] <- quantile(y, 0.25)
              gutdata.sim2$upper75[i] <- quantile(y, 0.875)
              gutdata.sim2$lower75[i] <- quantile(y, 0.125)
              gutdata.sim2$upper95[i] <- quantile(y, 0.975)
              gutdata.sim2$lower95[i] <- quantile(y, 0.025)
              gutdata.sim2$sample[i] <- sample(y, 1)
}

gutdata.sim2$exclude.fib <- mapvalues(gutdata.sim2$exclude.fib,
                                      from = levels(gutdata.sim2$exclude.fib),
                                      to = c("Fibres excluded",
                                             "Fibres not excluded"))

png('MPs Methodology Predictions Plot.png', width = 14, height = 8, 
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
  coord_cartesian(ylim = c(0, 25), xlim = c(0, 510)) +
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

#### Ingestion rate model ####

## Set up the data

ingestion <- subset(trophicfish2, !is.na(IR))

## Convert to successes/failures

ingestion$successes <- round(with(ingestion, N * IR), digits = 0)
ingestion$failures <- round(with(ingestion,  N * (1 - IR), digits = 0))

length(ingestion$species) # 639 data points
length(unique(ingestion$species)) # 478 species
length(unique(ingestion$family)) # from 157 families
length(unique(ingestion$study)) # 105 studies

ingestion$region <- as.character(ingestion$region)
ingestion$region <- as.factor(ingestion$region)

summary(ingestion)

## Specify model
ingmod <- function()
{
  # Likelihood
  for(i in 1:N)
  {
    y[i] ~ dbinom(p[i], n[i])
    logit(p[i]) <- alpha[region[i]] + beta[region[i]] * TL[i]
  }
  
  # Prior
  for(j in 1:Nregions)
  {
    alpha[j] ~ dnorm(mu_region, tau_region)
    beta[j] ~ dnorm(mu_TL, tau_TL)
  }
  mu_region ~ dnorm(0, 1)
  tau_region <- 1 / (sigma_region * sigma_region)
  sigma_region ~ dexp(1)
  mu_TL ~ dnorm(0, 1)
  tau_TL <- 1 / (sigma_TL * sigma_TL)
  sigma_TL ~ dexp(1)
}


## Initial values for MCMC chains
inginit <- function()
{
  list(
    "mu_region" = rnorm(1),
    "sigma_region" = 1,
    "mu_TL" = 1,
    "sigma_TL" = 1
  )
}

## Parameters to keep track of
ingparam <- c("alpha", "beta")


## Specify data
ingdata <- list(
  N = nrow(ingestion),
  Nregions = max(as.integer(ingestion$region)),
  y = as.numeric(ingestion$successes),
  n = as.numeric(ingestion$N),
  TL = as.numeric(scale(ingestion$TL, center = TRUE)),
  region = as.integer(ingestion$region)
)


## Run the model
ingrun1 <- jags.parallel(
  data = ingdata,
  inits = inginit,
  parameters.to.save = ingparam,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = ingmod
)

ingrun1
ingrun1mcmc <- as.mcmc(ingrun1)
xyplot(ingrun1mcmc, layout = c(6, ceiling(nvar(ingrun1mcmc)/6)))


# Update iterations to 10000
ingrun2 <- jags.parallel(
  data = ingdata,
  inits = inginit,
  parameters.to.save = ingparam,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 3,
  jags.seed = 123,
  model = ingmod
)

ingrun2
ingrun2mcmc <- as.mcmc(ingrun2)
xyplot(ingrun2mcmc, 
       layout = c(6, ceiling(nvar(ingrun2mcmc)/6)))

## Inference

ingrunHPD <- as.data.frame(summary(ingrun2mcmc)$quantiles)

ingrunMAP <- as.data.frame(summary(ingrun2mcmc)$statistics)

## Extract MAP and HPDIs for the parameters
Params2 <- data.frame(
  parameter = as.factor(rownames(ingrunMAP)),
  MAP = ingrunMAP[, 1],
  lower = ingrunHPD[, 1],
  upper = ingrunHPD[, 5]
)

with(ingestion, tapply(as.integer(region), region, mean))

Params2 <- Params2[c(1:36), ]
Params2$parameter <- as.character(Params2$parameter)
Params2$parameter <- as.factor(Params2$parameter)

ingrunparams <-
  c(
    "Africa - Inland Waters",
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
    "Atlantic, Southwest",
    "Atlantic, Western Central",
    "Europe - Inland Waters",
    "Indian Ocean, Antarctic",
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
    "Standardized trophic level:Indian Ocean, Antarctic"
  )

Params2$parameter <- mapvalues(
  Params2$parameter,
  from = levels(Params2$parameter),
  to = ingrunparams
)

Params2$sort <- c(nrow(Params2):1)

png('Ingestion HPDI Plot.png', 
    width = 14, 
    height = 12, 
    units = 'cm', 
    res = 500)

ggplot(Params2) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed',
             size = 0.5,
             colour = pal[2]) +
  geom_errorbar(aes(x = reorder(parameter, sort),
                    ymin = lower,
                    ymax = upper),
                size = 0.25) +
  geom_point(aes(x = reorder(parameter, sort),
                 y = MAP),
             size = 1.5,
             shape = 16,
             colour = pal[3]) +
  labs(x = 'Coefficient',
       y = '') +
  coord_flip() +
  theme1 +
  theme(plot.margin = margin(0.1, 0.5, 0.1, 0.5, unit = 'cm'))

dev.off()

## Posterior density plots

ingrun2long <- extract.post(ingrun2)

ingrun2long$variable <- mapvalues(ingrun2long$variable,
                                  from = levels(ingrun2long$variable),
                                  to = ingrunparams)

ingrun2long$order <- c(nrow(ingrun2long):1)

png(
  'Ingestion Posteriors Plot.png',
  width = 14,
  height = 14,
  units = 'cm',
  res = 500
)

ggplot(ingrun2long) +
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

## Rerun model to extract p
ingparam2 <- c("p")

ingrun3 <- jags.parallel(
  data = ingdata,
  inits = inginit,
  parameters.to.save = ingparam2,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 3,
  jags.seed = 123,
  model = ingmod
)

ingrun3mcmc <- as.mcmc(ingrun3)

ingestion$post.predict <-
  as.data.frame(summary(ingrun3mcmc)$statistics)[2:640, 1]
ingestion$lower95 <-
  as.data.frame(summary(ingrun3mcmc)$quantiles)[2:640, 1]
ingestion$upper95 <-
  as.data.frame(summary(ingrun3mcmc)$quantiles)[2:640, 5]

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


#### Body size MPs in guts model ####

## Set up the data
size <- subset(gutdata, total.length != 'NA')
size$life.stage <- as.factor(size$life.stage)

length(size$total.length)  # 394 data points remaining
length(unique(size$study))  # 61 studies
length(unique(size$species))  # 327 species
length(unique(size$study))  # 61 studies

## Specify model
sizemod <- function()
{
  # Likelihood
  for (i in 1:N)
  {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta_min.size * min.size[i] + beta_length * length[i]
  }
  # Prior
  alpha ~ dexp(1)
  beta_min.size ~ dnorm(-1, 1)
  beta_length ~ dnorm(0, 1)
  tau <- 1 / (sigma * sigma)
  sigma ~ dexp(1)
}

## Initial values
sizemod_init <- function()
{
  list(
    "alpha" = 1,
    "beta_min.size" = rnorm(1),
    "beta_length" = rnorm(1),
    "sigma" = 1
  )
}

## Parameters to keep track of
sizemod_params <- c("alpha", "beta_min.size", "beta_length", "sigma")

## Specify data

sizedata <- list(
  N = nrow(size),
  y = as.numeric(log(size$Mpsgut + 0.00538)),
  min.size = as.numeric(scale(size$min.size, center = TRUE)),
  length = as.numeric(scale(size$total.length, center = TRUE))
)

## Run the model
sizerun1 <- jags.parallel(
  data = sizedata,
  inits = sizemod_init,
  parameters.to.save = sizemod_params,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = sizemod
)

sizerun1
sizerun1mcmc <- as.mcmc(sizerun1)
xyplot(sizerun1mcmc)

## Increase to 5,000 interation

sizerun2 <- jags.parallel(
  data = sizedata,
  inits = sizemod_init,
  parameters.to.save = sizemod_params,
  n.chains = 3,
  n.iter = 5000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = sizemod
)

sizerun2  # DIC = 753
sizerun2mcmc <- as.mcmc(sizerun2)
summary(sizerun2mcmc)
xyplot(sizerun2mcmc)

## Inference

sizeHPDI <- as.data.frame(summary(sizerun2mcmc)$quantiles)

sizeMAP <- as.data.frame(summary(sizerun2mcmc)$statistics)

## Extract MAP and HPDIs for the parameters
sizeParams <- data.frame(
  parameter = as.factor(rownames(sizeMAP)),
  MAP = sizeMAP[, 1],
  lower = sizeHPDI[, 1],
  upper = sizeHPDI[, 5]
)

sizeParams <- sizeParams[-c(4:5), ]
sizeParams$parameter <- as.character(sizeParams$parameter)
sizeParams$parameter <- as.factor(sizeParams$parameter)

sizeParams$parameter <- mapvalues(
  sizeParams$parameter,
  from = levels(sizeParams$parameter),
  to = c(
    "Intercept",
    "Standardized Total Length (cm)",
    "Standardized lowest detectable particle size (microns)"
  )
)

png('Size Gut MPs HPDI Plot.png', 
    width = 9, 
    height = 2.25, 
    units = 'cm', 
    res = 500)

ggplot(sizeParams) +
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
             colour = pal[3]) +
  labs(x = 'Parameter',
       y = '') +
  coord_flip() +
  theme1

dev.off()

## Posterior density plots

sizerun2long <- extract.post(sizerun2)

sizerun2long$variable <- mapvalues(sizerun2long$variable,
                                   from = levels(sizerun2long$variable),
                                   to = c(
                                     "Intercept",
                                     "Standardized Total Length (cm)",
                                     "Standardized lowest detectable particle size (microns)"
                                   ))

sizerun2long$order <- c(nrow(sizerun2long):1)

png(
  'Body Size Model Posteriors Plot.png',
  width = 14,
  height = 5,
  units = 'cm',
  res = 500
)

ggplot(sizerun2long) +
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

## Rerun model and extract mu

sizemod_params1 <- c("mu")

sizerun3 <- jags.parallel(
  data = sizedata,
  inits = sizemod_init,
  parameters.to.save = sizemod_params1,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = sizemod
)

sizerun3
sizerun3mcmc <- as.mcmc(sizerun3)

size$post.predict <- 
  exp(as.data.frame(summary(sizerun3mcmc)$statistics)[-1, 1]) - 0.00538
size$lower95 <- 
  exp(as.data.frame(summary(sizerun3mcmc)$quantiles)[-1, 1]) - 0.00538
size$upper95 <- 
  exp(as.data.frame(summary(sizerun3mcmc)$quantiles)[-1, 5]) - 0.00538

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

## Simulate results for 1, 100, and 500 micron filter sizes

size.post <- data.frame(sizerun2$BUGSoutput$sims.list)

set.seed(5251)

generate.size.sim <- function(min.size, post, data) {
  size.sim <- 
    data.frame(length = seq(from = 1, to = 100, length.out = 5000),
               min.size = rep(min.size, 5000))
  
  size.sim$stand.length <-
    (size.sim$length - mean(data$total.length)) /
    sd(data$total.length - mean(data$total.length))
  
  size.sim$stand.min.size <-
    (size.sim$min.size - mean(data$min.size)) /
    sd(data$min.size - mean(data$min.size))
  
  for (i in 1:5000) {
    mu <-
      post$alpha + post$beta_length * size.sim$stand.length[i] +
      post$beta_min.size * size.sim$stand.min.size[i]
    y <- exp(rnorm(5000, mu, post$sigma)) - 0.00538
    size.sim$sigma[i] <- mean(post$sigma)
    size.sim$mu[i] <- mean(mu)
    size.sim$mean[i] <- exp(mean(mu)) - 0.00538
    size.sim$upper25[i] <- quantile(y, 0.625)
    size.sim$lower25[i] <- quantile(y, 0.375)
    size.sim$upper50[i] <- quantile(y, 0.75)
    size.sim$lower50[i] <- quantile(y, 0.25)
    size.sim$upper75[i] <- quantile(y, 0.875)
    size.sim$lower75[i] <- quantile(y, 0.125)
    size.sim$upper95[i] <- quantile(y, 0.975)
    size.sim$lower95[i] <- quantile(y, 0.025)
  }
  size.sim$sample <- with(size.sim, exp(rnorm(mu, mu, sigma)) - 0.00538)
  size.sim
}

sizeplot <- function(simdata){
  ggplot(simdata) +
    geom_point(aes(x = length,
                   y = sample),
               size = 0.25) +
    geom_ribbon(aes(x = length,
                    ymin = lower25,
                    ymax = upper25),
                alpha = 0.75,
                fill = pal[1]) +
    geom_ribbon(aes(x = length,
                    ymin = lower50,
                    ymax = upper50),
                alpha = 0.5,
                fill = pal[1]) +
    geom_ribbon(aes(x = length,
                    ymin = lower75,
                    ymax = upper75),
                alpha = 0.25,
                fill = pal[1]) +
    geom_ribbon(aes(x = length,
                    ymin = lower95,
                    ymax = upper95),
                alpha = 0.05,
                fill = pal[1]) +
    geom_line(aes(x = length,
                  y = mean),
              colour = pal[5]) +
    labs(x = 'Total Length (cm)',
         y = expression(paste(
           'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
         ))) +
    coord_cartesian(ylim = c(0, 1000)) +
    scale_y_continuous(trans = 'log1p',
                       expand = c(0, 0.05),
                       breaks = c(0, 1, 10, 100, 1000)) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, 100)) +
    theme1
}

size.sim1 <- generate.size.sim(1, size.post, size)

sizeA <- sizeplot(size.sim1)

size.sim2 <- generate.size.sim(100, size.post, size)

sizeB <- sizeplot(size.sim2)

size.sim3 <- generate.size.sim(500, size.post, size)

sizeC <- sizeplot(size.sim3)

png(
  'Size Effects Predictions Plot.png',
  width = 19,
  height = 7,
  units = 'cm',
  res = 500
)

plot_grid(sizeA, sizeB, sizeC, ncol = 3, labels = "AUTO")

dev.off()

#### Body size model for ingestion rates ####

## Set up the data

sizeing <- subset(allo, !is.na(IR))

## Convert to successes/failures

sizeing$successes <- round(with(sizeing, N * IR), digits = 0)
sizeing$failures <- round(with(sizeing,  N * (1 - IR), digits = 0))

length(sizeing$species) # 257 data points
length(unique(sizeing$species)) # 215 species
length(unique(sizeing$family)) # from 93 families
length(unique(sizeing$study)) # 47 studies

sizeing$region <- as.character(sizeing$region)
sizeing$region <- as.factor(sizeing$region)

summary(sizeing)

## Specify model
sizeingmod <- function()
{
  # Likelihood
  for(i in 1:N)
  {
    y[i] ~ dbinom(p[i], n[i])
    logit(p[i]) <- alpha[region[i]] + beta[region[i]] * length[i]
  }
  
  # Prior
  for(j in 1:Nregions)
  {
    alpha[j] ~ dnorm(mu_region, tau_region)
    beta[j] ~ dnorm(mu_length, tau_length)
  }
  mu_region ~ dnorm(0, 1)
  tau_region <- 1 / (sigma_region * sigma_region)
  sigma_region ~ dexp(1)
  mu_length ~ dnorm(0, 1)
  tau_length <- 1 / (sigma_length * sigma_length)
  sigma_length ~ dexp(1)
}


## Initial values for MCMC chains
sizeingmodinit <- function()
{
  list(
    "mu_region" = rnorm(1),
    "sigma_region" = 1,
    "mu_length" = rnorm(1),
    "sigma_length" = 1
  )
}

## Parameters to keep track of
sizeingmodparam <- c("alpha", "beta")


## Specify data
sizeingmoddata <- list(
  N = nrow(sizeing),
  Nregions = max(as.integer(sizeing$region)),
  y = as.numeric(sizeing$successes),
  n = as.numeric(sizeing$N),
  length = as.numeric(scale(sizeing$total.length, center = TRUE)),
  region = as.integer(sizeing$region)
)


## Run the model
sizeingrun1 <- jags.parallel(
  data = sizeingmoddata,
  inits = sizeingmodinit,
  parameters.to.save = sizeingmodparam,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = sizeingmod
)

sizeingrun1
sizeingrun1mcmc <- as.mcmc(sizeingrun1)
xyplot(sizeingrun1mcmc, layout = c(6, ceiling(nvar(sizeingrun1mcmc)/6)))

## Increase iterations to 10000

sizeingrun2 <- jags.parallel(
  data = sizeingmoddata,
  inits = sizeingmodinit,
  parameters.to.save = sizeingmodparam,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 3,
  jags.seed = 123,
  model = sizeingmod
)

sizeingrun2
sizeingrun2mcmc <- as.mcmc(sizeingrun2)
xyplot(sizeingrun2mcmc, layout = c(6, ceiling(nvar(sizeingrun2mcmc)/6)))

## Inference

sizeingHPDI <- as.data.frame(summary(sizeingrun2mcmc)$quantiles)

sizeingMAP <- as.data.frame(summary(sizeingrun2mcmc)$statistics)

## Extract MAP and HPDIs for the parameters
sizeingParams <- data.frame(
  parameter = as.factor(rownames(sizeingMAP)),
  MAP = sizeingMAP[, 1],
  lower = sizeingHPDI[, 1],
  upper = sizeingHPDI[, 5]
)

sizeingParams <- sizeingParams[-33, ]
sizeingParams$parameter <- as.character(sizeingParams$parameter)
sizeingParams$parameter <- as.factor(sizeingParams$parameter)

with(sizeing, tapply(as.integer(region), region, mean))

sizeing_paramNames <-
  c(
    "Africa - Inland Waters",
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
    "Atlantic, Southwest",
    "Atlantic, Western Central",
    "Europe - Inland Waters",
    "Indian Ocean, Antarctic",
    "Standardized total length:Africa - Inland Waters",
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
    "Standardized total length:Atlantic, Southwest",
    "Standardized total length:Atlantic, Western Central",
    "Standardized total length:Europe - Inland Waters",
    "Standardized total length:Indian Ocean, Antarctic"
  )

sizeingParams$parameter <- mapvalues(sizeingParams$parameter,
                                     from = levels(sizeingParams$parameter),
                                     to = sizeing_paramNames)

sizeingParams$sort <- c(nrow(sizeingParams):1)

png('Ingestion Body Size HPDI Plot.png', 
    width = 14, 
    height = 11, 
    units = 'cm', 
    res = 500)

ggplot(sizeingParams) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed',
             size = 0.5,
             colour = pal[2]) +
  geom_errorbar(aes(x = reorder(parameter, sort),
                    ymin = lower,
                    ymax = upper),
                size = 0.25) +
  geom_point(aes(x = reorder(parameter, sort),
                 y = MAP),
             size = 1,
             shape = 16,
             colour = pal[3]) +
  labs(x = 'Parameter',
       y = '') +
  coord_flip() +
  theme1 +
  theme(plot.margin = margin(0.1, 0.5, 0.1, 0.5, unit = 'cm'))

dev.off()

## Body size doesn't seem to affect ingestion
## Backs up similar patterns by region

## Posterior density plots

sizeingrun2long <- extract.post(sizeingrun2)

sizeingrun2long$variable <- mapvalues(sizeingrun2long$variable,
                                      from = levels(sizeingrun2long$variable),
                                      to = sizeing_paramNames)

sizeingrun2long$order <- c(nrow(sizeingrun2long):1)

png(
  'Body Size Ingestion Model Posteriors Plot.png',
  width = 14,
  height = 10,
  units = 'cm',
  res = 500
)

ggplot(sizeingrun2long) +
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

## Rerun model to extract p
sizeingparam2 <- c("p")

sizeingrun3 <- jags.parallel(
  data = sizeingmoddata,
  inits = sizeingmodinit,
  parameters.to.save = sizeingparam2,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 1000,
  n.thin = 3,
  jags.seed = 123,
  model = sizeingmod
)

sizeingrun3mcmc <- as.mcmc(sizeingrun3)

sizeing$post.predict <-
  as.data.frame(summary(sizeingrun3mcmc)$statistics)[2:258, 1]
sizeing$lower95 <-
  as.data.frame(summary(sizeingrun3mcmc)$quantiles)[2:258, 1]
sizeing$upper95 <-
  as.data.frame(summary(sizeingrun3mcmc)$quantiles)[2:258, 5]

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


#### Family model ####

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

## Model

## Specify model

fammod <- function()
{
  # Likelihood
  for (i in 1:N)
  {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <-
      alpha[family[i]] + beta_min.size * min.size[i] +
      beta_length[family[i]] * length[i]
  }
  
  # Prior
  for (j in 1:Nfamily)
  {
    alpha[j] ~ dnorm(mu_family, tau_family)
    beta_length[j] ~ dnorm(mu_length, tau_length)
  }
  sigma ~ dexp(1)
  tau <- 1 / (sigma * sigma)
  beta_min.size ~ dnorm(-1, 1)
  mu_family ~ dnorm(0, 1)
  sigma_family ~ dexp(1)
  tau_family <- 1 / (sigma_family * sigma_family)
  mu_length ~ dnorm(0, 1)
  sigma_length ~ dexp(1)
  tau_length <- 1 / (sigma_length * sigma_length)
}

## Generate initial values for MCMC

faminit <- function()
{
  list(
    "sigma" = 1,
    "beta_min.size" = rnorm(1),
    "mu_family" = rnorm(1),
    "sigma_family" = 1,
    "mu_length" = rnorm(1),
    "sigma_length" = 1
  )
}

## Keep track of parameters

famparam <- c("alpha", "beta_min.size", "beta_length", "sigma")

## Specify data

famdata <-
  list(
    y = log(fam$Mpsgut + 0.00538),
    length = as.numeric(scale(fam$total.length, center = TRUE)),
    min.size = as.numeric(scale(fam$min.size, center = TRUE)),
    family = as.integer(fam$family),
    N = nrow(fam),
    Nfamily = max(as.integer(fam$family))
  )

## Run the model
famrun1 <- jags.parallel(
  data = famdata,
  inits = faminit,
  parameters.to.save = famparam,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 1,
  jags.seed = 123,
  model = fammod
)

famrun1
famrun1mcmc <- as.mcmc(famrun1)
xyplot(famrun1mcmc, layout = c(6, ceiling(nvar(famrun1mcmc)/6)))

## Increase number of iterations to 10,000 

famrun2 <- jags.parallel(
  data = famdata,
  inits = faminit,
  parameters.to.save = famparam,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 2,
  jags.seed = 123,
  model = fammod
)

famrun2
famrun2mcmc <- as.mcmc(famrun2)
xyplot(famrun2mcmc, layout = c(6, ceiling(nvar(famrun1mcmc)/6)))

## Inference

famHPD <- as.data.frame(summary(famrun2mcmc)$quantiles)

famMAP <- as.data.frame(summary(famrun2mcmc)$statistics)

## Extract MAP and HPDIs for the parameters
fam.par.est <- data.frame(
  parameter = as.factor(rownames(famMAP)),
  MAP = famMAP[, 1],
  lower = famHPD[, 1],
  upper = famHPD[, 5]
)

with(fam, tapply(as.integer(family), family, mean))

fam.par.est <- fam.par.est[-c(26:27), ]
fam.par.est$parameter <- as.character(fam.par.est$parameter)
fam.par.est$parameter <- as.factor(fam.par.est$parameter)

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
    "Standardized lowest detectable particle size (microns)")


fam.par.est$parameter <-
  mapvalues(fam.par.est$parameter,
            from = levels(fam.par.est$parameter),
            to = famrun_paramNames)

png('Gut Content Family HPDI Plot.png', 
    width = 14, 
    height = 9, 
    units = 'cm', 
    res = 500)

ggplot(fam.par.est) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed',
             size = 0.5,
             colour = pal[2]) +
  geom_errorbar(aes(x = reorder(parameter, MAP),
                    ymin = lower,
                    ymax = upper),
                size = 0.25) +
  geom_point(aes(x = reorder(parameter, MAP),
                 y = MAP),
             size = 1,
             shape = 16,
             colour = pal[3]) +
  labs(x = 'Coefficient',
       y = '') +
  coord_flip() +
  theme1

dev.off()

## Posterior density plots

famrun2long <- extract.post(famrun2)

famrun2long$variable <- mapvalues(famrun2long$variable,
                                  from = levels(famrun2long$variable),
                                  to = famrun_paramNames)

famrun2long$order <- c(nrow(famrun2long):1)

png(
  'Family Model Posteriors Plot.png',
  width = 14,
  height = 10,
  units = 'cm',
  res = 500
)

ggplot(famrun2long) +
  geom_density_ridges(
    aes(x = value,
        y = reorder(variable, value, mean)),
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

## Run MCMC again and estimate mu

famparam2 <- "mu"

famrun3 <- jags.parallel(
  data = famdata,
  inits = faminit,
  parameters.to.save = famparam2,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 2,
  jags.seed = 123,
  model = fammod
)

famrun3mcmc <- as.mcmc(famrun3)

fam$post.predict <-
  exp(as.data.frame(summary(famrun4mcmc)$statistics)[2:181, 1]) - 0.00538
fam$lower95 <-
  exp(as.data.frame(summary(famrun4mcmc)$quantiles)[2:181, 1]) - 0.00538
fam$upper95 <-
  exp(as.data.frame(summary(famrun4mcmc)$quantiles)[2:181, 5]) - 0.00538

png('Gut Content Family Bayesian Plot.png', width = 9, height = 10, units = 'cm', 
    res = 500)

ggplot(fam) +
  geom_ribbon(aes(x = min.size, ymin = lower95, ymax = upper95), 
              alpha = 0.75, fill = pal[3]) +
  geom_line(aes(x = min.size, y = post.predict),
            size = 0.5, alpha = 0.8) +
  geom_point(aes(x = min.size, y = Mpsgut),
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

## Simulate predictions for different families if LDPS is held to 100 microns

set.seed(5251)

fam.sim <- 
  data.frame(family = as.factor(sample(as.integer(fam$family), 
                                       5000,
                                       replace = TRUE)),
             length = rep(20, 5000),
             min.size = rep(100, 5000))
fam.sim$stand.length <-
  (fam.sim$length - mean(fam$total.length)) /
  sd(fam$total.length - mean(fam$total.length))

fam.sim$stand.min.size <-
  (fam.sim$min.size - mean(fam$min.size)) /
  sd(fam$min.size - mean(fam$min.size))

for (i in 1:5000) {
  mu <-
    famrun2$BUGSoutput$sims.list$alpha[, fam.sim$family[i]] +
    (famrun2$BUGSoutput$sims.list$beta_length[, fam.sim$family[i]] *
       fam.sim$stand.length[i]) +
    (famrun2$BUGSoutput$sims.list$beta_min.size *
       fam.sim$stand.min.size[i])
  fam.sim$predict[i] <-
    exp(rnorm(1, mu, famrun2$BUGSoutput$sims.list$sigma)) - 0.00538
}

fam.sim$family <- as.factor(fam.sim$family)

fam.sim$family <- mapvalues(fam.sim$family,
                            from = levels(fam.sim$family),
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
  height = 8,
  units = 'cm',
  res = 500
)

ggplot() +
  geom_violin(
    data = fam.sim,
    aes(x = family,
        y = predict),
    fill = pal[3],
    alpha = 0.8,
    size = 0.5
  ) +
  geom_jitter(
    data = fam,
    aes(x = family,
        y = Mpsgut),
    size = 0.5,
    colour = pal[5],
    shape = 1
  ) +
  labs(x = 'Family',
       y = expression(paste(
         'Microplastic Concentration (particles ' ~ ind ^ -1 * ')'
       ))) +
  scale_y_continuous(trans = 'log1p', 
                     expand = c(0, 0), 
                     limits = c(-0.1, 325),
                     breaks = c(0, 1, 10, 100, 300)) +
  theme1 +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

dev.off()

#### Additional plots ####

## Plot according to lower limit of detection

png('Lower Limit Plot.png', width = 19, height = 10, 
    units = 'cm', res = 500)

ggplot(gutdata) +
  geom_point(aes(
    x = reorder(region, Mpsgut, mean),
    y = Mpsgut,
    size = N,
    colour = min.size
  ),
  alpha = 0.5) +
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
    breaks = c(0, 1, 5, 10, 20, 30),
    expand = c(0, 0.05),
    trans = 'log1p'
  ) +
  coord_flip() +
  scale_colour_gradient2(low = pal[3], mid = pal[2], high = pal[1],
                         midpoint = 250) +
  theme1

dev.off()

## Plot MP concentrations by region

png('MP Concentration by Region.png', width = 14, height = 10, 
    units = 'cm', res = 500)

ggplot(gutdata) + 
  geom_density_ridges(aes(x = Mpsgut,
                          y = reorder(region, Mpsgut, mean)),
                      fill = pal[2],
                      alpha = 0.8,
                      scale = 0.9) +
  geom_point(aes(x = Mpsgut,
                 y = region),
             size = 1, alpha = 0.3, colour = pal[5]) +
  labs(x = expression(paste(
    "Microplastic Concentration (particles '" ~
      ind ^ -1 * ")")),
    y = "FAO Area") +
  theme1

dev.off()

