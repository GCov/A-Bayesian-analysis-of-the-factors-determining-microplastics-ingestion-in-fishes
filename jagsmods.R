library(coda)
library(R2jags)
library(mcmcplots)
library(arm)
library(ggplot2)
library(cowplot)
library(plyr)

pal <- c("#c52f01",  #rust
         "#b88100",  #dark goldenrod 
         "#0076a9",  #celadon blue
         "#5a8400",  #avocado
         "#011e44")  #oxford blue

## Specify model

jagsmod1 <- function()
{
  # Likelihood
  for (i in 1:N)
  {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + inprod(b, x[i, ]) + inprod(c, z[i, ])
  }
  for (j in 1:Nb)
  {
    b[j] ~ dnorm(mu_b, tau_b)
  }
  # Prior
  for (k in 1:Nc)
  {
    c[k] ~ dnorm(0, 1)
  }
  alpha ~ dexp(1)
  sigma ~ dexp(1)
  tau <- 1 / (sigma * sigma)
  mu_b ~ dnorm(0, 1)
  sigma_b ~ dexp(1)
  tau_b <- 1 / (sigma_b * sigma_b)
}

## Generate initial values for MCMC

init1 <- function()
{
  list(
    "alpha" = 1,
    "b" = rnorm(ncol(x[, c(2:19, 37:54)])),
    "c" = rnorm(ncol(x[, c(2, 20:36)])),
    'sigma' = 1,
    "mu_b" = 0,
    "sigma_b" = 1
  )
}

## Keep track of parameters

param1 <- c("alpha", "b", "c", "sigma", "mu_b", "sigma_b", "mu")

## Build model design matrix

x <- model.matrix(
  log(Mpsgut + 1) ~
    scale(TL, center = TRUE) * region +
    environment +
    scale(min.size, center = TRUE) +
    polymer.ID +
    blanks +
    exclude.fib +
    N,
  data = gutdata
)[,-1]

## Specify data

jagsdata1 <-
  list(
    y = log(gutdata$Mpsgut + 1),
    x = x[, c(2:19, 37:54)],
    z = x[, c(2, 20:36)],
    N = nrow(gutdata),
    Nb = ncol(x[, c(2:19, 37:54)]),
    Nc = ncol(x[, c(2, 20:36)])
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

run1mcmc <- as.mcmc(run1)
traceplot(run1mcmc)

## Extend burnin to 2000 and up number of iterations to 10,000

run2 <- jags(
  data = jagsdata1,
  inits = init1,
  parameters.to.save = param1,
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 2000,
  n.thin = 8,
  jags.seed = 123,
  model = jagsmod1
)

run2mcmc <- as.mcmc(run2)
traceplot(run2mcmc)

## Increase iteration to 50000

run3 <- jags(
  data = jagsdata1,
  inits = init1,
  parameters.to.save = param1,
  n.chains = 3,
  n.iter = 50000,
  n.burnin = 2000,
  n.thin = 24,
  jags.seed = 123,
  model = jagsmod1
)

run3mcmc <- as.mcmc(run3)
traceplot(run3mcmc)


## Inference

summary(run3mcmc) 
HPD1 <- as.data.frame(summary(run3mcmc)$quantiles)
MAP1 <- as.data.frame(summary(run3mcmc)$statistics)

mu1 <- as.data.frame(run3$BUGSoutput$sims.list$mu)
sigma1 <- as.data.frame(run3$BUGSoutput$sims.list$sigma)

## Extract MAP and HPDIs for the parameters
Params1 <- data.frame(parameter = colnames(x), stringsAsFactors = FALSE)
Params1$MAP <- NA
Params1$lower <- NA
Params1$upper <- NA

Params1[c(2:19, 37:54), 2] <- MAP1[2:37, 1]
Params1[c(1, 20:36), 2] <- MAP1[38:55, 1]
Params1[c(2:19, 37:54), c(3, 4)] <- HPD1[2:37, c(1, 5)]
Params1[c(1, 20:36), c(3,4)] <- HPD1[38:55, c(1, 5)]
Params1[55, 2] <- MAP1[1, 1]
Params1[55, c(3, 4)] <- HPD1[1, c(1,5)]
Params1[55, 1] <- "Intercept"
Params1$parameter <- as.factor(Params1$parameter)

Params1$parameter <- mapvalues(
  Params1$parameter,
  from = levels(Params1$parameter),
  to = c(
    "Blanks were used",
    "Bathydemersal (environment)",
    "Bathypelagic (environment)",
    "Benthopelagic (environment",
    "Demersal (environment)",
    "Freshwater benthopelagic (environment)",
    "Freshwater demersal (environment)",
    "Freshwater pelagic (environment)",
    "Freshwater pelagic-neritic (environment",
    "Pelagic (environment)",
    "Pelagic-neritic (environment",
    "Pelagic-oceanic (environment)",
    "Reef-associated (environment)",
    "Fibres excluded",
    "Intercept",
    "Standardized sample size (number of fish)",
    "Polymer ID used",
    "America, North - Inland Waters (region)",
    "America, South - Inland Waters (region)",
    "Asia - Inland Waters (region)",
    "Atlantic, Eastern Central (region)",
    "Atlantic, Northeast (region)",
    "Atlantic, Southwest (region)",
    "Atlantic, Western Central (region)",
    "Europe - Inland Waters (region)",
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
    "Standardized lower limit of detection (microns)",
    "Standardized trophic level",
    "Standardized trophic level:America, North - Inland Waters",
    "Standardized trophic level:America, South - Inland Waters",
    "Standardized trophic level:Asia - Inland Waters",
    "Standardized trophic level:Atlantic, Eastern Central",
    "Standardized trophic level:Atlantic, Northeast",
    "Standardized trophic level:Atlantic, Southwest",
    "Standardized trophic level:Atlantic, Western Central",
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
    "Standardized trophic level:Pacific, Western Central"
  )
)

png('Gut Content HPDI Plot.png', 
    width = 13, 
    height = 15, 
    units = 'cm', 
    res = 300)

ggplot(Params1) +
  geom_errorbar(aes(x = parameter,
                    ymin = lower,
                    ymax = upper),
                size = 0.5) +
  geom_point(aes(x = parameter,
                 y = MAP),
             size = 1,
             shape = 19) +
  geom_hline(aes(yintercept = 0),
             linetype = 'dashed',
             size = 0.5) +
  labs(x = 'Coefficient',
       y = '') +
  coord_flip() +
  theme1

dev.off()

# set.seed(1234)
# for(j in 1:nrow(gutdata)) {
#   for(i in 1:nrow(mu1)) {
#   predict[i] <- rnorm(1, mean = mu1[i, j], sd = sigma1[i, ])
#   }
#   gutdata$post.predict[j] <- exp(mean(predict)) - 1
#   gutdata$lower95[j] <- exp(quantile(predict, 0.025)) - 1
#   gutdata$upper95[j] <- exp(quantile(predict, 0.975)) - 1
# }

set.seed(1234)
for(j in 1:nrow(gutdata)) {
  for(i in 1:nrow(mu1)) {
    mu[i] <- mu1[i, j]
  }
  gutdata$post.predict[j] <- exp(mean(mu)) - 1
  gutdata$lower95[j] <- exp(quantile(mu, 0.025)) - 1
  gutdata$upper95[j] <- exp(quantile(mu, 0.975)) - 1
}

# gutdata$post.predict[gutdata$post.predict < 0] <- 0
# gutdata$lower95[gutdata$lower95 < 0] <- 0


freshplot <-
  ggplot(subset(gutdata, study.habitat == 'Freshwater')) +
  geom_line(aes(x = TL, y = post.predict, colour = exclude.fib),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95, fill = exclude.fib), 
              alpha = 0.3) +
  geom_point(aes(x = TL, y = Mpsgut, colour = exclude.fib),
             shape = 1, size = 0.75) +
  facet_wrap(~ region, scales = 'free_x', ncol = 4,
             labeller = label_wrap_gen(width = 20)) +
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

marineplot <-
  ggplot(subset(gutdata, study.habitat == 'Marine')) +
  geom_line(aes(x = TL, y = post.predict, colour = exclude.fib),
            size = 0.5, alpha = 0.8) +
  geom_ribbon(aes(x = TL, ymin = lower95, ymax = upper95, fill = exclude.fib), 
              alpha = 0.3) +
  geom_point(aes(x = TL, y = Mpsgut, colour = exclude.fib),
             shape = 1, size = 0.75) +
  facet_wrap(~ region, scales = 'free_x', ncol = 4,
             labeller = label_wrap_gen(width = 20)) +
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


png('Gut Content Bayesian Plot.png', width = 16.7, height = 20, units = 'cm', res = 300)

plot_grid(freshplot, marineplot, labels = c('A', 'B'), rel_heights = c(1,1.9),
          nrow = 2, align = 'v')

dev.off()
