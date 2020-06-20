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
    mu[i] <- inprod(b, x[i,]) + inprod(c, z[i,])
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
  sigma ~ dexp(2)
  tau <- 1/(sigma*sigma)
  mu_b ~ dnorm(0, 1)
  sigma_b ~ dexp(2)
  tau_b <- 1/(sigma_b*sigma_b)
}

## Generate initial values for MCMC

init1 <- function()
{
  list("b" = rnorm(ncol(x[ , c(1:20, 38:55)])),
       "c" = rnorm(ncol(x[ , c(21:37)])),
       'sigma' = 1,
       "mu_b" = 0,
       "sigma_b" = 1)
}

## Keep track of parameters

param1 <- c("b", "c", "sigma", "mu", "mu_b", "sigma_b")

## Build model design matrix

x <- model.matrix(log(Mpsgut + 1) ~ 
                    scale(TL, center = TRUE) * region +
                    environment +
                    scale(min.size, center = TRUE) +
                    polymer.ID +
                    blanks +
                    exclude.fib +
                    N,
                  data = gutdata)

## Specify data

jagsdata1 <-
  list(
    y = log(gutdata$Mpsgut + 1),
    x = x[ , c(1:20, 38:55)],
    z = x[ , c(21:37)],
    N = nrow(gutdata),
    Nb = ncol(x[ , c(1:20, 38:55)]),
    Nc = ncol(x[ , c(21:37)])
  )

## Run the model
run1 <- jags(
  data = jagsdata1,
  inits = init1,
  parameters.to.save = param1,
  n.chains = 3,
  n.iter = 30000,
  n.burnin = 1000,
  model = jagsmod1
)

print(run1)
plot(run1)

run1mcmc <- as.mcmc(run1)
summary(run1mcmc) 
HPD1 <- as.data.frame(summary(run1mcmc)$quantiles)
MAP1 <- as.data.frame(summary(run1mcmc)$statistics)

mu1 <- as.data.frame(run1$BUGSoutput$sims.list$mu)
sigma1 <- as.data.frame(run1$BUGSoutput$sims.list$sigma)

## Extract MAP and HPDIs for the parameters
Params1 <- data.frame(parameter = colnames(x))
Params1$parameter <- as.factor(Params1$parameter)
Params1$MAP <- NA
Params1$lower <- NA
Params1$upper <- NA

Params1[c(1:20, 38:55), 2] <- MAP1[1:38, 1]
Params1[21:37, 2] <- MAP1[39:55, 1]
Params1[c(1:20, 38:55), c(3, 4)] <- HPD1[1:38, c(1, 5)]
Params1[21:37, c(3,4)] <- HPD1[39:55, c(1, 5)]

Params1$parameter <- mapvalues(
  Params1$parameter,
  from = levels(Params1$parameter),
  to = c(
    "Intercept",
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
  scale_fill_manual(values = pal[1, 4],
                    name = 'Fibres Excluded?') +
  scale_colour_manual(values = pal[1, 4],
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
  scale_fill_manual(values = pal[1, 4],
                    name = 'Fibres Excluded?') +
  scale_colour_manual(values = pal[1, 4],
                      name = 'Fibres Excluded?') +
  scale_x_continuous(breaks = seq(from = 2, to = 5, by = 1)) +
  scale_y_continuous(trans = 'log1p', breaks = c(0, 1, 10, 30)) +
  theme1


png('Gut Content Bayesian Plot.png', width = 16.7, height = 20, units = 'cm', res = 300)

plot_grid(freshplot, marineplot, labels = c('A', 'B'), rel_heights = c(1,1.9),
          nrow = 2, align = 'v')

dev.off()
