library(R2jags)
library(lattice)
library(coda)
library(DHARMa)
library(ggplot2)

## Generate test data 

set.seed(15256)
x <- runif(1000, -2, 2)  # predictor
p <- plogis(1 + (0.5 * x))  # probability of structural zeros
plot(p ~ x)  # visualize

lambda <- exp(0.5 - (0.4 * x))  # generate lambda
plot(lambda ~ x)  # visualize

size <- runif(1000, 1, 2000)  # sample size

y <- (1 - rbinom(x, 1, p)) * rpois(x, lambda)

plot(y ~ x)  # visualize all data

## Specify model for JAGS

testmod <- function() {
  ## Likelihood
  for(i in 1:N) {
    y[i] ~ dpois(mu[i])
    mu[i] <- lambda[i] * (1 - z[i]) + 0.00001
    z[i] ~ dbern(p[i])
    logit(p[i]) <- alpha_p + beta_p * x[i]
    log(lambda[i]) <- alpha + beta * x[i]
  }
  
  ## Prior
  alpha_p ~ dnorm(0, 1)
  beta_p ~ dnorm(0, 1)
  alpha ~ dnorm(0, 1)
  beta ~ dnorm(0, 1)
  
  ## Model predictions
  for(j in 1:N) {
    fitted[j] ~ dpois(mu[j])
  }
}

## Set initial values

testmodinit <- function() {
  list(
    alpha_p = rnorm(1),
    beta_p = rnorm(1),
    alpha = rnorm(1),
    beta = rnorm(1)
  )
}

## Parameters to keep

testmodparams <- c("alpha_p",
                   "beta_p",
                   "alpha",
                   "beta")
## Specify data

testmoddata <-
  list(
    N = length(x),
    y = y,
    x = x
  )

## Run model (note that I tried this a few t imes and 10,000 iterations is best)

testmodrun1 <-
  jags.parallel(
    data = testmoddata,
    inits = testmodinit,
    parameters.to.save = testmodparams,
    n.iter = 2000,
    n.burnin = 1000,
    n.thin = 1,
    jags.seed = 1656,
    model = testmod
  )

testmodrun1.mcmc <- as.mcmc(testmodrun1)
xyplot(testmodrun1.mcmc)  # Check convergence
testmodrun1  # Output matched input parameters!!! (pretty much)

## Rerun to extract just model predictions

testmodparams1 <- c("fitted")

testmodrun2 <-
  jags.parallel(
    data = testmoddata,
    inits = testmodinit,
    parameters.to.save = testmodparams1,
    n.iter = 2000,
    n.burnin = 1000,
    n.thin = 1,
    jags.seed = 1656,
    model = testmod
  )

## Create DHARMa object

DHARMa.test <- createDHARMa(simulatedResponse = 
                              t(testmodrun2$BUGSoutput$sims.list$fitted),
                            observedResponse = y,
                            fittedPredictedResponse = 
                              apply(t(testmodrun2$BUGSoutput$sims.list$fitted),
                                    1,
                                    median),
                            integerResponse = F,
                            seed = 51561)  

plot(DHARMa.test) ## DHARMa says we have a good model

## Plot model predictions

predict <- apply(t(testmodrun2$BUGSoutput$sims.list$fitted),
                 1,
                 median)
lower95 <- apply(t(testmodrun2$BUGSoutput$sims.list$fitted),
                 1,
                 quantile,
                 probs = 0.025)
upper95 <- apply(t(testmodrun2$BUGSoutput$sims.list$fitted),
                 1,
                 quantile,
                 probs = 0.975)

ggplot() +
  geom_ribbon(aes(x = x,
                  ymin = lower95,
                  ymax = upper95),
              fill = 'red',
              alpha = 0.5) +
  geom_line(aes(x = x,
                y = predict)) +
  geom_point(aes(x = x,
                 y = y)) +
  theme_classic()

