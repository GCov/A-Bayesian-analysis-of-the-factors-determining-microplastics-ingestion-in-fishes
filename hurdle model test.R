library(R2jags)
library(lattice)
library(coda)
library(DHARMa)
library(ggplot2)

## Generate test data 

set.seed(15256)
x <- runif(1000, -2, 2)  # predictor affects both p zeros and outcome if not 0
pone <- plogis(1 + (2 * x + rnorm(1000, 0, 1)))  # p of not getting zero
plot(pone ~ x)  # visualize
zeros <- rbinom(1000, 1, pone)  # generate data, 0 for 0, 1 for not 0
y <- zeros
error <- rnorm(length(y[y == 1]), 0, 0.5)  # error for y
y[y == 1] <- exp(1.5 + (0.75 * x[y == 1] + error))  # generate y values when not 0
plot(y ~ x)  # visualize all data
     
## Specify model for JAGS

testmod <- function() {
  ## Likelihood of zero or not zero
  for(i in 1:N) {
    notzero[i] ~ dbinom(p[i], 1)
    logit(p[i]) <- alpha_zero + beta_zero * x[i]  # linear predictor
  }
  
  ## Log-normal likelihood when data is not zero
  for(j in 1:Ny) {
    y[j] ~ dlnorm(mu[j], tau)
    mu[j] <- alpha + beta * x2[j]
  }
  alpha_zero ~ dnorm(0, 1)
  beta_zero ~ dnorm(0, 1)
  alpha ~ dnorm(0, 1)
  beta ~ dnorm(0, 1)
  sigma ~ dexp(1)
  tau <- 1 / (sigma * sigma)
  
  ## Model predictions
  for(k in 1:N) {
    fittedzeros[k] ~ 
      dbinom(ilogit(alpha_zero + beta_zero * x[k]), 1)
    fittedall[k] ~ dlnorm((alpha + beta * x[k]), tau)
    fitted[k] <- fittedzeros[k] * fittedall[k]
  }
}

## Set initial values

testmodinit <- function() {
  list(
    alpha_zero = rnorm(1),
    beta_zero = rnorm(1),
    alpha = rnorm(1),
    beta = rnorm(1),
    sigma = rexp(1)
  )
}

## Parameters to keep

testmodparams <- c("alpha_zero",
                   "beta_zero",
                   "alpha",
                   "beta",
                   "sigma")
## Specify data

testmoddata <-
  list(
    N = length(x),
    Ny = length(y[y > 0]),
    notzero = zeros,
    y = y[y > 0],
    x = x,
    x2 = x[y > 0]
  )

## Run model (note that I tried this a few t imes and 10,000 iterations is best)

testmodrun1 <-
  jags.parallel(
    data = testmoddata,
    inits = testmodinit,
    parameters.to.save = testmodparams,
    n.iter = 10000,
    n.burnin = 1000,
    n.thin = 2,
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
    n.iter = 10000,
    n.burnin = 1000,
    n.thin = 2,
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

