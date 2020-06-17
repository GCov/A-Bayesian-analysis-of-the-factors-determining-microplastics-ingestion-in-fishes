library(coda)
library(R2jags)
library(mcmcplots)
library(arm)
  
## Specify model

jagsmod1 <- function()
{
  # Likelihood
  for (i in 1:N)
  {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(b, x[i,])
  }
  # Prior
  for (j in 1:nb)
  {
    b[j] ~ dnorm(0, 1)
  }
  sigma ~ dexp(2)
  tau <- 1/(sigma*sigma)
  # Predicted values
  for (i in 1:N)
  {
    fitted_mu[i] <- mu[i]
  }
}

## Generate initial values for MCMC

init1 <- function()
{
  list("b" = rnorm(ncol(x)),
       'sigma' = 1)
}

## Keep track of parameters

param1 <- c("b", "sigma")

## Build model design matrix

x <- 
  model.matrix(log(Mpsgut + 1) ~ scale(TL, center = TRUE)*region, 
               data = gutdata)

## Specify data

jagsdata1 <-
  list(
    y = log(gutdata$Mpsgut + 1),
    x = x,
    N = nrow(gutdata),
    nb = ncol(x)
  )

## Run the model
run1 <- jags(
  data = jagsdata1,
  inits = init1,
  parameters.to.save = param1,
  n.chains = 3,
  n.iter = 20000,
  n.burnin = 1000,
  model = jagsmod1
)

print(run1)
plot(run1)

run1mcmc <- as.mcmc(run1)
plot(run1mcmc)
run1mat <- as.matrix(run1mcmc)
run1dat <- as.data.frame(run1mat)

summary(run1mcmc) 

sim.MPs <- data.frame(TL = gutdata$TL)
sim.MPs$stan.TL <- scale(sim.MPs$TL, center = TRUE)

for(i in 1:nrows(x))
{
  for(j in 1:)
  mu <- run1$BUGSoutput$sims.list$b
  MPs <- rexp(lambda)
  sim.MPs$mean[1] <- mean(MPs)
  sim.MPs$lower95[1] <- quantile(MPs, 0.025)
  sim.MPs$upper95[1] <- quantile(MPs, 0.975)
  sim.MPs$lower50[1] <- quantile(MPs, 0.25)
  sim.MPs$upper50[1] <- quantile(MPs, 0.75)
}

ggplot(data = sim.MPs) +
  geom_ribbon(aes(x = TL, 
                  ymin = lower95, 
                  ymax = upper95), 
              fill = 'red', 
              alpha = 0.3) +
  geom_ribbon(aes(x = TL, 
                  ymin = lower50, 
                  ymax = upper50), 
              fill = 'red', 
              alpha = 0.3) +
  geom_line(aes(x = TL,
                y = mean),
            size = 1) +
  geom_point(aes(x = gutdata$TL,
                 y = gutdata$Mpsgut),
             size = 1) +
  labs(x = 'Trophic Level',
       y = 'MPs')
