library(coda)
library(rjags)

TM1matrix <- with(gutdata, model.matrix(Mpsgut ~ TL + region))

jagsdata1 <- 
  list(MPs = gutdata$Mpsgut,
       TL = as.numeric(scale(gutdata$TL, center = TRUE)),
       N = length(gutdata$Mpsgut),
       region = as.integer(gutdata$region))

cat(
  "model {
  
  # Likelihood
  for (i in 1:N) {
  MPs[i] ~ dexp(lambda[i])
  log(lambda[i]) <- int + a * TL[i] + b[region[i]]
  }
  
  # Priors
  int ~ dexp(0.1)  
  a ~ dnorm(0, 1)
  b[region] ~ dnorm(0, 2.3)
  }",
  file = "jagsmod1.txt"
)

params1 <- c("int", "a", "b")

jagsmod1 <- jags.model(file = "jagsmod1.txt",
                       data = jagsdata1,
                       n.chains = 4,
                       n.adapt = 500)

samps1 <- coda.samples(jagsmod1, params1, n.iter = 10000)

summary(samps1)

summary(window(samps1, start = 5001))

plot(samps1)

a <- as.numeric(samps1[[1]][, 1])
int <- as.numeric(samps1[[1]][, 2])

sim.MPs <- data.frame(TL = gutdata$TL)
sim.MPs$stan.TL <- scale(sim.MPs$TL, center = TRUE)

for(i in 1:length(sim.MPs$TL)) {
  lambda <- int + a*sim.MPs$stan.TL[i]
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
