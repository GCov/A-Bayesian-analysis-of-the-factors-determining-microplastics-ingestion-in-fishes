mod1 <-
  glmmTMB(
    log(Mpsgut + 1) ~ 
      scale(TL, center = TRUE) + 
      environment + 
      scale(min.size, center = TRUE) +
      polymer.ID +
      blanks +
      exclude.fib +
      (TL | region),
    weights = N,
    data = gutdata
  )

library(coda)
library(rjags)

jagsdata1 <- 
  list(MPs = log(gutdata$Mpsgut + 1),
       TL = as.numeric(scale(gutdata$TL, center = TRUE)),
       E = as.integer(gutdata$environment),
       LOD = as.numeric(scale(gutdata$min.size, center = TRUE)),
       ID = as.integer(gutdata$polymer.ID),
       B = as.integer(gutdata$blanks),
       EF = as.integer(gutdata$exclude.fib),
       R = as.integer(gutdata$region),
       SS = as.numeric(scale(gutdata$N, center = TRUE)))

cat(
  "model {
  
  # Likelihood
  for (i in 1:727) {
  MPs[i] ~ dnorm(mu[i], tau)
  mu[i] <- int + a * TL[i] + b[R[i]] + c[R[i]] * TL[i] + d[E[i]] + 
  e * LOD[i] + f[ID[i]] + g[B[i]] + h[EF[i]] + k * SS[i]
  }
  
  # Priors
  int ~ dexp(1)
  a ~ dnorm(0, 0.1)
  b ~ dnorm(0, 0.01)
  c ~ dnorm(0, 0.1)
  d ~ dnorm(0, 0.01)
  e ~ dnorm(0, 0.1)
  f ~ dnorm(0, 0.01)
  g ~ dnorm(0, 0.01)
  h ~ dnorm(0, 0.01)
  k ~ dnorm(0, 0.01)
  tau ~ dexp(0.1)
  }",
  file = "jagsmod1.txt"
)

jagsmod1 <- jags.model(file = "jagsmod1.txt",
                       data = jagsdata1,
                       n.chains = 4,
                       n.adapt = 500)
