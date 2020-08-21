library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(brms)
library(bayesplot)
library(DHARMa)
color_scheme_set("viridis")

model1 <- brm(bf(totalcount ~ 
                   0 + 
                   offset(log(N)) + 
                   scale(TL, center = TRUE) + 
                   scale(min.size, center = TRUE) +
                   polymer.ID + blanks + exclude.fib +
                   (1 + TL | region) + 
                   (1 | environment),
                 zi ~ scale(min.size, center = TRUE)),
              data = gutdata,
              family = zero_inflated_poisson(link = "log", 
                                             link_zi = "logit"),
              prior = c(set_prior("normal(0, 1)", 
                                  class = "b"),
                        set_prior("normal(0, 1)",
                                  class = "sd"),
                        set_prior("normal(0, 1)",
                                  class = "b",
                                  dpar = "zi")),
              warmup = 500, 
              iter = 2000,
              chains = 3,
              inits = "random",
              cores = 8,
              seed = 8913,
              control = list(max_treedepth = 15))

mcmc_trace(model1)

summary(model1)

response1 <- apply(t(posterior_predict(model1, draws = 1000)), 
                   2,
                   function(x) {
                     x / gutdata$N
                   })
fitted1 <- apply(t(posterior_predict(model1, draws = 1000)), 1, median)

diagnose1 <- createDHARMa(simulatedResponse = response1,
                          observedResponse = gutdata$Mpsgut,
                          fittedPredictedResponse = fitted1)
plot(diagnose1)

conditional_effects(model1)

mcmc_plot(model1)

