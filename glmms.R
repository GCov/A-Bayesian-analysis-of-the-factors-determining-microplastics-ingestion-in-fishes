library(glmmTMB)

glmm1 <- 
  glmmTMB(totalcount ~ 
            offset(log(N)) +
            scale(min.size, center = TRUE) - 1 + 
            (scale(TL, center = TRUE) | region) +
            polymer.ID +
            exclude.fib +
            blanks,
          data = gutdata,
          family = nbinom1(link = "log"))
summary(glmm1)

res1 <- simulateResiduals(glmm1)
plot(res1)

glmm2 <-
  glmmTMB(cbind(successes, failures) ~ 
            scale(min.size, center = TRUE) -1 +
            (TL | region),
          family = betabinomial(link = 'logit'),
          data = ingestion)
summary(glmm2)

glmm3 <-
  glmmTMB(cbind(successes, failures) ~ 
            scale(min.size, center = TRUE) -1 +
            exclude.fib +
            (TL | region),
          family = binomial(link = 'logit'),
          data = ingestion)
summary(glmm3)
anova(glmm2, glmm3)

res2 <- simulateResiduals(glmm3)
plot(res2)
testDispersion(res2)
testZeroInflation(res2)

glmm4 <- 
  glmmTMB(cbind(successes, failures) ~ 
            scale(min.size, center = TRUE) -1 +
            (TL | region),
          family = binomial(link = 'logit'),
          ziformula = ~ scale(min.size, center = TRUE),
          data = ingestion)
summary(glmm4)

anova(glmm2, glmm3, glmm4)

glmm5 <-
  glmmTMB(totalcount ~ 
            offset(log(N)) +
            exclude.fib +
            scale(min.size, center = TRUE) - 1 + 
            (scale(total.length, center = TRUE) | region),
          data = size,
          family = nbinom1(link = "log"),
          ziformula = ~ scale(min.size, center = TRUE))
summary(glmm5)

glmm6 <-
  glmmTMB(totalcount ~ 
            offset(log(N)) +
            exclude.fib +
            scale(min.size, center = TRUE) - 1 + 
            (scale(total.length, center = TRUE) | region),
          data = size,
          family = nbinom1(link = "log"))
summary(glmm6)
anova(glmm5, glmm6)

glmm7 <- 
  glmmTMB(totalcount ~ 
            offset(log(N)) +
            exclude.fib +
            scale(min.size, center = TRUE) - 1 + 
            (scale(total.length, center = TRUE) | family),
          data = fam,
          family = poisson(link = "log"))
summary(glmm7)
res7 <- simulateResiduals(glmm7)
plot(res7)
