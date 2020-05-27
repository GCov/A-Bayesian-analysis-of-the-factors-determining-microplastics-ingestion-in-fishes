library(MCMCglmm)

### MCMCglmm glm stuff

gutdata <- as.data.frame(gutdata)

prior.MCMCmod1 <-
  list(R=list(V=1,fix=1),
       G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

MCMCmod1 <-
  MCMCglmm(
    log(Mpsgut + 1) ~ TL * study.habitat,
    random = ~ (TL | region),
    rcov = ~ N,
    family = 'gaussian',
    prior = prior.MCMCmod1,
    data = gutdata,
    nitt = 1000,
    thin = 10,
    burnin = 100
  )
