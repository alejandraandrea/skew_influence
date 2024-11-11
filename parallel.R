rm(list=ls())

library(parallel)
source("sim.R")

link <- 'logit'
ncores <- 4
set.seed(3713)
NN <- c(100,200,400)
pp <- c(1,5,10,15)

lambda.ini <- .5
alpha.cov <- 0.1#0.10 0.50, 0.90
xcov <- function(x){(1-alpha.cov^2)^(1/2)*x + alpha.cov*zq}

for(i in NN){
  N <- i
  zq  <- rnorm(N,0,1)
  for(j in pp){
    p <- j
    mnormal <- matrix(rnorm(N*p, mean=0,sd=1), N, p)
    mxcov <- apply(mnormal,2,xcov)
    X <- mxcov
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {library(optimx);library(PresenceAbsence);library(matrixcalc)})
    clusterExport(cl, varlist = c('sim','N','p','lambda.ini','X','alpha.cov','link'))
    sim <- parLapply(cl, 1:10, function(x) sim(x, N, p, lambda.ini , X, alpha.cov, link))
    stopCluster(cl)
    saveRDS(sim, file = paste0('n',N,'_p',p,'_',link,'.rdata'))
  }
}

