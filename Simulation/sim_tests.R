###############################################################################
#Monte Carlo Simulation Studies: Simulation Study III
###############################################################################

#Time
start_time <- proc.time()

sim <- function(x,N, p, lambda.ini, X, alpha.cov, link){
  source('functions.R')

  betas<-matrix(1,p,1)
  eta <- X%*%betas
  prob <- sapply(eta, cdfsl, lambda.ini, link = link) #all.equal F<-exp(eta)/(1+exp(eta)) when lambda=0
  prob <- ifelse(prob >= 0.9999, 0.9999, ifelse(prob<=0.0001,0.0001,prob))

  #Generate replicates of response
  R <- 6#500 #number of replicates
  R.aux <- R*4#1000 #number of additional replicates
  r <- 1 #index of replicates
  Y <- matrix(NA,N,R.aux) #matrix of replicates for the response
  pl <- dim(X)[2] #dimension of regression coefficients vector of logistic model
  EPL <- matrix(NA,pl,R.aux) #matrix of parameter estimates of logistic model
  psl <- pl+1 #dimension of regression coefficients vector of skew-logistic model
  EPSL<-matrix(NA,psl,R)
  RL <- numeric(R)
  W <- numeric(R)
  pvalueRV <- numeric(R)
  pvalueW<-numeric(R)

  count<-0
  count.error<-0
  repeat{
    if (count == R.aux) break
    count <- count+1
    y <- rbinom(N,1,prob) #simulating Bernoulli random variables of skew-logistic model
    options(warn = 2) #to transform the warnings in error
    df <- data.frame(y,X)
    fit.l <- try (glm(y ~ X-1, family = binomial(link = 'logit'), data=df,
                      control = glm.control(maxit = 500)), TRUE)
    if(inherits(fit.l,"try-error"))
    {count.error <- count.error+1
    next
    }
    else
      Y[,count] <- y
    if(link == 'logit'){
      EPL[,count] <- fit.l$coefficients
    }else if(link == 'probit'){
      EPL[,count] <- fit.l$coefficients/1.7
    }else if(link == 'skewt'){
      EPL[,count] <- fit.l$coefficients/1.52
    }
  }

  Y <- Y[,!colSums(is.na(Y)) == nrow(Y), drop=FALSE]
  EPL <- EPL[,!colSums(is.na(EPL)) == nrow(EPL), drop=FALSE]
  Y <- Y[,1:R]
  EPL <- EPL[,1:R]
  Ni <- rep(1,N)

  for(r in 1:R){
    beta0 <- EPL[,r]
    theta0 <- c(0, beta0)
    fit.sl <- optimx(par = theta0,  y = Y[,r], X = X, Ni = Ni, link = link,
                     fn = negll.sl.group,
                     method = c('BFGS'),
                     control = list(maxit = 500))

    coef.opt<-coef(fit.sl)
    EPSL[1,r]<-coef.opt[,1]#lambda
    EPSL[2:psl,r]<-coef.opt[,2:psl]#betas
    #likelihood ratio test
    lambdaH0 <- 0
    theta.est.H0 <-c(lambdaH0, EPL[,r])
    l0 <- ll.sl.group(theta.est.H0, y = Y[,r], X, Ni, link)
    theta.est.H1 <- EPSL[,r]
    l1 <- ll.sl.group(theta.est.H1, y = Y[,r], X, Ni, link)
    RL[r] <- -2 * (l0 - l1)
    pvalueRV[r] <- pchisq(RL[r], 1, lower.tail = FALSE)

    #Wald test
    hessian.sl <- attributes(fit.sl)$details["BFGS", "nhatend"][[1]]
    oim.sl <- hessian.sl
    var.error <- diag(solve(oim.sl))
    W[r] <- (EPSL[1,r]-lambdaH0)^2/var.error[1]
    pvalueW[r] <- pchisq(W[r], 1, lower.tail=FALSE)
  }

  datalist <- list("N" = N, "p" = p, "R"=R, "lambda.ini" = lambda.ini,
                   "alpha.cov" = alpha.cov, "Y" = Y, "X" = X,
                   "EPL" = EPL, "count" = count, "count.error" = count.error,
                   "EPSL" = EPSL, 'RL' = RL, 'pvalueRV' = pvalueRV, 'W' = W, 'pvalueW' = pvalueW)

  return(datalist)
}


proc.time() - start_time