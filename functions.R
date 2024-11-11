pdfsl <- function(z, lambda, link = 'logit'){
  if(link == 'logit'){
    pdfsl <- 2*dlogis(z, log = F)*pnorm(lambda*z)
  }else if(link == 'probit'){
    pdfsl <-2*dnorm(z, log = F)*pnorm(lambda*z)
  }else if(link == 'skewt'){
    pdfsl <- 2*dt(z, df = 5, log = F)*pnorm(lambda*z)
  }
  return(pdfsl)
}

cdfsl<-function(z,lambda, link = 'logit'){
  integrate(pdfsl,lower=-Inf, upper=z, lambda = lambda, link = link,
            stop.on.error=FALSE)$value
  }

F2<-function(lambda, link){
  f2<-function(z){
    if(link == 'logit'){
      f2s <- sqrt(2/pi)*z*dlogis(z, log = F)*exp(-.5*lambda^2*z^2)
    }else if(link == 'probit'){
      f2s <- sqrt(2/pi)*z*dnorm(z, log = F)*exp(-.5*lambda^2*z^2)
    }else  if(link == 'skewt'){
      f2s <- sqrt(2/pi)*z*dt(z, df = 5, log = F)*exp(-.5*lambda^2*z^2)
    }
    return(f2s)
  }
  return(f2)
}

pii.lambda <- function(lambda,eta, link = 'logit'){
  if(abs(eta)< 10){
    eta0 <- 10
    a <- integrate(F2(lambda, link),lower=0,upper=eta0)$value
    F3 <- ifelse(eta > 0,
                 integrate(F2(lambda, link),lower = 0,upper=eta)$value - a,
                 integrate(F2(lambda, link),lower=-eta0,upper=eta)$value)
  }else F3 <- 0
  return(F3)
}


F20<-function(link){
  f200<-function(z){
    if(link == 'logit'){
      f20 <- sqrt(2/pi)*z*dlogis(z, log = F)
    }else if(link=='probit'){
      f20 <-sqrt(2/pi)*z*dnorm(z, log = F)
    }else  if(link == 'skewt'){
      f20 <- sqrt(2/pi)*z*dt(z, df = 5, log = F)
    }
    return(f20)
  }
  return(f200)
}

pii.0<-function(eta, link){
  F30<-integrate(F20(link),lower=-Inf,upper=eta)
  return(F30$value)
}

###############################################################################
#Estimation
# Optim lambda know
iterbeta <- function(par, y, X, Ni, link = 'logit', lambda){
  roptim <- optim(par = par, y = y, X = X, Ni = Ni, link = 'logit',
                  lambda = lambda,
                  fn = negll.sl.groupb,method=c("BFGS"),
                  gr = negscore.betas.groupl,
                  control = list(maxit = 300))
  ll <- roptim$value#negll.sl.groupb(par = roptim$par, y, X, Ni, 'logit',lambda)
  return(ll)
}

fbeta <- function(par, y, X, Ni, link = 'logit'){
  la <- c(-10,-5,-2,-1,0,1,2,5,10)
  ll <- sapply(la, function(l) iterbeta(par, y, X, Ni, link, lambda = l))
  lambda <- la[which.min(ll)]
  return(lambda)
}

ll.sl.group <- function(par,y, X, Ni, link = 'logit'){
  lambda <- par[1]
  betas <- par[2:length(par)]
  eta <- X%*%betas
  pii <- sapply(eta,cdfsl,lambda, link)
  pii <- ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  val <-  sum(log(choose(Ni,y)) + y*log(pii)+(Ni-y)*log(1-pii))
  return(val)
}

negll.sl.group<-function(par,y, X, Ni, link = 'logit'){
  val <- -ll.sl.group(par, y, X, Ni, link)
  return(val)
}

ll.sl.groupl <- function(par,y, X, Ni, link = 'logit', beta0){
  lambda <- par
  betas  <- beta0
  eta    <- X%*%betas
  pii <- sapply(eta,cdfsl,lambda, link)
  pii <- ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  val <-  sum(log(choose(Ni,y)) + y*log(pii)+(Ni-y)*log(1-pii))
  return(val)
}

negll.sl.groupl<-function(par,y, X, Ni, link = 'logit',beta0){
  val <- -ll.sl.groupl(par,y, X, Ni, link, beta0)
  return(val)
}

ll.sl.groupb <- function(par,y, X, Ni, link = 'logit', lambda){
  betas  <- par
  eta    <- X%*%betas
  pii <- sapply(eta,cdfsl,lambda, link)
  pii <- ifelse(pii>=0.999,0.999, ifelse(pii<=0.001,0.001,pii))
  val <-  sum(log(choose(Ni,y)) + y*log(pii)+(Ni-y)*log(1-pii))
  return(val)
}

negll.sl.groupb<-function(par,y, X, Ni, link = 'logit', lambda){
  val <- -ll.sl.groupb(par,y, X, Ni, link, lambda)
  return(val)
}


#Score lambda vector
score.lambda.group <- function(par, y, X, Ni, link = 'logit') {
  lambda <- par[1]
  betas  <- par[2:length(par)]
  eta    <- X %*% betas
  pii <- sapply(eta, cdfsl, lambda, link)
  pii <- ifelse(pii >= 0.99, 0.99, ifelse(pii <= 0.01, 0.01, pii))
  pil <- sapply(eta, pii.lambda, lambda, link)
  mui <- as.vector(Ni*pii)
  vi <- pii * (1 - pii)
  si <- (pil / vi)
  Ulambda <- t(si)%*%(y-mui)
  return(Ulambda)
}

score.lambda.groupl <- function(par, y, X, Ni, link = 'logit', beta0) {
  lambda <- par[1]
  eta    <- X %*% beta0
  pii <- sapply(eta, cdfsl, lambda, link)
  pii <- ifelse(pii >= 0.999, 0.999, ifelse(pii <= 0.001, 0.001, pii))
  pil <- sapply(eta, pii.lambda, lambda, link)
  mui <- as.vector(Ni*pii)
  vi <- pii * (1 - pii)
  si <- (pil / vi)
  Ulambda <- t(si)%*%(y-mui)
  return(Ulambda)
}

negscore.lambda.group <- function(par, y, X, Ni, link = 'logit'){
  val <- -score.lambda.group(par, y, X, Ni, link)
  return(val)
}

negscore.lambda.groupl <- function(par, y, X, Ni, link = 'logit', beta0){
  val <- -score.lambda.groupl(par, y, X, Ni, link, beta0)
  return(val)
}

score.betas.group <- function(par,y, X, Ni, link = 'logit'){
  lambda <- par[1]
  betas  <- par[2:length(par)]
  eta <- X%*%betas
  pii <- sapply(eta, cdfsl, lambda, link)
  pii <- ifelse(pii >= 0.99, 0.99, ifelse(pii <= 0.01, 0.01, pii))
  fg <- sapply(eta, pdfsl, lambda, link)
  vi <- pii*(1-pii)
  wi<- as.vector(fg/vi)
  W <- diag(wi)
  mu <- as.vector(Ni*pii)
  Ub <- t(X)%*%W%*%(y-mu)
  return(Ub)
}

negscore.betas.group <- function(par,y, X, Ni, link = 'logit'){
  val <- -score.betas.group(par, y, X, Ni, link)
  return(val)
}


score.betas.groupl <- function(par,y, X, Ni, link = 'logit', lambda){
  betas  <- par
  eta <- X%*%betas
  pii <- sapply(eta, cdfsl, lambda, link)
  pii <- ifelse(pii >= 0.9999, 0.9999, ifelse(pii <= 0.0001, 0.0001, pii))
  fg <- sapply(eta, pdfsl, lambda, link)
  vi <- pii*(1-pii)
  wi<- as.vector(fg/vi)
  W <- diag(wi)
  mu <- as.vector(Ni*pii)
  Ub <- t(X)%*%W%*%(y-mu)
  return(as.vector(Ub))
}

negscore.betas.groupl <- function(par,y, X, Ni, link = 'logit', lambda){
  val <- -score.betas.groupl(par,y, X, Ni, link, lambda)
  return(val)
}

score.group <- function(par, y, X, Ni, link = 'logit'){
  s <- c(score.lambda.group(par, y, X, Ni, link),
         score.betas.group(par, y, X, Ni, link))
  return(s)
}

negscore.group <- function(par, y, X, Ni, link = 'logit'){
  val <- -score.group(par, y, X, Ni, link = 'logit')
  return(val)
}


FisherScoring <- function(par, y, X, max_iter = 200, tol = 1e-5){
  theta <- par
  for (iter in 1:max_iter) {
    Ul <- score.lambda(theta, y, X)
    Ubetas <- score.betas(theta, y, X)
    U <- c(Ul, Ubetas)
    K <- EIM(par, y)

    theta_new <- theta + solve(K) %*% U
    if (max(abs(theta_new - theta)) < tol) {
      cat("Convergencia alcanzada en", iter, "iteraciones\n")
      break
    }
    print(theta_new)
    theta <- theta_new
  }
  return(theta)
}

FisherScoring.group <- function(par, y, X, Ni, link = 'logit', max_iter = 500, tol = 1e-8){
  theta <- par
  for (iter in 1:max_iter) {
    U <- score.group(theta, y, X, Ni, link)
    K <- EIM.group(theta, y, X, Ni, link)

    theta_new <- theta + solve(K) %*% U

    if (max(abs(theta_new - theta)) < tol) {
      cat("Convergencia alcanzada en", iter, "iteraciones\n")
      break
    }
    print(theta_new[1])
    theta <- theta_new
  }
  return(theta)
}

FisherScoring.group.beta <- function(par, y, X, Ni, lambda, link = 'logit', max_iter = 500, tol = 1e-8){
  beta <- par
  for (iter in 1:max_iter) {
    U <- score.betas.group(c(lambda,beta), y, X, Ni, link)
    K <- Kbb.group(N, X, Ni, beta, lambda, link)
    beta_new <- beta + solve(K) %*% U
    if (max(abs(beta_new - beta)) < tol) {
      cat("Convergencia alcanzada en", iter, "iteraciones\n")
      break
    }
    print(beta_new[1])
    beta <- beta_new
  }
  return(beta)
}


Delta <- function(par, y, max_iter = 500, tol = 1e-3){
  N <- length(y)
  Ni <- rep(1,N)
  piil <- matrix(NA, N, 1)
  z <- numeric(N)
  pii <- matrix(NA, N, 1)
  dii <- numeric(N)
  theta <- theta.ini
  iter <- 1
  for (iter in 1:max_iter) {
    eta <- X%*%theta[2:length(theta)]
    summary(eta)
    pii <- apply(eta,1,function(e) cdf.sl(e, theta[1]))
    pii <- ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
    piil <- apply(eta, 1, function(e) F3(theta[1],e))

    zi <- y-Ni*pii
    dii <- Ni/(pii*(1-pii))
    Xtheta <- cbind(piil,X)
    D <- diag(dii)
    b <- lm(zi ~ -1 + Xtheta, weights = dii)
    b <- coef(b)
    theta_new <- theta + b
    print(theta_new)
    if (max(abs(theta_new - theta)) < tol) {
      cat("Convergencia alcanzada en", iter, "iteraciones\n")
      break
    }
    theta <- theta_new
  }
  print(iter)
  return(theta)
}



