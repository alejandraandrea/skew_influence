source('functions.R')

Kll.group <- function(N,X,Ni,betas,lambda,link){
  eta <- X%*%betas
  pii <- sapply(eta, cdfsl, lambda, link)
  pii <- ifelse(pii >= 0.9999, 0.9999, ifelse(pii <= 0.0001, 0.0001, pii))
  pil <- sapply(eta, pii.lambda, lambda, link)
  vi <- pii*(1-pii)
  s <- sum(Ni*pil^2/vi)
  return(s)
}

K00.group <- function(N,X,Ni,betas, link){
  eta <- X%*%betas
  pi0 <- sapply(eta, pii.0, link)
  if(link == 'logit'){
    Gi <- plogis(eta,log=FALSE)
  }else if(link == 'probit'){
    Gi <- pnorm(eta,log=FALSE)
  }else Gi <- pt(eta,df = 5, log=FALSE)

  vi <- Gi*(1-Gi)
  s <- sum(Ni*pi0^2/vi)
  return(s)
}

Klb.group<-function(N,X,Ni,betas,lambda, link){
  eta <- X%*%betas
  pii <- sapply(eta, cdfsl, lambda, link)
  pii <- ifelse(pii >= 0.9999, 0.9999, ifelse(pii <= 0.0001, 0.0001, pii))
  pil <- sapply(eta, pii.lambda, lambda, link)

  vi <- pii*(1-pii)
  fgi <- sapply(eta,pdfsl,lambda, link)
  gammaii<- Ni*fgi*pil/vi
  s <- t(gammaii)%*%X
  return(s)
}


K0b.group<-function(N,X,Ni,betas, link){
  eta <- X%*%betas
  if(link == 'logit'){
    gi <- dlogis(eta, log=FALSE)
    Gi <- plogis(eta,log=FALSE)
  }else if(link == 'probit'){
    gi <- dnorm(eta, log=FALSE)
    Gi <- pnorm(eta,log=FALSE)
  }else {
    gi <- dt(eta, df= 5,log=FALSE)
    Gi <- pt(eta, df = 5, log=FALSE)
  }
  pi0<- sapply(eta,pii.0, link)
  gammaii0<- Ni*gi*pi0/(Gi*(1-Gi))

  s <- t(gammaii0)%*%X
  return(s)
}


Kbb.group<-function(N,X,Ni,betas,lambda,link){
  eta <- X%*%betas
  pii<-sapply(eta,cdfsl,lambda,link)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.00001,0.00001,pii))
  vi<-pii*(1-pii)
  fgi <- sapply(eta,pdfsl,lambda,link)
  wii<- Ni*fgi^2/vi
  Wa <- diag(as.vector(wii))
  s <- t(X)%*%Wa%*%X
  return(s)
}

Kbb0.group<-function(N,X,Ni,betas,link){
  eta <- X%*%betas
  if(link == 'logit'){
    gi <- dlogis(eta, log=FALSE)
    Gi <- plogis(eta,log=FALSE)
  }else if(link == 'probit'){
    gi <- dnorm(eta, log=FALSE)
    Gi <- pnorm(eta,log=FALSE)
  }else {
    gi <- dt(eta, df= 5,log=FALSE)
    Gi <- pt(eta, df = 5, log=FALSE)
  }
  vi0 <- Gi*(1-Gi)
  wii0<- Ni*gi^2/vi0

  W0 <- diag(as.vector(wii0))
  s <- t(X)%*%W0%*%X
  return(s)
}


EIM.group <- function(par, y, X, Ni, link = 'logit'){
  N <- nrow(X)
  lambda <- par[1]
  betas  <- par[2:length(par)]
  p <- length(betas)
  eta    <- X%*%betas
  K11 <- Kll.group(N,X,Ni,betas,lambda,link)
  K12 <- Klb.group(N,X,Ni,betas,lambda,link)
  K22 <- Kbb.group(N,X,Ni,betas,lambda,link)
  K21 <- t(K12)
  K1 <- as.matrix(cbind(K11,K12))
  K2 <- as.matrix(cbind(K21,K22))

  EIM<- -rbind(K1,K2)
  return(EIM)
}


EIM.group0 <- function(par, y, X, Ni, link = 'logit'){
  N <- nrow(X)
  lambda <- par[1]
  betas  <- par[2:length(par)]
  p <- length(betas)
  eta    <- X%*%betas
  if(lambda != 0){
    K11 <- Kll.group(N,X,Ni,betas,lambda,link)
    K12 <- Klb.group(N,X,Ni,betas,lambda,link)
    K22 <- Kbb.group(N,X,Ni,betas,lambda,link)
  }else{
    K11 <- K00.group(N,X,Ni,betas,link)
    K12 <- K0b.group(N,X,Ni,betas,link)
    K22 <- Kbb0.group(N,X,Ni,betas,link)
  }
  K21 <- t(K12)
  K1 <- as.matrix(cbind(K11,K12))
  K2 <- as.matrix(cbind(K21,K22))

  EIM <- -rbind(K1,K2)
  return(EIM)
}

