###############################################################################
#Monte Carlo Simulation Studies: Simulation Study I
###############################################################################

#Time
start_time <- proc.time()

#rm(list=ls())

#Packages
library("optimx")
library("PresenceAbsence")
library("matrixcalc")

#Initial values
N<-250 #numbers of observations
print(N)
set.seed(3711)
x1<-runif(N,0,1)
beta1<-1 #initial value beta0 defined by fit skew-t model
eta<-beta1*x1

#Function to calculate the CDF skew-t
nu<-5 #Suggested by Gómez et al. (2007) 
pdf.st<-function(z,lambda){2*dt(z,nu)*pnorm(lambda*z)}
cdf.st<-function(z,lambda){integrate(pdf.st,lower=-Inf,upper=z,lambda=lambda,stop.on.error=FALSE)$value}

lambda.ini<-0.5 #lambda value
print(lambda.ini)
prob<-sapply(eta,cdf.st,lambda.ini) 
prob<-ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))

###############################################################################
#Parameters estimation
#Negative log-likelihood function
negll.st<-function(par,y){
  lambda <- par[1]
  beta1  <- par[2]
  yres   <- y
  eta    <- beta1*x1 
  pii<-sapply(eta,cdf.st,lambda)
  #Set high values for 0 < pi < 1
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  val <- -sum(yres * log(pii) + (1 - yres) * log(1 - pii))
  return(val)
}

#Gradient of negative log-likelihood function
grad.negll.st<-function(par,y){
  lambda<-par[1]
  beta1<-par[2]
  yres<-y
  eta<-beta1*x1 
  pii<-sapply(eta,cdf.st,lambda)
  #Set high values for 0 < pi < 1
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  fl <- function(z,lambda){2*z*dt(z,nu)*dnorm(lambda*z)}
  Fl <- function(z,lambda){integrate(fl,lower=-Inf,upper=z,lambda=lambda,stop.on.error=FALSE)$value}
  fg <- function(z,lambda){2*dt(z,nu)*pnorm(lambda*z)}
  pil<-sapply(eta,Fl,lambda)
  fgl<-sapply(eta,fg,lambda)
  pst <- length(par)
  val <- numeric(pst)
  val[1] <- -sum((yres-pii)/(pii*(1-pii))*pil)
  val[2] <- -sum((yres-pii)/(pii*(1-pii))*fgl*x1)            
  return(val)
}

#Score lambda vector
score.lambda <- function(par,y){
  lambda <- par[1]
  beta1 <- par[2]
  yres <- y
  size <- 1
  eta <- beta1*x1 
  pii<-sapply(eta,cdf.st,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  fl <- function(z,lambda){2*z*dt(z,nu)*dnorm(lambda*z)}
  Fl <- function(z,lambda){integrate(fl,lower=-Inf,upper=z,lambda=lambda,stop.on.error=FALSE)$value}
  pil <- sapply(eta,Fl,lambda)
  mui <- as.vector(size*pii)
  si <- as.vector(pil/(pii*(1-pii)))
  Ul<- t(si)%*%(yres-mui)
  return(Ul)
}

score.betas <- function(par,y){
  lambda <- par[1]
  beta1 <- par[2]
  yres <- y
  size <- 1
  eta <- beta1*x1 
  pii<-sapply(eta,cdf.st,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  fg <- sapply(1:length(y),function(j){pdf.st(eta[j],lambda)})
  wi<- fg/(pii*(1-pii))
  W <- diag(wi)
  mu <- as.vector(size*pii)
  Ub <- t(x1)%*%W%*%(yres-mu)
  return(Ub)
}


#t model
#Function to calculate the CDF t-student
nu<-5 #Sugerido por Gómez et al. (2007) 
pdf.t<-function(z){dt(z,nu)}
cdf.t<-function(z){pt(z,nu)}

#Negative log-likelihood function
negll.t<-function(par,y){
  beta1  <- par
  yres   <- y
  eta    <- beta1*x1 
  pii<-sapply(eta,cdf.t)
  #Set high values for 0 < pi < 1
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  val <- -sum(yres * log(pii) + (1 - yres) * log(1 - pii))
  return(val)
}

#Generate replicates of response
R<-500 #number of replicates
R.aux<-550#1000 #number of additional replicates
r<-1 #index of replicates
Y<-matrix(NA,N,R.aux) #matrix of replicates for the response
YY<-matrix(0,N,R.aux)
pt<-1 #dimension of regression coefficients vector of t model
EPT<-matrix(NA,pt,R.aux) #matrix of parameter estimates of t model
pst<-2 #dimension of regression coefficients vector of skew-t model
EPST<-matrix(NA,pst,R)
SENESPPCC.t<-matrix(NA,3,R)
SENESPPCC.st<-matrix(NA,3,R)


count<-0
count.error<-0
repeat{
  if (count == R.aux) break
  count <-count+1
  y<-rbinom(N,1,prob) #simulating Bernoulli random variables of skew-t model
  YY[,count]<-y
  #options(warn = 2) #to transform the warnings in error
  # try() to avoid the function failing
  df<-data.frame(y,x1)
  fit.l<-try(glm(y~x1-1,family=binomial(link="logit"),data=df),TRUE)
  theta.ini<-fit.l$coefficients
  if(inherits(fit.l,"try-error")) {
    count.error<-count.error+1 
    next
  }
  else 
  fit.t<-try(optimx(par=theta.ini,y=y,fn=negll.t,method="BFGS"),TRUE)
  coef.opt<-coef(fit.t)
  EPT[,r]<-coef.opt[1,1]#BFGS #to store the parameter estimates
  if(inherits(fit.t,"try-error")) {
    count.error<-count.error+1 
  next
  }
  else 
    Y[,count]<-y #to store the responses
  EPT[,count]<-coef.opt[1,1]#BFGS #to store the parameter estimates
}

Y<-Y[,!colSums(is.na(Y)) == nrow(Y),drop=FALSE]
EPT<-EPT[,!colSums(is.na(EPT)) == nrow(EPT),drop=FALSE]
print(dim(Y))

Y<-Y[,1:R]
EPT<-EPT[,1:R]  

countY<-sum(Y) #calculating the proportion of ones
print(countY)
print(N*R)

for(r in 1:R){
  theta.ini<-c(0,EPT[r])
  fit.st<-optimx(par=theta.ini,y=Y[,r],fn=negll.st,method ="BFGS")
  coef.opt<-coef(fit.st)
  EPST[1,r]<-coef.opt[,1]#lambda
  EPST[2,r]<-coef.opt[,2]#beta1
  #sens spec acc
  eta.t<- EPT[r]*x1
  lambda0<-0
  prob.t<-sapply(eta.t,cdf.st,lambda0)
  prob.t<-ifelse(prob.t>=0.9999,0.9999, ifelse(prob.t<=0.0001,0.0001,prob.t))
  id<-c(1:N)
  df.t<-data.frame(id,Y[,r],prob.t)
  cuttof<-optimal.thresholds(DATA = df.t, threshold = 101, which.model=1,
                             model.names = NULL, na.rm = TRUE, opt.methods=2, 0.85, 0.85,
                             obs.prev = NULL, smoothing = 1, FPC, FNC)
  SENESPPCC.t[1,r]<-sensitivity(cmx(df.t,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.t[2,r]<-specificity(cmx(df.t,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.t[3,r]<-pcc(cmx(df.t,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  eta.st<- EPST[2,r]*x1
  prob.st<-sapply(eta.st,cdf.st,EPST[1,r])
  prob.st<-ifelse(prob.st>=0.9999,0.9999,ifelse(prob.st<=0.0001,0.0001,prob.st))
  df.st<-data.frame(id,Y[,r],prob.st)
  cuttof<-optimal.thresholds(DATA = df.st, threshold = 101, which.model=1,
                             model.names = NULL, na.rm = TRUE, opt.methods=2, 0.85, 0.85,
                             obs.prev = NULL, smoothing = 1, FPC, FNC)
  SENESPPCC.st[1,r]<-sensitivity(cmx(df.st,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.st[2,r]<-specificity(cmx(df.st,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.st[3,r]<-pcc(cmx(df.st,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
}


###############################################################################
#Estimation performance measures
#Mean
meanRbeta1est <- mean(EPT,na.rm = TRUE)

#Bias
biasRbeta1est <- meanRbeta1est-beta1

#Root mean squared error RMSE
RMSERbeta1 <- sqrt(1/R*sum((EPT-beta1)^2,na.rm = TRUE))

#Percent error PE
PERbeta1est <- ((meanRbeta1est-beta1)/beta1)*100

L.t<-list(meanb1=meanRbeta1est,
          biasb1=biasRbeta1est,
          RMSERb1=RMSERbeta1,
          PERb1=PERbeta1est)
L.t

###
#Mean
meanRlambdaestST<-mean(EPST[1,],na.rm = TRUE)
meanRbeta1estST<-mean(EPST[2,],na.rm = TRUE)

#Bias
biasRlambdaestST<-meanRlambdaestST-lambda.ini
biasRbeta1estST<-meanRbeta1estST-beta1

#Root mean squared error RMSE
RMSERlambdaST<-sqrt(1/R*sum((EPST[1,]-lambda.ini)^2,na.rm = TRUE))
RMSERbeta1ST<-sqrt(1/R*sum((EPST[2,]-beta1)^2,na.rm = TRUE))

#Percent error PE
PERlambdaestST<- ((meanRlambdaestST-lambda.ini)/lambda.ini)*100
PERbeta1estST<- ((meanRbeta1estST-beta1)/beta1)*100

L.st<-list(meanl=meanRlambdaestST,
           meanb1=meanRbeta1estST,
           biasl=biasRlambdaestST,
           biasb1=biasRbeta1estST,
           RMSERl=RMSERlambdaST,
           RMSERb1=RMSERbeta1ST,
           PERl=PERlambdaestST,
           PERb1=PERbeta1estST)

L.st


##Prediction performance measures
#Mean
SEN.t<-mean(SENESPPCC.t[1,],na.rm = TRUE)
ESP.t<-mean(SENESPPCC.t[2,],na.rm = TRUE)
PCC.t<-mean(SENESPPCC.t[3,],na.rm = TRUE)

L.t.2<-list(sens=SEN.t,spec=ESP.t,pcc=PCC.t)
L.t.2

SEN.st<-mean(SENESPPCC.st[1,],na.rm = TRUE)
ESP.st<-mean(SENESPPCC.st[2,],na.rm = TRUE)
PCC.st<-mean(SENESPPCC.st[3,],na.rm = TRUE)

L.st.2<-list(sens=SEN.st,spec=ESP.st,pcc=PCC.st)
L.st.2

proc.time() - start_time


