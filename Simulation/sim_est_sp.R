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
beta1<-1 #initial value beta0 defined by fit skew-probit model
eta<-beta1*x1

#Function to calculate the CDF skew-probit
pdf.sp<-function(z,lambda){2*dnorm(z)*pnorm(lambda*z)}
cdf.sp<-function(z,lambda){integrate(pdf.sp,lower=-Inf,upper=z,lambda =lambda,stop.on.error=FALSE)$value}
lambda.ini<-0.5 #lambda value
print(lambda.ini)
prob<-sapply(eta,cdf.sp,lambda.ini) 
prob<-ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))

###############################################################################
#Parameters estimation
#Negative log-likelihood function
negll.sp<-function(par,y){
  lambda <- par[1]
  beta1  <- par[2]
  yres   <- y
  eta    <- beta1*x1 
  #pdf.sp <- function(z,lambda){2*dnorm(z,log=FALSE)*pnorm(lambda*z)}
  #cdf.sp <- function(z,lambda){integrate(pdf.sp,lower=-Inf, upper=z,lambda
  #                                       =lambda,stop.on.error=FALSE)$value}
  pii<-sapply(eta,cdf.sp,lambda)
  #Set high values for 0 < pi < 1
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  val <- -sum(yres * log(pii) + (1 - yres) * log(1 - pii))
  return(val)
}

#Gradient of negative log-likelihood function
grad.negll.sp<-function(par,y){
  lambda<-par[1]
  beta1<-par[2]
  yres<-y
  eta<-beta1*x1 
  #pdf.sp <- function(z,lambda){2*dnorm(z,log=FALSE)*pnorm(lambda*z)}
  #cdf.sp <- function(z,lambda){integrate(pdf.sp,lower=-Inf, upper=z,lambda
  #                                       =lambda,stop.on.error=FALSE)$value}
  pii<-sapply(eta,cdf.sp,lambda)
  #Set high values for 0 < pi < 1
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  fl <- function(z,lambda){2*z*dnorm(z)*dnorm(lambda*z)}
  Fl <- function(z,lambda){integrate(fl,lower=-Inf,upper=z,lambda=lambda,stop.on.error=FALSE)$value}
  fg <- function(z,lambda){2*dnorm(z)*pnorm(lambda*z)}
  pil<-sapply(eta,Fl,lambda)
  fgl<-sapply(eta,fg,lambda)
  psl <- length(par)
  val <- numeric(psl)
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
  #fi <- function(z,lambda){2*dnorm(z,log=FALSE)*pnorm(lambda*z)}
  #Fi <- function(z,lambda){integrate(fi,lower=-Inf, upper=z,lambda
  #                                   =lambda,stop.on.error=FALSE)$value}
  pii<-sapply(eta,cdf.sp,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  fl <- function(z,lambda){2*z*dnorm(z,log=FALSE)*dnorm(lambda*z)}
  Fl <- function(z,lambda){integrate(fl,lower=-Inf, upper=z,lambda=lambda,stop.on.error=FALSE)$value}
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
  pii<-sapply(eta,cdf.sp,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  fg <- sapply(1:length(y),function(j){pdf.sp(eta[j],lambda)})
  wi<- fg/(pii*(1-pii))
  W <- diag(wi)
  mu <- as.vector(size*pii)
  Ub <- t(x1)%*%W%*%(yres-mu)
  return(Ub)
}


#Generate replicates of response
R<-500 #number of replicates
R.aux<-1000 #number of additional replicates
r<-1 #index of replicates
Y<-matrix(NA,N,R.aux) #matrix of replicates for the response
pp<-1 #dimension of regression coefficients vector of logistic model
EPP<-matrix(NA,pp,R.aux) #matrix of parameter estimates of logistic model
psp<-2 #dimension of regression coefficients vector of skew-logistic model
EPSP<-matrix(NA,psp,R)
SENESPPCC.p<-matrix(NA,3,R)
SENESPPCC.sp<-matrix(NA,3,R)

count<-0
count.error<-0
repeat{
  if (count == R.aux) break
  count <-count+1
  y<-rbinom(N,1,prob) #simulating Bernoulli random variables of skew-probit model
  options(warn = 2) #to transform the warnings in error
  # try() to avoid the function failing
  df<-data.frame(y,x1)
  fit.p<-try(glm(y~x1-1,family=binomial(link="probit"),data=df),TRUE)
  if(inherits(fit.p,"try-error")) 
  {count.error<-count.error+1 
  next
  }
  else 
    Y[,count]<-y #to store the responses
  EPP[,count]<-fit.p$coefficients #to store the parameter estimates
}

Y<-Y[,!colSums(is.na(Y)) == nrow(Y),drop=FALSE]
EPP<-EPP[,!colSums(is.na(EPP)) == nrow(EPP),drop=FALSE]
print(dim(Y))

Y<-Y[,1:R]
EPP<-EPP[,1:R]  

countY<-sum(Y) #calculating the proportion of ones
print(countY)
print(N*R)

for(r in 1:R){
  theta.ini<-c(0,EPP[r])
  fit.sp<-optimx(par=theta.ini,y=Y[,r],fn=negll.sp,method ="BFGS")
  coef.opt<-coef(fit.sp)
  EPSP[1,r]<-coef.opt[,1]#lambda
  EPSP[2,r]<-coef.opt[,2]#beta1
  #sens spec acc
  eta.p<- EPP[r]*x1
  lambda0<-0
  prob.p<-sapply(eta.p,cdf.sp,lambda0)
  prob.p<-ifelse(prob.p>=0.9999,0.9999, ifelse(prob.p<=0.0001,0.0001,prob.p))
  id<-c(1:N)
  df.p<-data.frame(id,Y[,r],prob.p)
  cuttof<-optimal.thresholds(DATA = df.p, threshold = 101, which.model=1,
                             model.names = NULL, na.rm = TRUE, opt.methods=2, 0.85, 0.85,
                             obs.prev = NULL, smoothing = 1, FPC, FNC)
  SENESPPCC.p[1,r]<-sensitivity(cmx(df.p,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.p[2,r]<-specificity(cmx(df.p,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.p[3,r]<-pcc(cmx(df.p,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  eta.sp<- EPSP[2,r]*x1
  prob.sp<-sapply(eta.sp,cdf.sp,EPSP[1,r])
  prob.sp<-ifelse(prob.sp>=0.9999,0.9999, ifelse(prob.sp<=0.0001,0.0001,prob.sp))
  df.sp<-data.frame(id,Y[,r],prob.sp)
  cuttof<-optimal.thresholds(DATA = df.sp, threshold = 101, which.model=1,
                             model.names = NULL, na.rm = TRUE, opt.methods=2, 0.85, 0.85,
                             obs.prev = NULL, smoothing = 1, FPC, FNC)
  SENESPPCC.sp[1,r]<-sensitivity(cmx(df.sp,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.sp[2,r]<-specificity(cmx(df.sp,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.sp[3,r]<-pcc(cmx(df.sp,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
}


###############################################################################
#Estimation performance measures
#Mean
meanRbeta1est <- mean(EPP,na.rm = TRUE)

#Bias
biasRbeta1est <- meanRbeta1est-beta1

#Root mean squared error RMSE
RMSERbeta1 <- sqrt(1/R*sum((EPP-beta1)^2,na.rm = TRUE))

#Percent error PE
PERbeta1est <- ((meanRbeta1est-beta1)/beta1)*100

L.p<-list(meanb1=meanRbeta1est,
          biasb1=biasRbeta1est,
          RMSERb1=RMSERbeta1,
          PERb1=PERbeta1est)
L.p

###
#Mean
meanRlambdaestSP<-mean(EPSP[1,],na.rm = TRUE)
meanRbeta1estSP<-mean(EPSP[2,],na.rm = TRUE)

#Bias
biasRlambdaestSP<-meanRlambdaestSP-lambda.ini
biasRbeta1estSP<-meanRbeta1estSP-beta1

#Root mean squared error RMSE
RMSERlambdaSP<-sqrt(1/R*sum((EPSP[1,]-lambda.ini)^2,na.rm = TRUE))
RMSERbeta1SP<-sqrt(1/R*sum((EPSP[2,]-beta1)^2,na.rm = TRUE))

#Percent error PE
PERlambdaestSP<- ((meanRlambdaestSP-lambda.ini)/lambda.ini)*100
PERbeta1estSP<- ((meanRbeta1estSP-beta1)/beta1)*100

L.sp<-list(meanl=meanRlambdaestSP,
           meanb1=meanRbeta1estSP,
           biasl=biasRlambdaestSP,
           biasb1=biasRbeta1estSP,
           RMSERl=RMSERlambdaSP,
           RMSERb1=RMSERbeta1SP,
           PERl=PERlambdaestSP,
           PERb1=PERbeta1estSP)

L.sp


##Prediction performance measures
#Mean
SEN.p<-mean(SENESPPCC.p[1,],na.rm = TRUE)
ESP.p<-mean(SENESPPCC.p[2,],na.rm = TRUE)
PCC.p<-mean(SENESPPCC.p[3,],na.rm = TRUE)

L.p.2<-list(sens=SEN.p,spec=ESP.p,pcc=PCC.p)
L.p.2

SEN.sp<-mean(SENESPPCC.sp[1,],na.rm = TRUE)
ESP.sp<-mean(SENESPPCC.sp[2,],na.rm = TRUE)
PCC.sp<-mean(SENESPPCC.sp[3,],na.rm = TRUE)

L.sp.2<-list(sens=SEN.sp,spec=ESP.sp,pcc=PCC.sp)
L.sp.2

proc.time() - start_time


