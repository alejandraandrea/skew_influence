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
beta1<-1 #initial value beta0 defined by fit skew-logistic model
eta<-beta1*x1

#Function to calculate the CDF skew-logistic
pdf.sl<-function(z,lambda){2*dlogis(z,log=FALSE)*pnorm(lambda*z)}
cdf.sl<-function(z,lambda){integrate(pdf.sl,lower=-Inf, upper=z,
                                     lambda =lambda,stop.on.error=FALSE)$value}
lambda.ini<-0.5 #lambda value
print(lambda.ini)
prob<-sapply(eta,cdf.sl,lambda.ini) #all.equal F<-exp(eta)/(1+exp(eta)) when lambda=0
prob<-ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))

###############################################################################
#Parameters estimation
#Negative log-likelihood function
negll.sl<-function(par,y){
  lambda <- par[1]
  beta1  <- par[2]
  yres   <- y
  eta    <- beta1*x1 
  #pdf.sl <- function(z,lambda){2*dlogis(z,log=FALSE)*pnorm(lambda*z)}
  #cdf.sl <- function(z,lambda){integrate(pdf.sl,lower=-Inf, upper=z,lambda
  #                                       =lambda,stop.on.error=FALSE)$value}
  pii<-sapply(eta,cdf.sl,lambda)
  #Set high values for 0 < pi < 1
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  val <- -sum(yres * log(pii) + (1 - yres) * log(1 - pii))
  return(val)
}

#Gradient of negative log-likelihood function
grad.negll.sl<-function(par,y){
  lambda<-par[1]
  beta1<-par[2]
  yres<-y
  eta<-beta1*x1 
  #pdf.sl <- function(z,lambda){2*dlogis(z,log=FALSE)*pnorm(lambda*z)}
  #cdf.sl <- function(z,lambda){integrate(pdf.sl,lower=-Inf, upper=z,lambda
  #                                       =lambda,stop.on.error=FALSE)$value}
  pii<-sapply(eta,cdf.sl,lambda)
  #Set high values for 0 < pi < 1
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  fl <- function(z,lambda){2*z*dlogis(z,log=FALSE)*dnorm(lambda*z)}
  Fl <- function(z,lambda){integrate(fl,lower=-Inf, upper=z,lambda=lambda,
                                     stop.on.error=FALSE)$value}
  fg <- function(z,lambda){2*dlogis(z,log=FALSE)*pnorm(lambda*z)}
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
  #fi <- function(z,lambda){2*dlogis(z,log=FALSE)*pnorm(lambda*z)}
  #Fi <- function(z,lambda){integrate(fi,lower=-Inf, upper=z,lambda
  #                                   =lambda,stop.on.error=FALSE)$value}
  pii<-sapply(eta,cdf.sl,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  fl <- function(z,lambda){2*z*dlogis(z,log=FALSE)*dnorm(lambda*z)}
  Fl <- function(z,lambda){integrate(fl,lower=-Inf, upper=z,lambda=lambda,
                                     stop.on.error=FALSE)$value}
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
  pii<-sapply(eta,cdf.sl,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  fg <- sapply(1:length(y),function(j){pdf.sl(eta[j],lambda)})
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
pl<-1 #dimension of regression coefficients vector of logistic model
EPL<-matrix(NA,pl,R.aux) #matrix of parameter estimates of logistic model
psl<-2 #dimension of regression coefficients vector of skew-logistic model
EPSL<-matrix(NA,psl,R)
SENESPPCC.l<-matrix(NA,3,R)
SENESPPCC.sl<-matrix(NA,3,R)


count<-0
count.error<-0
repeat{
  if (count == R.aux) break
  count <-count+1
  y<-rbinom(N,1,prob) #simulating Bernoulli random variables of skew-logistic model
  options(warn = 2) #to transform the warnings in error
  # try() to avoid the function failing
  df<-data.frame(y,x1)
  fit.l<-try(glm(y~x1-1,family=binomial(link="logit"),data=df),TRUE)
  if(inherits(fit.l,"try-error")) 
  {count.error<-count.error+1 
  next
  }
  else 
    Y[,count]<-y #to store the responses
  EPL[,count]<-fit.l$coefficients #to store the parameter estimates
}

Y<-Y[,!colSums(is.na(Y)) == nrow(Y),drop=FALSE]
EPL<-EPL[,!colSums(is.na(EPL)) == nrow(EPL),drop=FALSE]
print(dim(Y))

Y<-Y[,1:R]
EPL<-EPL[,1:R]  

countY<-sum(Y) #calculating the proportion of ones
print(countY)
print(N*R)

for(r in 1:R){
  theta.ini<-c(0,EPL[r])
  fit.sl<-optimx(par=theta.ini,y=Y[,r],fn=negll.sl,
                 method ="BFGS")
  coef.opt<-coef(fit.sl)
  EPSL[1,r]<-coef.opt[,1]#lambda
  EPSL[2,r]<-coef.opt[,2]#beta1
  #sens spec acc
  eta.l<- EPL[r]*x1
  lambda0<-0
  prob.l<-sapply(eta.l,cdf.sl,lambda0)
  prob.l<-ifelse(prob.l>=0.9999,0.9999, ifelse(prob.l<=0.0001,0.0001,prob.l))
  id<-c(1:N)
  df.l<-data.frame(id,Y[,r],prob.l)
  cuttof<-optimal.thresholds(DATA = df.l, threshold = 101, which.model=1,
                             model.names = NULL, na.rm = TRUE, opt.methods=2, 0.85, 0.85,
                             obs.prev = NULL, smoothing = 1, FPC, FNC)
  SENESPPCC.l[1,r]<-sensitivity(cmx(df.l,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.l[2,r]<-specificity(cmx(df.l,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.l[3,r]<-pcc(cmx(df.l,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  eta.sl<- EPSL[2,r]*x1
  prob.sl<-sapply(eta.sl,cdf.sl,EPSL[1,r])
  prob.sl<-ifelse(prob.sl>=0.9999,0.9999, ifelse(prob.sl<=0.0001,0.0001,prob.sl))
  df.sl<-data.frame(id,Y[,r],prob.sl)
  cuttof<-optimal.thresholds(DATA = df.sl, threshold = 101, which.model=1,
                             model.names = NULL, na.rm = TRUE, opt.methods=2, 0.85, 0.85,
                             obs.prev = NULL, smoothing = 1, FPC, FNC)
  SENESPPCC.sl[1,r]<-sensitivity(cmx(df.sl,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.sl[2,r]<-specificity(cmx(df.sl,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
  SENESPPCC.sl[3,r]<-pcc(cmx(df.sl,threshold = cuttof[1,2],which.model=1),st.dev=FALSE)
}


###############################################################################
#Estimation performance measures
#Mean
meanRbeta1est <- mean(EPL,na.rm = TRUE)

#Bias
biasRbeta1est <- meanRbeta1est-beta1

#Root mean squared error RMSE
RMSERbeta1 <- sqrt(1/R*sum((EPL-beta1)^2,na.rm = TRUE))

#Percent error PE
PERbeta1est <- ((meanRbeta1est-beta1)/beta1)*100

L.l<-list(meanb1=meanRbeta1est,
          biasb1=biasRbeta1est,
          RMSERb1=RMSERbeta1,
          PERb1=PERbeta1est)
L.l

###
#Mean
meanRlambdaestSL<-mean(EPSL[1,],na.rm = TRUE)
meanRbeta1estSL<-mean(EPSL[2,],na.rm = TRUE)

#Bias
biasRlambdaestSL<-meanRlambdaestSL-lambda.ini
biasRbeta1estSL<-meanRbeta1estSL-beta1

#Root mean squared error RMSE
RMSERlambdaSL<-sqrt(1/R*sum((EPSL[1,]-lambda.ini)^2,na.rm = TRUE))
RMSERbeta1SL<-sqrt(1/R*sum((EPSL[2,]-beta1)^2,na.rm = TRUE))

#Percent error PE
PERlambdaestSL<- ((meanRlambdaestSL-lambda.ini)/lambda.ini)*100
PERbeta1estSL<- ((meanRbeta1estSL-beta1)/beta1)*100

L.sl<-list(meanl=meanRlambdaestSL,
           meanb1=meanRbeta1estSL,
           biasl=biasRlambdaestSL,
           biasb1=biasRbeta1estSL,
           RMSERl=RMSERlambdaSL,
           RMSERb1=RMSERbeta1SL,
           PERl=PERlambdaestSL,
           PERb1=PERbeta1estSL)

L.sl


##Prediction performance measures
#Mean
SEN.l<-mean(SENESPPCC.l[1,],na.rm = TRUE)
ESP.l<-mean(SENESPPCC.l[2,],na.rm = TRUE)
PCC.l<-mean(SENESPPCC.l[3,],na.rm = TRUE)

L.l.2<-list(sens=SEN.l,spec=ESP.l,pcc=PCC.l)
L.l.2

SEN.sl<-mean(SENESPPCC.sl[1,],na.rm = TRUE)
ESP.sl<-mean(SENESPPCC.sl[2,],na.rm = TRUE)
PCC.sl<-mean(SENESPPCC.sl[3,],na.rm = TRUE)

L.sl.2<-list(sens=SEN.sl,spec=ESP.sl,pcc=PCC.sl)
L.sl.2

proc.time() - start_time


