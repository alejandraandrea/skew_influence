###############################################################################
#Monte Carlo Simulation Studies: Simulation Study II
###############################################################################

#Time
start_time <- proc.time()

#Packages
library(optimx)
library(matrixcalc)

#Initial values
N<-250 #numbers of observations
print(N)
addTaskCallback(function(...) {set.seed(3711);TRUE})
x1<-runif(N,0,1)
x1<-sort(x1)
beta1<-1 #initial value beta1 defined by fit skew-probit model
eta<-beta1*x1

#t model
nu<-5 #Suggested by Gómez et al. (2007) 

#Function to calculate the CDF of skew-t
pdf.st<-function(z,lambda){2*dt(z,nu)*pnorm(lambda*z)}
cdf.st<-function(z,lambda){integrate(pdf.st,lower=-Inf, upper=z,lambda =lambda,stop.on.error=FALSE)$value}

lambda<-0
#print(lambda)
prob<-sapply(eta,cdf.st,lambda) 
prob<-ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))

#t model
#Function to calculate the CDF t-student
nu<-5 #Suggested by Gómez et al. (2007) 
pdf.t<-function(z){dt(z,nu)}
cdf.t<-function(z){pt(z,nu)}

#Negative log-likelihood function
negll.t <- function(par,data){
  formula <- cbind(data$y,data$n-data$y)~data$x #Adaptar para el dataset
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  beta <- par
  eta <- X[,-1]*beta
  prob <- sapply(1:length(y),function(j){cdf.t(eta[j])})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val<- - sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  return(val)
}


#Generate replicates of response
addTaskCallback(function(...) {set.seed(3711);TRUE})
R<-500 #number of replicates
r<-1 #index of replicates
Y<-matrix(NA,N,R) #matrix of replicates for the response
Yper<-matrix(NA,N,R)
pt<-1 #dimension of regression coefficients vector of probit model
EPT<-matrix(NA,pt,R) #matrix of parameter estimates of probit model
EPper<-matrix(NA,pt,R)

count1<-0 #count of total replicates
count2<-0 #count of replicates without error in the fit

while(r<=R){
  count1<-count1+1
  y<-rbinom(N,1,prob) #simulating Bernoulli random variables of skew-probit model
  options(warn = 2) #to transform the warnings in error
  # try() to avoid the function failing
  df1<-data.frame(y,x1)
  #estimation
  fit.l<-try(glm(y~x1-1,family=binomial(link="logit"),data=df1),TRUE) #TRUE  the report of warnings is suppressed
  if(!inherits(fit.l,"try-error")){#If there is no error
    theta.ini<-fit.l$coefficients
    s<-1
    n<-rep(s,N)
    x<-x1
    df2<-data.frame(y,n,x)
    fit.t<-optimx(par=theta.ini,data=df2,fn=negll.t,method="BFGS")
    coef.opt<-coef(fit.t)
    EPT[,r]<-coef.opt[1,1]#BFGS #to store the parameter estimates
    Y[,r]<-y #to store the responses
    delta<-0.5 #delta value
    k1<-8
    x1delta<-x1
    x1delta[k1]<-x1[k1]+delta
    y<-Y[,r]
    s<-1
    n<-rep(s,N)
    x<-x1delta
    dfdelta<-data.frame(y,n,x)
    options(warn = 2) 
    fitper<-try(optimx(par=theta.ini,data=dfdelta,fn=negll.t,method="BFGS"),TRUE)
    if(!inherits(fitper,"try-error")){#If there is no error
      coef.opt<-coef(fitper)
      EPper[,r]<-coef.opt[1,1]  #to store the parameter estimates
    }
  }
  r<-r+1
  count2<-count2+1
}

count1
count2


#Packages
library(matrixcalc)

#t
nu<-5 #Suggested by Gómez et al. (2007) 

#Observed Information Matrix
g<-function(z){
  g<-dt(z,nu); #PDF t
  return(g)
}

dg<-function(z){
  dg<--gamma((nu+1)/2)/(sqrt(pi*nu)*gamma(nu/2))*
    (1+1/nu)*z*(1+z^2/nu)^(-(nu+3)/2) #derivative PDF t 
  return(dg)
}

pdf.st<-function(z,lambda){2*g(z)*pnorm(lambda*z)}
cdf.st<-function(z,lambda){integrate(pdf.st,lower=-Inf, upper=z,lambda =lambda,stop.on.error=FALSE)$value}

F2<-function(lambda){
  f2<-function(z){
    f2<-2*z*g(z)*dnorm(lambda*z,0,1)
  }
  return(f2)
}

F3<-function(lambda,eta){
  F3<-integrate(F2(lambda),lower=-Inf,upper=eta)
  return(F3$value)
}

F4<-function(lambda){
  f4<-function(z){
    f4<- -2*lambda*z^3*g(z)*dnorm(lambda*z,0,1)
  }
  return(f4)
}

F5<-function(lambda,eta){
  F5<-integrate(F4(lambda),lower=-Inf,upper=eta)
  return(F5$value)
}

Llli<-function(Ni,yi,etai,lambda){
  pii<-cdf.st(etai,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  piil<-F3(lambda,etai)
  si<-piil/(pii*(1-pii))
  vi<-pii*(1-pii)
  piill<-F5(lambda,etai)
  sil<-vi^(-1)*(piill-(1-2*pii)*piil^2*vi^(-1))
  Llli<-(yi-Ni*pii)*sil-Ni*si*piil
  return(Llli)
}

Lll<-function(N,data,betas,lambda){
  s<-matrix(0,1,1)
  for(i in 1:N){
    yi<-data[i,1]
    Ni<-data[i,2]
    xi<-data[i,-(1:2)]
    etai<-t(betas)%*%xi
    a<-Llli(Ni,yi,etai,lambda)
    s<-s+a
  }
  return(s)
}

F6<-function(lambda,eta){
  as.numeric(2*g(eta)*pnorm(lambda*eta,0,1))
}

F7<-function(lambda,eta){
  as.numeric(2*eta*g(eta)*dnorm(lambda*eta,0,1))
}

Llbi<-function(Ni,yi,xi,etai,lambda){
  pii<-cdf.st(etai,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  vi<-pii*(1-pii)
  piil<-F3(lambda,etai)
  si<-piil/(pii*(1-pii))
  piib<-F6(lambda,etai)*xi
  piilb<-F7(lambda,etai)*xi
  sib<-vi^(-1)*(piilb-(1-2*pii)*piil*piib*vi^(-1))
  Llbi<-(yi-Ni*pii)*sib-Ni*si*piib
  return(Llbi)
}

Llb<-function(p,N,data,betas,lambda){
  s<-matrix(0,1,p)
  for(i in 1:N){
    yi<-data[i,1]
    Ni<-data[i,2]
    xi<-data[i,-(1:2)]
    etai<-t(betas)%*%xi
    a<-Llbi(Ni,yi,xi,etai,lambda)
    s<-s+a
  }
  return(s)
}

F8<-function(lambda,eta){
  as.numeric(2*(g(eta)*dnorm(lambda*eta,0,1)*lambda+pnorm(lambda*eta,0,1)*dg(eta)))
}

Lbbi<-function(Ni,yi,xi,etai,lambda){
  pii<-cdf.st(etai,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  vi<-pii*(1-pii)
  wii<-F6(lambda,etai)/(pii*(1-pii))
  vii<-vi^(-1)*(yi-Ni*pii)*(F8(lambda,etai)-F6(lambda,etai)^2*vi^(-1)*(1-2*pii))-F6(lambda,etai)*Ni*wii
  Lbbi<-vii*(xi%*%t(xi))
  return(Lbbi)
}

Lbb<-function(p,N,data,betas,lambda){
  s<-matrix(0,p,p)
  for(i in 1:N){
    yi<-data[i,1]
    Ni<-data[i,2]
    xi<-as.matrix(data[i,-(1:2)])
    etai<-t(betas)%*%xi
    a<-Lbbi(Ni,yi,xi,etai,lambda)
    s<-s+a
  }
  return(s)
}

OIM<-function(p,N,data,betas,lambda){
  if(lambda!=0){
    L11 <- -Lll(N,data,betas,lambda)
    L12 <- -Llb(p,N,data,betas,lambda)
    L21 <- t(L12)
    L22 <- -Lbb(p,N,data,betas,lambda)
    L1 <- as.matrix(cbind(L11,L12))
    L2 <- as.matrix(cbind(L21,L22))
    OIM<-rbind(L1,L2)
  }
  else{
    L22 <- -Lbb(p,N,data,betas,lambda)
    OIM<-L22
  }
  return(OIM)
}


#Perturbation matrix under skew perturbation scheme
nu<-5 

g<-function(z){
  g<-dt(z,nu); #PDF t
  return(g)
}

cdf.t<-function(z){ #CDF standard normal
  cdf.t<-integrate(g,lower=-Inf,upper=z)
  return(cdf.t$value)
}

g2<-function(z){
  g2<-(2/sqrt(2*pi))*z*g(z)
  return(g2)
}

g3<-function(z){
  g3<-z*g(z)
  return(g3)
}

F9<-function(z){
  F9<-integrate(g2,lower=-Inf,upper=z)
  return(F9$value)
}

F10<-function(z){
  F10<-integrate(g3,lower=-Inf,upper=z)
  return(F10$value)
}

PMi<-function(i,data,betas){
  yi <- data[i,1]
  Ni <- data[i,2]
  xi <- as.matrix(data[i,-(1:2)])
  etai <- as.numeric(t(betas)%*%xi)
  pii <- cdf.t(etai)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  vi0 <- pii*(1-pii)
  gii <- (Ni*F9(etai))^2/(Ni*pii*(1-pii))
  ki <- 1/sqrt(gii)
  kib <- -(1/2)*ki*gii^(-1)*Ni*F9(etai)*g(etai)*vi0^(-1)*(4*etai/sqrt(2*pi)-vi0^(-1)*F9(etai)*(1-2*pii))*xi
  pii0bl <- (1/sqrt(2*pi))*(2*etai*g(etai)*ki*xi + 2*kib*F10(etai))
  pii0l <- ki*(2/sqrt(2*pi))*F10(etai)
  ri0b <- -vi0^(-1)*(vi0^(-1)*(yi-Ni*pii)*(1-2*pii)+Ni)*g(etai)*xi
  PMi<- vi0^(-1)*(yi-Ni*pii)*pii0bl+pii0l*ri0b
  return(PMi)
}

PM<-function(p,N,data,betas){
  PM<-PMi(1,data,betas)
  for(i in 2:N){
    PM<-cbind(PM,PMi(i,data,betas))
  }
  PM<-matrix(PM,nrow=p,ncol=N,dimnames = NULL)
  return(PM)
}

#Conformal normal curvature 
CNC<-function(M1,M2){
  C1<-2*(t(M2)%*%solve(M1)%*%M2)
  C2<-C1/matrix.trace(C1)
  if(isSymmetric(C2)){
    CNC<-diag(C2)
  }
  else{
    C2<-1/2*(C2+t(C2))
    CNC<-diag(C2)
  }
  return(CNC)
}

#With d different from 0
ipcncd<-function(k1,q,cuttof,L){
  cnc<-L
  c<-numeric(q)
  if(cnc[k1]>cuttof){
    c[k1]<-1
  }
  return(c)
}

q<-N
IPSE<-matrix(0,q,R)
for(r in 1:R){
  #r=3
  size<-rep(1,N)
  #data<-matrix(cbind(Y[,r],size,x1),nrow=N,ncol=3,dimnames=NULL)
  data<-matrix(cbind(Y[,r],size,x1delta),nrow=N,ncol=3,dimnames=NULL)
  N<-nrow(data)
  p<-ncol(data)-2
  lambdaH0<-0
  lambda<-lambdaH0
  #betas<-matrix(EPT[,r],nrow=p,ncol=1,dimnames = NULL)
  betas<-matrix(EPper[,r],nrow=p,ncol=1,dimnames = NULL)
  OIMT<-OIM(p,N,data,betas,lambda)
  PMT<-PM(p,N,data,betas)
  CNCT<-CNC(OIMT,PMT)
  cuttofSE<-mean(CNCT)+2*sd(CNCT)/sqrt(q)
  k1<-8
  IPSE[,r]<-ipcncd(k1,q,cuttofSE,CNCT)
}

countY<-sum(Y)
print(countY/(N*R))

countip<-rowSums(IPSE)[k1]
print(countip)
perip<-countip/R*100 #calculating percentages of correct detection 
print(perip)

proc.time() - start_time

