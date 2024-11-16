###############################################################################
#Local influence under skew appropriate perturbation: 
#Age of menarche of Warsaw girls
###############################################################################

###############################################################################
#Read data
df<-read.table("li_girls.txt", sep="\t", header=TRUE) 

#Get the data structure
str(df)

#Recognize each column of data
attach(df)

###########################################################################
#Packages
library(matrixcalc)

#Logistic model

#Observed Information Matrix
g<-function(z){
  #g<-exp(z)/(1+exp(z))^2; #PDF logistic
  g<-dlogis(z,log=FALSE); #PDF logistic
  return(g)
}

dg<-function(z){
  dg<-(exp(z)*(1+exp(z))-2*exp(2*z))/(1+exp(z))^3 #derivative PDF logistic
  return(dg)
}

pdf.sl<-function(z,lambda){2*g(z)*pnorm(lambda*z)}
cdf.sl<-function(z,lambda){integrate(pdf.sl,lower=-Inf, upper=z,
                                     lambda =lambda,stop.on.error=FALSE)$value}

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
  pii<-cdf.sl(etai,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.00001,0.00001,pii))
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
  pii<-cdf.sl(etai,lambda)
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
  pii<-cdf.sl(etai,lambda)
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
g<-function(z){
  #g<-exp(z)/(1+exp(z))^2; #PDF logistic
  g<-dlogis(z); #PDF logistic
  return(g)
}

cdf.l<-function(z){ #CDF logistic
  cdf.l<-integrate(g,lower=-Inf,upper=z)
  return(cdf.l$value)
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
  pii <- cdf.l(etai)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.0001,0.0001,pii))
  vi0 <- pii*(1-pii)
  gii <- (Ni*F9(etai))^2/(Ni*pii*(1-pii))
  ki <- 1/sqrt(gii)
  kib <- -(1/2)*ki*gii^(-1)*Ni*F9(etai)*g(etai)*vi0^(-1)*(4*etai/sqrt(2*pi)
                                                          -vi0^(-1)*F9(etai)*(1-2*pii))*xi
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


#Function to calculate the influential points with CNC
ipcnc<-function(q,cuttof,L){
  cnc<-L
  index<-vector(length=q)
  c<-0
  for(i in 1:q){
    if(cnc[i]>cuttof){ 
      index[i]<-i
      c<-c+1
    }
  }
  L<-list(cnc,cuttof,c,index)
  return(L)
}

#Transforming a data.frame to matrix
girls<-as.matrix(df);

#Do not assign names to matrix columns (dimnames NULL)
data<-matrix(girls, nrow = nrow(girls), ncol=ncol(girls), dimnames = NULL);

#Total number of observations
N<-nrow(data)

#Number of predictors
p<-ncol(data)-2

#Parameter estimates for Logistic model 
fit.l<- glm(cbind(Mens,Entr-Mens)~Ida,family=binomial("logit"),data = df)
lambda <- 0
beta0.est.l <- fit.l$coefficients[1]
beta1.est.l <- fit.l$coefficients[2]
betas<-matrix(c(beta0.est.l,beta1.est.l),nrow=p,ncol=1,dimnames = NULL)
OIML<-OIM(p,N,data,betas,lambda)
OIML
PML<-PM(p,N,data,betas)
PML
CNCL<-CNC(OIML,PML)
q<-N
cuttofSE<-mean(CNCL)+2*sd(CNCL)/sqrt(q)
plot(c(1:N),CNCL,ylim=c(0,1),xlab = "Index",ylab="Conformal normal curvature",
     col="black",pch=1,lwd=1)
plot(c(1:N),ylim=c(0,0.1),CNCL,xlab = "Index",ylab="Conformal normal curvature",
     col="black",pch=1,lwd=1)
abline(h=cuttofSE,lty=2,lwd=1)

IPL<-ipcnc(q,cuttofSE,CNCL)
IPL

#Odds ratio Logistic model 
PHIP<-function(lambda,xr,xf,er,ef,M1,M2){
  PHI1<-cdf.sl(ef,lambda)*(1-cdf.sl(er,lambda))
  PHI2<-cdf.sl(er,lambda)*(1-cdf.sl(ef,lambda))
  PHI<-PHI1/PHI2
  DPHI1<- - cdf.sl(ef,lambda)*pdf.sl(er,lambda)*xr+pdf.sl(ef,lambda)*(1-cdf.sl(er,lambda))*xf
  DPHI2<- - cdf.sl(er,lambda)*pdf.sl(ef,lambda)*xf+pdf.sl(er,lambda)*(1-cdf.sl(ef,lambda))*xr  
  DPHI<- 1/PHI2*(DPHI1-DPHI2*PHI)
  P<- -t(M2)%*%solve(M1)%*%DPHI
  return(P)
}

f=1 #Age 9.21
xf <- as.matrix(data[f,-(1:2)])
ef <- as.numeric(t(betas)%*%xf)
r=2 #Age 10.21
xr <- as.matrix(data[r,-(1:2)])
er <- as.numeric(t(betas)%*%xr)
PHIPL<-PHIP(lambda,xr,xf,er,ef,OIML,PML)
cuttofSEPHIPL<-mean(PHIPL)+2*sd(PHIPL)/sqrt(q)
plot(c(1:q),PHIPL,xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttofSEPHIPL,lty=2,lwd=1)

cuttofSEPHIPLabs<-mean(abs(PHIPL))+2*sd(abs(PHIPL))/sqrt(q)
plot(c(1:q),abs(PHIPL),xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttofSEPHIPLabs,lty=2,lwd=1)

IPPHIPL<-ipcnc(N,cuttofSEPHIPLabs,abs(PHIPL))
IPPHIPL



###########################################################################
###########################################################################
###########################################################################
###########################################################################

ls() # Print names of objects in R Environment
remove(list=ls()) # Remove all objects from R Environment

###########################################################################
###########################################################################
###########################################################################
###########################################################################

###########################################################################
#Probit model

#Observed Information Matrix
g<-function(z){
  g<-dnorm(z); #PDF standard normal
  return(g)
}

dg<-function(z){
  dg<--1/sqrt(2*pi)*z*exp(-(z^2)/2) #derivative PDF standard normal 
  return(dg)
}

pdf.sp<-function(z,lambda){2*g(z)*pnorm(lambda*z)}
cdf.sp<-function(z,lambda){integrate(pdf.sp,lower=-Inf, upper=z,
                                     lambda =lambda,stop.on.error=FALSE)$value}

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
  pii<-cdf.sp(etai,lambda)
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
  pii<-cdf.sp(etai,lambda)
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
  pii<-cdf.sp(etai,lambda)
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
g<-function(z){
  g<-dnorm(z); #PDF standard normal
  return(g)
}

cdf.p<-function(z){ #CDF standard normal
  cdf.p<-integrate(g,lower=-Inf,upper=z)
  return(cdf.p$value)
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
  pii <- cdf.p(etai)
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

#Function to calculate the influential points with CNC
ipcnc<-function(q,cuttof,L){
  cnc<-L
  index<-vector(length=q)
  c<-0
  for(i in 1:q){
    if(cnc[i]>cuttof){ 
      index[i]<-i
      c<-c+1
    }
  }
  L<-list(cnc,cuttof,c,index)
  return(L)
}

#Parameter estimates for Probit model 
fit.p<- glm(cbind(Mens,Entr-Mens)~Ida,family=binomial("probit"),data=df)

lambda <- 0
beta0.est.p <- fit.p$coefficients[1]
beta1.est.p <- fit.p$coefficients[2]
betas<-matrix(c(beta0.est.p,beta1.est.p),nrow=p,ncol=1,dimnames = NULL)
OIMP<-OIM(p,N,data,betas,lambda)
OIMP
PMP<-PM(p,N,data,betas)
PMP
CNCP<-CNC(OIMP,PMP)
q<-N
cuttofSE<-mean(CNCP)+2*sd(CNCP)/sqrt(q)

#plot(c(1:N),CNCP,ylim=c(0,1),xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
plot(c(1:N),CNCP,ylim=c(0,0.1),xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
abline(h=cuttofSE,lty=2,lwd=1)

IPP<-ipcnc(q,cuttofSE,CNCP)
IPP

#Odds ratio Probit model 
PHIP<-function(lambda,xr,xf,er,ef,M1,M2){
  PHI1<-cdf.sp(ef,lambda)*(1-cdf.sp(er,lambda))
  PHI2<-cdf.sp(er,lambda)*(1-cdf.sp(ef,lambda))
  PHI<-PHI1/PHI2
  DPHI1<- - cdf.sp(ef,lambda)*pdf.sp(er,lambda)*xr+pdf.sp(ef,lambda)*(1-cdf.sp(er,lambda))*xf
  DPHI2<- - cdf.sp(er,lambda)*pdf.sp(ef,lambda)*xf+pdf.sp(er,lambda)*(1-cdf.sp(ef,lambda))*xr  
  DPHI<- 1/PHI2*(DPHI1-DPHI2*PHI)
  P<- -t(M2)%*%solve(M1)%*%DPHI
  return(P)
}

f=1 #Age 9.21
xf <- as.matrix(data[f,-(1:2)])
ef <- as.numeric(t(betas)%*%xf)
r=2 #Age 10.21
xr <- as.matrix(data[r,-(1:2)])
er <- as.numeric(t(betas)%*%xr)
PHIPP<-PHIP(lambda,xr,xf,er,ef,OIMP,PMP)
cuttofSEPHIPP<-mean(PHIPP)+2*sd(PHIPP)/sqrt(q)
plot(c(1:q),PHIPP,xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttofSEPHIPP,lty=2,lwd=1)

cuttofSEPHIPPabs<-mean(abs(PHIPP))+2*sd(abs(PHIPP))/sqrt(q)
plot(c(1:q),abs(PHIPP),xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttofSEPHIPPabs,lty=2,lwd=1)

IPPHIPP<-ipcnc(N,cuttofSEPHIPPabs,abs(PHIPP))
IPPHIPP


###########################################################################
###########################################################################
###########################################################################
###########################################################################

ls() # Print names of objects in R Environment
remove(list=ls()) # Remove all objects from R Environment

###########################################################################
###########################################################################
###########################################################################
###########################################################################

###########################################################################
#t model
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
cdf.st<-function(z,lambda){integrate(pdf.st,lower=-Inf, upper=z,
                                     lambda =lambda,stop.on.error=FALSE)$value}

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

#Function to calculate the influential points with CNC
ipcnc<-function(q,cuttof,L){
  cnc<-L
  index<-vector(length=q)
  c<-0
  for(i in 1:q){
    if(cnc[i]>cuttof){ 
      index[i]<-i
      c<-c+1
    }
  }
  L<-list(cnc,cuttof,c,index)
  return(L)
}

#Parameter estimates for t model 
lambda <- 0
beta0.est.t <- -14.574919
beta1.est.t <- 1.121047
betas<-matrix(c(beta0.est.t,beta1.est.t),nrow=p,ncol=1,dimnames = NULL)
OIMT<-OIM(p,N,data,betas,lambda) 
OIMT

PMT<-PM(p,N,data,betas)
PMT

CNCT<-CNC(OIMT,PMT)
q<-N
cuttofSE<-mean(CNCT)+2*sd(CNCT)/sqrt(q)
plot(c(1:N),CNCT,ylim=c(0,1),xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
plot(c(1:N),CNCT,ylim=c(0,0.1),xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
abline(h=cuttofSE,lty=2,lwd=1)

IPT<-ipcnc(q,cuttofSE,CNCT)
IPT

#Odds ratio t model 
PHIP<-function(lambda,xr,xf,er,ef,M1,M2){
  PHI1<-cdf.st(ef,lambda)*(1-cdf.st(er,lambda))
  PHI2<-cdf.st(er,lambda)*(1-cdf.st(ef,lambda))
  PHI<-PHI1/PHI2
  DPHI1<- - cdf.st(ef,lambda)*pdf.st(er,lambda)*xr+pdf.st(ef,lambda)*(1-cdf.st(er,lambda))*xf
  DPHI2<- - cdf.st(er,lambda)*pdf.st(ef,lambda)*xf+pdf.st(er,lambda)*(1-cdf.st(ef,lambda))*xr  
  DPHI<- 1/PHI2*(DPHI1-DPHI2*PHI)
  P<- -t(M2)%*%solve(M1)%*%DPHI
  return(P)
}

f=1 #Age 9.21
xf <- as.matrix(data[f,-(1:2)])
etaf <- as.numeric(t(betas)%*%xf)
r=2 #Age 10.21
xr <- as.matrix(data[r,-(1:2)])
etar <- as.numeric(t(betas)%*%xr)
PHIPT<-PHIP(lambda,xr,xf,etar,etaf,OIMT,PMT)
cuttofSEPHIPT<-mean(PHIPT)+2*sd(PHIPT)/sqrt(q)
plot(c(1:q),PHIPT,xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttofSEPHIPT,lty=2,lwd=1)

cuttofSEPHIPTabs<-mean(abs(PHIPT))+2*sd(abs(PHIPT))/sqrt(q)
plot(c(1:q),abs(PHIPT),xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttofSEPHIPTabs,lty=2,lwd=1)

IPPHIPT<-ipcnc(N,cuttofSEPHIPTabs,abs(PHIPT))
IPPHIPT


