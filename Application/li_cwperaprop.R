###############################################################################
#Local influence under case-weight appropriate perturbation: 
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

#Skew-logistic model

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
cdf.sl<-function(z,lambda){integrate(pdf.sl,lower=-Inf,upper=z,
                                     lambda=lambda,stop.on.error=FALSE)$value}

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
  pii<-cdf.sl(etai,lambda)
  pii<-ifelse(pii>=0.9999,0.9999, ifelse(pii<=0.00001,0.00001,pii))
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


#Perturbation matrix under case-weight perturbation scheme
pdf.sl<-function(z,lambda){2*dlogis(z,log=FALSE)*pnorm(lambda*z)}
cdf.sl<-function(z,lambda){integrate(pdf.sl,lower=-Inf,upper=z,
                                     lambda=lambda,stop.on.error=FALSE)$value}

litheta<-function(Ni,yi,etai,lambda){
  litheta<-lchoose(Ni,yi)+yi*log(cdf.sl(etai,lambda))+(Ni-yi)*log(1-cdf.sl(etai,lambda))
  return(litheta)
}

dsc<-function(Ni,yi,etai,lambda){
  s<-0 
  yseq<-seq(0,Ni,1)
  for(j in 1:(Ni+1)){
    yaux<-yseq[j]
    a<-litheta(Ni,yaux,etai,lambda)^2*exp(litheta(Ni,yaux,etai,lambda))
    s<-s+a
  }
  return(s)
}

dpc<-function(Ni,yi,etai,lambda){
  s<-0 
  yseq<-seq(0,Ni,1)
  for(j in 1:(Ni+1)){
    yaux<-yseq[j]
    a<-litheta(Ni,yaux,etai,lambda)*exp(litheta(Ni,yaux,etai,lambda))
    s<-s+a
  }
  return(s)
}

Vi0<-function(Ni,yi,etai,lambda){
  Vi0<-dsc(Ni,yi,etai,lambda)-dpc(Ni,yi,etai,lambda)^2
  return(Vi0)
}

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

dlithetalambda<-function(Ni,yi,etai,lambda){
  pii<-cdf.sl(etai,lambda)
  piil<-F3(lambda,etai)
  si<-piil/(pii*(1-pii))
  dlithetalambda<-si*(yi-Ni*pii)
  return(dlithetalambda)
}

F6<-function(lambda,eta){
  as.numeric(2*g(eta)*pnorm(lambda*eta,0,1))
}

dlithetabeta<-function(Ni,yi,xi,etai,lambda){
  pii<-cdf.sl(etai,lambda)
  wii<-F6(lambda,etai)/(pii*(1-pii))
  dlithetabeta<-wii*(yi-Ni*pii)*xi
  return(dlithetabeta)
}

dlitheta<-function(Ni,yi,xi,etai,lambda){
  dlitheta<-matrix(0,p+1,1)
  L1<-as.matrix(dlithetalambda(Ni,yi,etai,lambda))
  L2<-as.matrix(dlithetabeta(Ni,yi,xi,etai,lambda))
  dlitheta<-rbind(L1,L2)
  return(dlitheta)
}

dVi0theta<-function(Ni,yi,xi,etai,lambda){
  s<-0 
  yseq<-seq(0,Ni,1)
  for(j in 1:(Ni+1)){
    yaux<-yseq[j]
    l1<-litheta(Ni,yaux,etai,lambda)
    l2<-dpc(Ni,yi,etai,lambda)
    l3<-dlitheta(Ni,yi,xi,etai,lambda)
    a<-((2+l1)*l1-2*l2*(1+l1))*exp(l1)*l3
    s<-s+a
  }
  return(s)
}

PMi<-function(i,data,betas,lambda){
  yi <- data[i,1]
  Ni <- data[i,2]
  xi <- as.matrix(data[i,-(1:2)])
  etai <- as.numeric(t(betas)%*%xi)
  l1<-Vi0(Ni,yi,etai,lambda)
  l2<-dlitheta(Ni,yi,xi,etai,lambda)
  l3<-litheta(Ni,yi,etai,lambda)
  l4<-dVi0theta(Ni,yi,xi,etai,lambda)
  PMi<-l1^(-1/2)*(l2-(1/2)*l1^(-1)*l3*l4)
  return(PMi);
}

PM<-function(p,N,data,betas,lambda){
  PM<-PMi(1,data,betas,lambda)
  for(i in 2:N){
    PM<-cbind(PM,PMi(i,data,betas,lambda))
  }
  PM<-matrix(PM,nrow=p+1,ncol=N,dimnames = NULL)
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

#Parameter estimates for skew-logistic model 
lambda<-0.506505
beta0.est<- -17.361260
beta1.est<-1.406158
betas<-matrix(c(beta0.est,beta1.est),nrow=p,ncol=1,dimnames = NULL)
OIMSL<-OIM(p,N,data,betas,lambda)
PMSL<-PM(p,N,data,betas,lambda)
CNCSL<-CNC(OIMSL,PMSL)
q<-N
cuttofSE<-mean(CNCSL)+2*sd(CNCSL)/sqrt(q)
plot(c(1:N),CNCSL,ylim=c(0,1),xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
#plot(c(1:N),CNCSL,xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
abline(h=cuttofSE,lty=2,lwd=1)

IPSL<-ipcnc(q,cuttofSE,CNCSL)
IPSL

#Odds ratio Skew-logistic model
PHIP<-function(lambda,xr,xf,er,ef,M1,M2){
  PHI1<-cdf.sl(ef,lambda)*(1-cdf.sl(er,lambda))
  PHI2<-cdf.sl(er,lambda)*(1-cdf.sl(ef,lambda))
  PHI<-PHI1/PHI2
  DPHI1betas<- - cdf.sl(ef,lambda)*pdf.sl(er,lambda)*xr+pdf.sl(ef,lambda)*(1-cdf.sl(er,lambda))*xf
  DPHI2betas<- - cdf.sl(er,lambda)*pdf.sl(ef,lambda)*xf+pdf.sl(er,lambda)*(1-cdf.sl(ef,lambda))*xr  
  DPHIbetas<- as.matrix(1/PHI2*(DPHI1betas-DPHI2betas*PHI))
  DPHI1lambda<-F3(lambda,ef)*(1-cdf.sl(er,lambda))-cdf.sl(ef,lambda)*F3(lambda,er)
  DPHI2lambda<-F3(lambda,er)*(1-cdf.sl(ef,lambda))-cdf.sl(er,lambda)*F3(lambda,ef)
  DPHIlambda<-as.matrix(1/PHI2*(DPHI1lambda-DPHI2lambda*PHI))
  DPHI<-rbind(DPHIlambda,DPHIbetas)
  P<- -t(M2)%*%solve(M1)%*%DPHI
  return(P)
}

r=1 #Age 9.21
xr <- as.matrix(data[r,-(1:2)])
er <- as.numeric(t(betas)%*%xr)
f=2 #Age 10.21
xf <- as.matrix(data[f,-(1:2)])
ef <- as.numeric(t(betas)%*%xf)
PHIPSL<-PHIP(lambda,xr,xf,er,ef,OIMSL,PMSL)
q<-N
cuttof<-mean(abs(PHIPSL))+2*sd(abs(PHIPSL))/sqrt(q)

plot(c(1:q),abs(PHIPSL),xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
plot(c(1:q),PHIPSL,xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttof,lty=2,lwd=1)

IPPHIPSL<-ipcnc(N,cuttof,abs(PHIPSL))
IPPHIPSL



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
#Skew-probit model

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
cdf.sp<-function(z,lambda){integrate(pdf.sp,lower=-Inf,upper=z,
                                     lambda=lambda,stop.on.error=FALSE)$value}

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

#Perturbation matrix under case-weight perturbation scheme
pdf.sp<-function(z,lambda){2*g(z)*pnorm(lambda*z)}
cdf.sp<-function(z,lambda){integrate(pdf.sp,lower=-Inf, upper=z,
                                     lambda=lambda,stop.on.error=FALSE)$value}

litheta<-function(Ni,yi,etai,lambda){
  litheta<-lchoose(Ni,yi)+yi*log(cdf.sp(etai,lambda))+(Ni-yi)*log(1-cdf.sp(etai,lambda))
  return(litheta)
}

dsc<-function(Ni,yi,etai,lambda){
  s<-0 
  yseq<-seq(0,Ni,1)
  for(j in 1:(Ni+1)){
    yaux<-yseq[j]
    a<-litheta(Ni,yaux,etai,lambda)^2*exp(litheta(Ni,yaux,etai,lambda))
    s<-s+a
  }
  return(s)
}

dpc<-function(Ni,yi,etai,lambda){
  s<-0 
  yseq<-seq(0,Ni,1)
  for(j in 1:(Ni+1)){
    yaux<-yseq[j]
    a<-litheta(Ni,yaux,etai,lambda)*exp(litheta(Ni,yaux,etai,lambda))
    s<-s+a
  }
  return(s)
}

Vi0<-function(Ni,yi,etai,lambda){
  Vi0<-dsc(Ni,yi,etai,lambda)-dpc(Ni,yi,etai,lambda)^2
  return(Vi0)
}

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

dlithetalambda<-function(Ni,yi,etai,lambda){
  pii<-cdf.sp(etai,lambda)
  piil<-F3(lambda,etai)
  si<-piil/(pii*(1-pii))
  dlithetalambda<-si*(yi-Ni*pii)
  return(dlithetalambda)
}


F6<-function(lambda,eta){
  as.numeric(2*g(eta)*pnorm(lambda*eta,0,1))
}

dlithetabeta<-function(Ni,yi,xi,etai,lambda){
  pii<-cdf.sp(etai,lambda)
  wii<-F6(lambda,etai)/(pii*(1-pii))
  dlithetabeta<-wii*(yi-Ni*pii)*xi
  return(dlithetabeta)
}

dlitheta<-function(Ni,yi,xi,etai,lambda){
  if(lambda!=0){
    dlitheta<-matrix(0,p+1,1)
    L1<-as.matrix(dlithetalambda(Ni,yi,etai,lambda))
    L2<-as.matrix(dlithetabeta(Ni,yi,xi,etai,lambda))
    dlitheta<-rbind(L1,L2)
  }
  else{
    dlitheta<-matrix(0,p,1)
    L2<-as.matrix(dlithetabeta(Ni,yi,xi,etai,lambda))
    dlitheta<-L2
  }
  return(dlitheta)
}

dVi0theta<-function(Ni,yi,xi,etai,lambda){
  s<-0 
  yseq<-seq(0,Ni,1)
  for(j in 1:(Ni+1)){
    yaux<-yseq[j]
    l1<-litheta(Ni,yaux,etai,lambda)
    l2<-dpc(Ni,yi,etai,lambda)
    l3<-dlitheta(Ni,yi,xi,etai,lambda)
    a<-((2+l1)*l1-2*l2*(1+l1))*exp(l1)*l3
    s<-s+a
  }
  return(s)
}

PMi<-function(i,data,betas,lambda){
  yi <- data[i,1]
  Ni <- data[i,2]
  xi <- as.matrix(data[i,-(1:2)])
  etai <- as.numeric(t(betas)%*%xi)
  l1<-Vi0(Ni,yi,etai,lambda)
  l2<-dlitheta(Ni,yi,xi,etai,lambda)
  l3<-litheta(Ni,yi,etai,lambda)
  l4<-dVi0theta(Ni,yi,xi,etai,lambda)
  PMi<-l1^(-1/2)*(l2-(1/2)*l1^(-1)*l3*l4)
  return(PMi);
}

PM<-function(p,N,data,betas,lambda){
  PM<-PMi(1,data,betas,lambda)
  for(i in 2:N){
    PM<-cbind(PM,PMi(i,data,betas,lambda))
  }
  #PM<-matrix(PM,nrow=p+1,ncol=N,dimnames = NULL)
  return(PM)
}

#Parameter estimates for Skew-probit model 
lambda<- -0.002190207
beta0.est<- -11.817312719
beta1.est<-  0.907556599
betas<-matrix(c(beta0.est,beta1.est),nrow=p,ncol=1,dimnames = NULL)
OIMSP<-OIM(p,N,data,betas,lambda)
PMSP<-PM(p,N,data,betas,lambda)
CNCSP<-CNC(OIMSP,PMSP)
q<-N
cuttofSE<-mean(CNCSP)+2*sd(CNCSP)/sqrt(q)
plot(c(1:N),CNCSP,ylim=c(0,1),xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
#plot(c(1:N),CNCSP,xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
abline(h=cuttofSE,lty=2,lwd=1)

IPSP<-ipcnc(q,cuttofSE,CNCSP)
IPSP

#Odds ratio Skew-probit model
PHIP<-function(lambda,xr,xf,er,ef,M1,M2){
  PHI1<-cdf.sp(ef,lambda)*(1-cdf.sp(er,lambda))
  PHI2<-cdf.sp(er,lambda)*(1-cdf.sp(ef,lambda))
  PHI<-PHI1/PHI2
  DPHI1betas<- - cdf.sp(ef,lambda)*pdf.sp(er,lambda)*xr+pdf.sp(ef,lambda)*(1-cdf.sp(er,lambda))*xf
  DPHI2betas<- - cdf.sp(er,lambda)*pdf.sp(ef,lambda)*xf+pdf.sp(er,lambda)*(1-cdf.sp(ef,lambda))*xr  
  DPHIbetas<- as.matrix(1/PHI2*(DPHI1betas-DPHI2betas*PHI))
  DPHI1lambda<-F3(lambda,ef)*(1-cdf.sp(er,lambda))-cdf.sp(ef,lambda)*F3(lambda,er)
  DPHI2lambda<-F3(lambda,er)*(1-cdf.sp(ef,lambda))-cdf.sp(er,lambda)*F3(lambda,ef)
  DPHIlambda<-as.matrix(1/PHI2*(DPHI1lambda-DPHI2lambda*PHI))
  DPHI<-rbind(DPHIlambda,DPHIbetas)
  P<- -t(M2)%*%solve(M1)%*%DPHI
  return(P)
}

f=1 #Age 9.21
xf <- as.matrix(data[f,-(1:2)])
ef <- as.numeric(t(betas)%*%xf)
r=2 #Age 10.21
xr <- as.matrix(data[r,-(1:2)])
er <- as.numeric(t(betas)%*%xr)
PHIPSP<-PHIP(lambda,xr,xf,er,ef,OIMSP,PMSP)
q<-N
cuttofabs<-mean(abs(PHIPSP))+2*sd(abs(PHIPSP))/sqrt(q)
cuttof<-mean(PHIPSP)+2*sd(PHIPSP)/sqrt(q)

plot(c(1:q),abs(PHIPSP),xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttofabs,lty=2,lwd=1)

plot(c(1:q),PHIPSP,xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttof,lty=2,lwd=1)

IPPHIPSP<-ipcnc(N,cuttof,abs(PHIPSP))
IPPHIPSP


#Parameter estimates for Probit model 
lambda<- 0
beta0.est<- -11.8156
beta1.est<-  0.9075
betas<-matrix(c(beta0.est,beta1.est),nrow=p,ncol=1,dimnames = NULL)
OIMP<-OIM(p,N,data,betas,lambda)
PMP<-PM(p,N,data,betas,lambda)
CNCP<-CNC(OIMP,PMP)
q<-N
cuttofSE<-mean(CNCP)+2*sd(CNCP)/sqrt(q)
plot(c(1:N),CNCP,ylim=c(0,1),xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
#plot(c(1:N),CNCSP,xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
abline(h=cuttofSE,lty=2,lwd=1)

IPP<-ipcnc(q,cuttofSE,CNCP)
IPP

#Odds ratio Probit model
PHIP<-function(lambda,xr,xf,er,ef,M1,M2){
  PHI1<-cdf.sp(ef,lambda)*(1-cdf.sp(er,lambda))
  PHI2<-cdf.sp(er,lambda)*(1-cdf.sp(ef,lambda))
  PHI<-PHI1/PHI2
  DPHI1betas<- - cdf.sp(ef,lambda)*pdf.sp(er,lambda)*xr+pdf.sp(ef,lambda)*(1-cdf.sp(er,lambda))*xf
  DPHI2betas<- - cdf.sp(er,lambda)*pdf.sp(ef,lambda)*xf+pdf.sp(er,lambda)*(1-cdf.sp(ef,lambda))*xr  
  DPHIbetas<- as.matrix(1/PHI2*(DPHI1betas-DPHI2betas*PHI))
  DPHI<-DPHIbetas
  P<- -t(M2)%*%solve(M1)%*%DPHI
  return(P)
}

r=1 #Age 9.21
xr <- as.matrix(data[r,-(1:2)])
er <- as.numeric(t(betas)%*%xr)
f=2 #Age 10.21
xf <- as.matrix(data[f,-(1:2)])
ef <- as.numeric(t(betas)%*%xf)
PHIPP<-PHIP(lambda,xr,xf,er,ef,OIMP,PMP)
q<-N
cuttofabs<-mean(abs(PHIPP))+2*sd(abs(PHIPP))/sqrt(q)
cuttof<-mean(PHIPP)+2*sd(PHIPP)/sqrt(q)

plot(c(1:q),abs(PHIPP),xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttofabs,lty=2,lwd=1)

plot(c(1:q),PHIPP,xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttof,lty=2,lwd=1)

IPPHIPP<-ipcnc(N,cuttofabs,abs(PHIPP))
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

#Observed Information Matrix
g<-function(z,nu){
  g<-dt(z,nu); #PDF t
  return(g)
}

dg<-function(z,nu){
  dg<-gamma((nu+1)/2)/(sqrt(pi*nu)*gamma(nu/2))*(-(1/2)*(nu+1))*((2*z)/nu)*(1+z^2/nu)^(-(1/2)*(nu+1)-1) #derivative PDF t 
  return(dg)
}

pdf.st<-function(z,nu,lambda){2*g(z,nu)*pnorm(lambda*z)}
cdf.st<-function(z,nu,lambda){integrate(pdf.st,lower=-Inf,upper=z,nu=nu,
                                        lambda=lambda,stop.on.error=FALSE)$value}


#F2<-function(z,nu,lambda){2*z*g(z,nu)*dnorm(lambda*z,0,1)}
#F3<-function(z,nu,lambda){integrate(F2,lower=-Inf,upper=z,nu=nu,
#                                        lambda=lambda,stop.on.error=FALSE)$value}

F2<-function(nu,lambda){
  f2<-function(z){
    f2<-2*z*g(z,nu)*dnorm(lambda*z,0,1)
  }
  return(f2)
}

F3<-function(nu,lambda,eta){ 
  F3<-integrate(F2(nu,lambda),lower=-Inf,upper=eta)
  return(F3$value) 
}

F4<-function(nu,lambda){
  f4<-function(z){
    f4<- -2*lambda*z^3*g(z,nu)*dnorm(lambda*z,0,1)
  }
  return(f4)
}

F5<-function(nu,lambda,eta){ 
  F5<-integrate(F4(nu,lambda),lower=-Inf,upper=eta)
  return(F5$value) 
}

Llli<-function(Ni,yi,etai,nu,lambda){
  pii<-cdf.st(etai,nu,lambda)
  piil<-F3(nu,lambda,etai)
  si<-piil/(pii*(1-pii))
  vi<-pii*(1-pii)
  piill<-F5(nu,lambda,etai)
  sil<-vi^(-1)*(piill-(1-2*pii)*piil^2*vi^(-1))
  Llli<-(yi-Ni*pii)*sil-Ni*si*piil
  return(Llli)
}

Lll<-function(N,data,nu,betas,lambda){
  s<-matrix(0,1,1)
  for(i in 1:N){
    yi<-data[i,1]
    Ni<-data[i,2]
    xi<-data[i,-(1:2)]
    etai<-t(betas)%*%xi
    a<-Llli(Ni,yi,etai,nu,lambda)
    s<-s+a
  }
  return(s)
}

F6<-function(nu,lambda,eta){
  as.numeric(2*g(eta,nu)*pnorm(lambda*eta,0,1))
}

F7<-function(nu,lambda,eta){
  as.numeric(2*eta*g(eta,nu)*dnorm(lambda*eta,0,1))
}

Llbi<-function(Ni,yi,xi,etai,nu,lambda){
  pii<-cdf.st(etai,nu,lambda)
  vi<-pii*(1-pii)
  piil<-F3(nu,lambda,etai)
  si<-piil/(pii*(1-pii))
  piib<-as.numeric(F6(nu,lambda,etai))*xi
  piilb<-F7(nu,lambda,etai)*xi
  sib<-vi^(-1)*(piilb-(1-2*pii)*piil*piib*vi^(-1))
  Llbi<-(yi-Ni*pii)*sib-Ni*si*piib
  return(Llbi)
}

Llb<-function(p,N,data,nu,betas,lambda){
  s<-matrix(0,1,p)
  for(i in 1:N){
    yi<-data[i,1]
    Ni<-data[i,2]
    xi<-data[i,-(1:2)]
    etai<-t(betas)%*%xi
    a<-Llbi(Ni,yi,xi,etai,nu,lambda)
    s<-s+a
  }
  return(s)
}

F8<-function(nu,lambda,eta){
  as.numeric(2*(g(eta,nu)*dnorm(lambda*eta,0,1)*lambda+pnorm(lambda*eta,0,1)*dg(eta,nu)))
}

Lbbi<-function(Ni,yi,xi,etai,nu,lambda){
  pii<-cdf.st(etai,nu,lambda)
  vi<-pii*(1-pii)
  wii<-F6(nu,lambda,etai)/(pii*(1-pii))
  vii<-vi^(-1)*(yi-Ni*pii)*(F8(nu,lambda,etai)-F6(nu,lambda,etai)^2*vi^(-1)*(1-2*pii))-F6(nu,lambda,etai)*Ni*wii
  Lbbi<-as.numeric(vii)*(xi%*%t(xi))
  return(Lbbi)
}

Lbb<-function(p,N,data,nu,betas,lambda){
  s<-matrix(0,p,p)
  for(i in 1:N){
    yi<-data[i,1]
    Ni<-data[i,2]
    xi<-as.matrix(data[i,-(1:2)])
    etai<-t(betas)%*%xi
    a<-Lbbi(Ni,yi,xi,etai,nu,lambda)
    s<-s+a
  }
  return(s)
}

OIM<-function(p,N,data,nu,betas,lambda){
  if(lambda!=0){
    L11 <- -Lll(N,data,nu,betas,lambda)
    L12 <- -Llb(p,N,data,nu,betas,lambda)
    L21 <- t(L12)
    L22 <- -Lbb(p,N,data,nu,betas,lambda)
    L1 <- as.matrix(cbind(L11,L12))
    L2 <- as.matrix(cbind(L21,L22))
    OIM<-rbind(L1,L2)
  }
  else{
    L22 <- -Lbb(p,N,data,nu,betas,lambda)
    OIM<-L22
  }
  return(OIM)
}

#Perturbation matrix under case-weight perturbation scheme
pdf.st<-function(z,nu,lambda){2*g(z,nu)*pnorm(lambda*z)}
cdf.st<-function(z,nu,lambda){integrate(pdf.st,lower=-Inf,upper=z,nu=nu,
                                        lambda=lambda,stop.on.error=FALSE)$value}

litheta<-function(Ni,yi,etai,nu,lambda){
  litheta<-lchoose(Ni,yi)+yi*log(cdf.st(etai,nu,lambda))+(Ni-yi)*log(1-cdf.st(etai,nu,lambda))
  return(litheta)
}

dsc<-function(Ni,yi,etai,nu,lambda){
  s<-0 
  yseq<-seq(0,Ni,1)
  for(j in 1:(Ni+1)){
    yaux<-yseq[j]
    a<-litheta(Ni,yaux,etai,nu,lambda)^2*exp(litheta(Ni,yaux,etai,nu,lambda))
    s<-s+a
  }
  return(s)
}

dpc<-function(Ni,yi,etai,nu,lambda){
  s<-0 
  yseq<-seq(0,Ni,1)
  for(j in 1:(Ni+1)){
    yaux<-yseq[j]
    a<-litheta(Ni,yaux,etai,nu,lambda)*exp(litheta(Ni,yaux,etai,nu,lambda))
    s<-s+a
  }
  return(s)
}

Vi0<-function(Ni,yi,etai,nu,lambda){
  Vi0<-dsc(Ni,yi,etai,nu,lambda)-dpc(Ni,yi,etai,nu,lambda)^2
  return(Vi0)
}

F2<-function(nu,lambda){
  f2<-function(z){
    f2<-2*z*g(z,nu)*dnorm(lambda*z,0,1)
  }
  return(f2)
}

F3<-function(nu,lambda,eta){ 
  F3<-integrate(F2(nu,lambda),lower=-Inf,upper=eta)
  return(F3$value) 
}

dlithetalambda<-function(Ni,yi,etai,nu,lambda){
  pii<-cdf.st(etai,nu,lambda)
  piil<-F3(nu,lambda,etai)
  si<-piil/(pii*(1-pii))
  dlithetalambda<-si*(yi-Ni*pii)
  return(dlithetalambda)
}

F6<-function(nu,lambda,eta){
  2*g(eta,nu)*pnorm(lambda*eta,0,1)
}

dlithetabeta<-function(Ni,yi,xi,etai,nu,lambda){
  pii<-cdf.st(etai,nu,lambda)
  wii<-F6(nu,lambda,etai)/(pii*(1-pii))
  dlithetabeta<-wii*(yi-Ni*pii)*xi
  return(dlithetabeta)
}

dlitheta<-function(Ni,yi,xi,etai,nu,lambda){
  dlitheta<-matrix(0,p+1,1)
  L1<-as.matrix(dlithetalambda(Ni,yi,etai,nu,lambda))
  L2<-as.matrix(dlithetabeta(Ni,yi,xi,etai,nu,lambda))
  dlitheta<-rbind(L1,L2)
  return(dlitheta)
}

dVi0theta<-function(Ni,yi,xi,etai,nu,lambda){
  s<-0 
  yseq<-seq(0,Ni,1)
  for(j in 1:(Ni+1)){
    yaux<-yseq[j]
    l1<-litheta(Ni,yaux,etai,nu,lambda)
    l2<-dpc(Ni,yi,etai,nu,lambda)
    l3<-dlitheta(Ni,yi,xi,etai,nu,lambda)
    a<-((2+l1)*l1-2*l2*(1+l1))*exp(l1)*l3
    s<-s+a
  }
  return(s)
}

PMi<-function(i,data,nu,betas,lambda){
  yi <- data[i,1]
  Ni <- data[i,2]
  xi <- as.matrix(data[i,-(1:2)])
  etai <- as.numeric(t(betas)%*%xi)
  l1<-Vi0(Ni,yi,etai,nu,lambda)
  l2<-dlitheta(Ni,yi,xi,etai,nu,lambda)
  l3<-litheta(Ni,yi,etai,nu,lambda)
  l4<-dVi0theta(Ni,yi,xi,etai,nu,lambda)
  PMi<-l1^(-1/2)*(l2-(1/2)*l1^(-1)*l3*l4)
  return(PMi);
}

PM<-function(p,N,data,nu,betas,lambda){
  PM<-PMi(1,data,nu,betas,lambda)
  for(i in 2:N){
    PM<-cbind(PM,PMi(i,data,nu,betas,lambda))
  }
  PM<-matrix(PM,nrow=p+1,ncol=N,dimnames = NULL)
  return(PM)
}

#Parameter estimates for Skew-t model 
lambda<-  0.5799017
beta0.est<- -12.6600775
beta1.est<-  1.0158864
betas<-matrix(c(beta0.est,beta1.est),nrow=p,ncol=1,dimnames = NULL)
nu<-5

OIMST<-OIM(p,N,data,nu,betas,lambda) 
PMST<-PM(p,N,data,nu,betas,lambda)
CNCST<-CNC(OIMST,PMST)
q<-N
cuttofSE<-mean(CNCST)+2*sd(CNCST)/sqrt(q)
plot(c(1:N),CNCST,ylim=c(0,1),xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
plot(c(1:N),CNCST,xlab = "Index",ylab="Conformal normal curvature",col="black",pch=1,lwd=1)
abline(h=cuttofSE,lty=2,lwd=1)

IPST<-ipcnc(q,cuttofSE,CNCST)
IPST

#Odds ratio Skew-t model
PHIP<-function(nu,lambda,xr,xf,er,ef,M1,M2){
  PHI1<-cdf.st(ef,nu,lambda)*(1-cdf.st(er,nu,lambda))
  PHI2<-cdf.st(er,nu,lambda)*(1-cdf.st(ef,nu,lambda))
  PHI<-PHI1/PHI2
  DPHI1betas<- - cdf.st(ef,nu,lambda)*pdf.st(er,nu,lambda)*xr+pdf.st(ef,nu,lambda)*(1-cdf.st(er,nu,lambda))*xf
  DPHI2betas<- - cdf.st(er,nu,lambda)*pdf.st(ef,nu,lambda)*xf+pdf.st(er,nu,lambda)*(1-cdf.st(ef,nu,lambda))*xr  
  DPHIbetas<- as.matrix(1/PHI2*(DPHI1betas-DPHI2betas*PHI))
  DPHI1lambda<-F3(nu,lambda,ef)*(1-cdf.st(er,nu,lambda))-cdf.st(ef,nu,lambda)*F3(nu,lambda,er)
  DPHI2lambda<-F3(nu,lambda,er)*(1-cdf.st(ef,nu,lambda))-cdf.st(er,nu,lambda)*F3(nu,lambda,ef)
  DPHIlambda<-as.matrix(1/PHI2*(DPHI1lambda-DPHI2lambda*PHI))
  DPHI<-rbind(DPHIlambda,DPHIbetas)
  P<- -t(M2)%*%solve(M1)%*%DPHI
  return(P)
}

r=1 #Age 9.21
xr <- as.matrix(data[r,-(1:2)])
er <- as.numeric(t(betas)%*%xr)
f=2 #Age 10.21
xf <- as.matrix(data[f,-(1:2)])
ef <- as.numeric(t(betas)%*%xf)
PHIPST<-PHIP(nu,lambda,xr,xf,er,ef,OIMST,PMST)
q<-N
cuttofabs<-mean(abs(PHIPST))+2*sd(abs(PHIPST))/sqrt(q)
cuttof<-mean(PHIPST)+2*sd(PHIPST)/sqrt(q)

plot(c(1:q),abs(PHIPST),xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttofabs,lty=2,lwd=1)

plot(c(1:q),PHIPST,xlab = "Index",ylab="Slope",col="black",pch=1,lwd=1)
abline(h=cuttof,lty=2,lwd=1)

IPPHIPST<-ipcnc(N,cuttofabs,abs(PHIPST))
IPPHIPST



