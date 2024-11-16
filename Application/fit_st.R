###############################################################################
#Parameter estimation: Age of menarche of Warsaw girls
###############################################################################

###############################################################################
#Read data
df<-read.table("girls.txt", sep="\t", header=TRUE)

#Get the data structure
str(df)

#Recognize each column of data
attach(df)

#Formula
form<-cbind(df$Mens,df$Entr-df$Mens)~df$Ida

################################################################################
#Skew-t model

#Function to calculate the CDF skew-t
nu<-5 #Suggest by Gómez et al. (2007) 
pdf.st<-function(z,lambda){2*dt(z,nu)*pnorm(lambda*z)}
cdf.st<-function(z,lambda){integrate(pdf.st,lower=-Inf, upper=z,
                                     lambda =lambda,stop.on.error=FALSE)$value}

#Negative log-likelihood function
negll.st <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.st(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val<- - sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  return(val)
}

library(optimx)
EPT<-c(-14.5749,1.1210)
theta.ini<-c(0,EPT)
pst<-length(theta.ini)
EPST<-matrix(NA,pst,1)
fit.st<-optimx(par=theta.ini,formula=form,data=df,fn=negll.st,method=c("BFGS","L-BFGS-B"))
coef.opt<-coef(fit.st)
EPST[1,1]<-coef.opt[1,1]#lambda BFGS
EPST[2,1]<-coef.opt[1,2]#beta0 BFGS
EPST[3,1]<-coef.opt[1,3]#beta1 BFGS
EPST

#P-value 
coef.opt<-coef(fit.st)
par.est.st<-coef.opt[1,]
p1 <- length(par.est.st)
hessian.st<- attributes(fit.st)$details["BFGS", "nhatend"][[1]]
oim.st<- hessian.st 
oim.st
std.error<-sqrt(diag(solve(oim.st)))
std.error
zvalue<-par.est.st/std.error
zvalue

pvalue<-numeric(length(zvalue))
for(j in 1:3){
  if(zvalue[j]>0){
    pvalue[j]<-pnorm(zvalue[j],lower.tail = FALSE)
  }
  else{
    pvalue[j]<-pnorm(zvalue[j],lower.tail = TRUE)
  }
}
2*pvalue

IC095<-matrix(NA,3,2)
for(j in 1:3){
  IC095[j,]<-par.est.st[j]+c(-1,1)*qnorm(0.975)*std.error[j]
}
IC095

#Fitted probabilities 
fit.prob.st<-function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.st(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  return(prob)
}

fit.prob<-fit.prob.st(par.est.st,form,df)
fit.prob

#Normalized randomized quantile residuals with envelopes 
#Based on Stasinopoulos et al. (2017) and Paula (2013).
a<-pbinom(df$Mens-1,df$Entr,fit.prob) 
b<-pbinom(df$Mens,df$Entr,fit.prob)
N<-nrow(df)
u<-runif(N,a,b)
res<-sort(qnorm(u))
#Envelopes
m<-100
e<-matrix(0,N,m);
yi<-numeric(N)
Ni<-df$Entr
xi<-df$Ida
EPST.envel<-matrix(NA,pst,1)
for(i in 1:m){
  for(j in 1:N){
    dif<-runif(Ni[j])-fit.prob[j];
    dif[dif>=0]<-0;
    dif[dif<0]<-1;
    yi[j]<-sum(dif);
  }
  theta.ini<-c(0,EPT)  
  dfi<-data.frame(yi,Ni,xi)
  formi<-cbind(dfi$yi,dfi$Ni-dfi$yi)~dfi$xi 
  fiti<-optimx(par=theta.ini,formula=formi,data=dfi,fn=negll.st,method="BFGS")  
  par.est.i<-coef(fiti)[1,]
  fit.prob.i<-fit.prob.st(par.est.i,formi,dfi)
  a<-pbinom(yi-1,Ni,fit.prob.i)
  b<-pbinom(yi,Ni,fit.prob.i)
  u<-runif(N,a,b)
  e[,i]<-sort(qnorm(u))
}

e1<-numeric(N)
e2<-numeric(N)

for(i in 1:N){
  eo<-sort(e[i,])
  e1[i]<-(eo[2]+eo[3])/2
  e2[i]<-(eo[97]+eo[98])/2
}

med<-apply(e,1,mean)
ran<-range(res,e1,e2)

par(pty="s")
qqnorm(res, main=" ",xlab="Theoretical Quantiles",
       ylab="Quantile residuals", ylim=ran, pch=16)
par(new=T)
qqnorm(e1, main=" ", axes=F,xlab="",ylab="",type="l",ylim=ran,lty=1)
par(new=T)
qqnorm(e2, main=" ", axes=F,xlab="",ylab="", type="l",ylim=ran,lty=1)
par(new=T)
qqnorm(med, main=" ", axes=F,xlab="", ylab="", type="l",ylim=ran,lty=2)

#Deviance
#Satured model
fac.sat<-factor(1:length(df$Mens)) 
df.sat<-df
df.sat$fac.sat<-fac.sat
form.sat<- cbind(df.sat$Mens,df.sat$Entr-df.sat$Mens)~df.sat$fac.sat

#Function to calculate the CDF t-student
nu<-5 #Suggest by Gómez et al. (2007) 
pdf.t<-function(z){dt(z,nu)}
cdf.t<-function(z){pt(z,nu)}

#Negative log-likelihood function
negll.t <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  beta <- par
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.t(eta[j,1])})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val<- - sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  return(val)
}

#options(warn = 1)
fit.l.sat<- glm(formula=form.sat,family=binomial(link="logit"),data=df.sat,maxit=1000)
theta.ini<-fit.l.sat$coefficients
pt<-length(theta.ini)
EPT<-matrix(NA,pt,1)
fit.t<-optimx(par=theta.ini,formula=form.sat,data=df,fn=negll.t,method=c("BFGS","L-BFGS-B"))
coef.opt<-coef(fit.t)
EPT[,1]<-coef.opt[1,]#betas BFGS
betas.ini.sat<-EPT

#Optimization
theta.ini.sat<-c(0,betas.ini.sat)
fit.st.sat<-optimx(par=theta.ini.sat,formula=form.sat,data=df.sat,fn=negll.st, method = "BFGS")
par.est.st.sat<-coef(fit.st.sat)
loglike.st.sat<- -negll.st(par.est.st.sat,form.sat,df.sat) 
loglike.st.sat

#Null model
form.null<- cbind(df$Mens,df$Entr-df$Mens)~1

#Negative log-likelihood function
negll.t <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  beta <- par
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.t(eta[j,1])})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val<- - sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  return(val)
}

fit.l.null<- glm(formula=form.null,family=binomial(link="logit"),data=df,maxit=1000)
theta.ini<-fit.l.null$coefficients
pt<-length(theta.ini)
EPT<-matrix(NA,pt,1)
fit.t<-optimx(par=theta.ini,formula=form.null,data=df,fn=negll.t,method=c("BFGS","L-BFGS-B"))
coef.opt<-coef(fit.t)
EPT[1,1]<-coef.opt[1,]#betas BFGS
betas.ini.null<-EPT

#Optimization
theta.ini.null<-c(0,betas.ini.null)
fit.st.null<-optimx(par=theta.ini.null,formula=form.null,data=df,fn=negll.st,method =c("BFGS","L-BFGS-B"))
fit.st.null 
coef.opt.null<-coef(fit.st.null)
par.est.st.null<-coef.opt.null[1,]
loglike.st.null<- -negll.st(par.est.st.null,form.null,df) 
loglike.st.null

#Proposed model
loglike.prop<- -negll.st(par.est.st,form,df) 
loglike.prop

#Residual deviance
res.dev<-2*(loglike.st.sat - loglike.prop) 
res.dev

#Null deviance
res.null<-2*(loglike.st.sat - loglike.st.null) 
res.null

#AIC
AIC.st <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.st(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- 2*p1-2*val1
  return(val2)
}

AIC.st(par.est.st,form,df) 

AICc.st <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.st(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- (2*p1-2*val1)+((2*p1^(2)+2*p1)/(length(y)-p1-1))
  return(val2)
}

AICc.st(par.est.st,form,df) 

#BIC
BIC.st <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.st(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- 2*log(length(y))-2*val1
  return(val2)
}

BIC.st(par.est.st,form,df) 


