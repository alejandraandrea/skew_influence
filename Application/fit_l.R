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
#Logistic model

#Function to calculate the CDF Logistic
pdf.l<-function(z){dlogis(z,log=FALSE)}
cdf.l<-function(z){plogis(z,log=FALSE)}

#Negative log-likelihood function
negll.l <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  beta <- par
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.l(eta[j,1])})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val<- - sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  return(val)
}

library(optimx)
theta.ini<-c(-21,1.6) 
pl<-length(theta.ini)
EPL<-matrix(NA,pl,1)
fit.l<-optimx(par=theta.ini,formula=form,data=df,fn=negll.l,method=c("BFGS","L-BFGS-B"))
coef.opt<-coef(fit.l)
EPL[1,1]<-coef.opt[1,1]#beta0 BFGS
EPL[2,1]<-coef.opt[1,2]#beta1 BFGS
EPL

#P-value 
coef.opt<-coef(fit.l)
par.est.l<-coef.opt[1,]
p1<-length(par.est.l)
hessian.l<- attributes(fit.l)$details["BFGS", "nhatend"][[1]]
oim.l<- hessian.l 
oim.l
std.error<-sqrt(diag(solve(oim.l)))
std.error

zvalue<-par.est.l/std.error
pvalue<-numeric(length(zvalue))

for(j in 1:2){
  if(zvalue[j]>0){
    pvalue[j]<-pnorm(zvalue[j],lower.tail = FALSE)
  }
  else{
    pvalue[j]<-pnorm(zvalue[j],lower.tail = TRUE)
  }
}
2*pvalue

#IC
IC095<-matrix(NA,2,2)
for(j in 1:2){
  IC095[j,]<-par.est.l[j]+c(-1,1)*qnorm(0.975)*std.error[j]
}
IC095

#Fitted probabilities
fit.prob.l<-function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  beta <- par
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.l(eta[j,1])})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  return(prob)
}

fit.prob<-fit.prob.l(par.est.l,form,df)
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
EPL.envel<-matrix(NA,pl,1)
for(i in 1:m){
  for(j in 1:N){
    dif<-runif(Ni[j])-fit.prob[j];
    dif[dif>=0]<-0;
    dif[dif<0]<-1;
    yi[j]<-sum(dif);
  }
  theta.ini<-c(-21,1.6) 
  dfi<-data.frame(yi,Ni,xi)
  formi<-cbind(dfi$yi,dfi$Ni-dfi$yi)~dfi$xi
  fiti<-optimx(par=theta.ini,formula=formi,data=dfi,fn=negll.l,method="BFGS")  
  par.est.i<-coef(fiti)[1,]
  fit.prob.i<-fit.prob.l(par.est.i,formi,dfi)
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
form.sat<- cbind(df$Mens,df$Entr-df$Mens)~df.sat$fac.sat
fit.l.sat<- glm(formula=form.sat,family=binomial(link="logit"),data=df.sat,maxit=1000)
betas.ini.sat<-fit.l.sat$coefficients

#Optimization
theta.ini.sat<-betas.ini.sat
fit.l.sat<-optimx(par=theta.ini.sat,formula=form.sat,data=df.sat,fn=negll.l, method = "BFGS")
par.est.l.sat<-coef(fit.l.sat)[1,] 
loglike.l.sat<- -negll.l(par.est.l.sat,form.sat,df.sat) 
loglike.l.sat

#Null model
form.null<- cbind(df$Mens,df$Entr-df$Mens)~1
fit.l.null<- glm(formula=form.null,family=binomial(link="logit"),data=df)
betas.ini.null<-fit.l.null$coefficients

#Optimization
theta.ini.null<-betas.ini.null
fit.l.null<-optimx(par=theta.ini.null,formula=form.null,data=df,fn=negll.l,method =c("BFGS","L-BFGS-B"))
fit.l.null 
par.est.l.null<-coef(fit.l.null)[1,]
loglike.l.null<- -negll.l(par.est.l.null,form.null,df) 
loglike.l.null

#Proposed model
loglike.prop<- -negll.l(par.est.l,form,df) 
loglike.prop

#Residual deviance
res.dev<-2*(loglike.l.sat - loglike.prop) 
res.dev

#Null deviance
res.null<-2*(loglike.l.sat - loglike.l.null) 
res.null

#AIC
AIC.l <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  beta <- par
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.l(eta[j,1])})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- 2*p1-2*val1
  return(val2)
}

AIC.l(par.est.l,form,df) 

AICc.l <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  beta <- par
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.l(eta[j,1])})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- (2*p1-2*val1)+((2*p1^(2)+2*p1)/(length(y)-p1-1))
  return(val2)
}

AICc.l(par.est.l,form,df) 

#BIC
BIC.l <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  beta <- par
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.l(eta[j,1])})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- 2*log(length(y))-2*val1
  return(val2)
}

BIC.l(par.est.l,form,df) 




