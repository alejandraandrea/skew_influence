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
#Skew-logistic model
 
#Function to calculate the CDF skew-logistic
pdf.sl<-function(z,lambda){2*dlogis(z,log=FALSE)*pnorm(lambda*z)}
cdf.sl<-function(z,lambda){integrate(pdf.sl,lower=-Inf, upper=z,
                                     lambda =lambda,stop.on.error=FALSE)$value}

#Negative log-likelihood function
negll.sl <- function(par,formula,data){
  #formula <- cbind(data$Mens,data$Entr-data$Mens)~data$Ida 
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.sl(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val<- - sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  return(val)
}

library(optimx)
fit.l<-glm(form,family=binomial(link="logit"),data=df)
EPL<-fit.l$coefficients
theta.ini<-c(0,EPL)
psl<-length(theta.ini)
EPSL<-matrix(NA,psl,1)
fit.sl<-optimx(par=theta.ini,formula=form,data=df,fn=negll.sl,method=c("BFGS","L-BFGS-B"))
coef.opt<-coef(fit.sl)
EPSL[1,1]<-coef.opt[1,1]#lambda BFGS
EPSL[2,1]<-coef.opt[1,2]#beta0 BFGS
EPSL[3,1]<-coef.opt[1,3]#beta1 BFGS
EPSL

#P-value 
coef.opt<-coef(fit.sl)
par.est.sl<-coef.opt[1,]
p1 <- length(par.est.sl)
hessian.sl<- attributes(fit.sl)$details["BFGS", "nhatend"][[1]]
oim.sl<- hessian.sl 
oim.sl
std.error<-sqrt(diag(solve(oim.sl)))
std.error
zvalue<-par.est.sl/std.error
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

#IC
IC095<-matrix(NA,3,2)
for(j in 1:3){
  IC095[j,]<-par.est.sl[j]+c(-1,1)*qnorm(0.975)*std.error[j]
}
IC095

#Fitted Probabilities 
fit.prob.sl<-function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.sl(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  return(prob)
}

fit.prob<-fit.prob.sl(par.est.sl,form,df)
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
EPSL.envel<-matrix(NA,psl,1)
for(i in 1:m){
  for(j in 1:N){
    dif<-runif(Ni[j])-fit.prob[j];
    dif[dif>=0]<-0;
    dif[dif<0]<-1;
    yi[j]<-sum(dif);
  }
  theta.ini<-c(0,EPL)  
  dfi<-data.frame(yi,Ni,xi)
  formi<-cbind(dfi$yi,dfi$Ni-dfi$yi)~dfi$xi 
  fiti<-optimx(par=theta.ini,formula=formi,data=dfi,fn=negll.sl,method="BFGS")  
  par.est.i<-coef(fiti)[1,]
  fit.prob.i<-fit.prob.sl(par.est.i,formi,dfi)
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
df<-df2
fac.sat<-factor(1:length(df$Mens)) 
df.sat<-df
df.sat$fac.sat<-fac.sat
form.sat<- cbind(df$Mens,df$Entr-df$Mens)~df.sat$fac.sat
fit.l.sat<- glm(formula=form.sat,family=binomial(link="logit"),data=df.sat,maxit=1000)
betas.ini.sat<-fit.l.sat$coefficients

#Optimization
theta.ini.sat<-c(0,betas.ini.sat)
fit.sl.sat<-optimx(par=theta.ini.sat,formula=form.sat,data=df.sat,fn=negll.sl,method = "BFGS")
par.est.sl.sat<-coef(fit.sl.sat)
loglike.sl.sat<- -negll.sl(par.est.sl.sat,form.sat,df.sat) 
loglike.sl.sat

#Null model
form.null<- cbind(df$Mens,df$Entr-df$Mens)~1
fit.l.null<- glm(formula=form.null,family=binomial(link="logit"),data=df)
betas.ini.null<-fit.l.null$coefficients

#Optimization
theta.ini.null<-c(0,betas.ini.null)
fit.sl.null<-optimx(par=theta.ini.null,formula=form.null,data=df,fn=negll.sl,method =c("BFGS","L-BFGS-B"))
fit.sl.null  
coef.opt.null<-coef(fit.sl.null)
par.est.sl.null<-coef.opt.null[1,]
loglike.sl.null<- -negll.sl(par.est.sl.null,form.null,df) 
loglike.sl.null

#Proposed model
form<-form2
df<-df2
loglike.prop<- -negll.sl(par.est.sl,form,df) 
loglike.prop

#Residual deviance
res.dev<-2*(loglike.sl.sat - loglike.prop) 
res.dev

#Null deviance
res.null<-2*(loglike.sl.sat - loglike.sl.null) 
res.null

#AIC
AIC.sl <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.sl(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- 2*p1-2*val1
  return(val2)
}

AIC.sl(par.est.sl,form,df) 

AICc.sl <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.sl(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- (2*p1-2*val1)+((2*p1^(2)+2*p1)/(length(y)-p1-1))
  return(val2)
}

AICc.sl(par.est.sl,form,df) 

#BIC
BIC.sl <- function(par,formula,data){
  mf <- model.frame(formula,data)
  mr <- model.extract(mf,"response")
  size <- mr[,2]+mr[,1]
  y <- mr[,1]
  X  <- model.matrix(formula,data=mf)
  p1 <- length(par)
  lambda <- par[1]
  beta <- par[2:p1]
  eta <- X%*%beta
  prob <- sapply(1:length(y),function(j){cdf.sl(eta[j,1],lambda)})
  prob<- ifelse(prob>=0.9999,0.9999, ifelse(prob<=0.0001,0.0001,prob))
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- 2*log(length(y))-2*val1
  return(val2)
}

BIC.sl(par.est.sl,form,df) 



