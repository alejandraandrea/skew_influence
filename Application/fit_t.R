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
#t model

#Function to calculate the CDF t-student
nu<-5 #Suggest by Gómez et al. (2007) 
pdf.t<-function(z){dt(z,nu)}
cdf.t<-function(z){pt(z,nu)}

#Negative log-likelihood function
negll.t <- function(par,formula,data){
  #formula <- cbind(data$Mens,data$Entr-data$Mens)~data$Ida 
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

theta.ini<-c(-21,1.6) #glm logistic 
pt<-length(theta.ini)
EPT<-matrix(NA,pt,1)
fit.t<-optimx(par=theta.ini,formula=form,data=df,fn=negll.t,method=c("BFGS","L-BFGS-B"))
coef.opt<-coef(fit.t)
EPT[1,1]<-coef.opt[1,1]#beta0 BFGS
EPT[2,1]<-coef.opt[1,2]#beta1 BFGS
EPT

#P-value 
coef.opt<-coef(fit.t)
par.est.t<-coef.opt[1,]
p1<-length(par.est.t)
hessian.t<- attributes(fit.t)$details["BFGS", "nhatend"][[1]]
oim.t<- hessian.t
oim.t
std.error<-sqrt(diag(solve(oim.t)))
std.error
zvalue<-par.est.t/std.error
zvalue
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
  IC095[j,]<-par.est.t[j]+c(-1,1)*qnorm(0.975)*std.error[j]
}
IC095

#Fitted probabilities
fit.prob.t<-function(par,formula,data){
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
  return(prob)
}

fit.prob<-fit.prob.t(par.est.t,form,df)
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
EPT.envel<-matrix(NA,pt,1)
for(i in 1:m){
  for(j in 1:N){
    dif<-runif(Ni[j])-fit.prob[j];
    dif[dif>=0]<-0;
    dif[dif<0]<-1;
    yi[j]<-sum(dif);
  }
  theta.ini<-c(-21,1.6) #glm logistic 
  dfi<-data.frame(yi,Ni,xi)
  formi<-cbind(dfi$yi,dfi$Ni-dfi$yi)~dfi$xi
  fiti<-optimx(par=theta.ini,formula=formi,data=dfi,fn=negll.t,method="BFGS")  
  par.est.i<-coef(fiti)[1,]
  fit.prob.i<-fit.prob.t(par.est.i,formi,dfi)
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
#options(warn = 1) #Error: (converted from warning) glm.fit: fitted probabilities numerically 0 or 1 occurred
fac.sat<-factor(1:length(df$Mens)) 
df.sat<-df
df.sat$fac.sat<-fac.sat
form.sat<- cbind(df$Mens,df$Entr-df$Mens)~df.sat$fac.sat
fit.l.sat<- glm(formula=form.sat,family=binomial(link="logit"),data=df.sat,maxit=1000)
betas.ini.sat<-fit.l.sat$coefficients

#Optimization
theta.ini.sat<-betas.ini.sat
fit.t.sat<-optimx(par=theta.ini.sat,formula=form.sat,data=df.sat,fn=negll.t,method = "BFGS")
par.est.t.sat<-coef(fit.t.sat)[1,]
loglike.t.sat<- -negll.t(par.est.t.sat,form.sat,df.sat) 
loglike.t.sat

#Null model
form.null<- cbind(df$Mens,df$Entr-df$Mens)~1
fit.t.null<- glm(formula=form.null,family=binomial(link="logit"),data=df)
betas.ini.null<-fit.t.null$coefficients

#Optimization
theta.ini.null<-betas.ini.null
fit.t.null<-optimx(par=theta.ini.null,formula=form.null,data=df,fn=negll.t,method =c("BFGS","L-BFGS-B"))
fit.t.null 
coef.opt.null<-coef(fit.t.null)
par.est.t.null<-coef.opt.null[1,]
loglike.t.null<- -negll.t(par.est.t.null,form.null,df) 
loglike.t.null

#Proposed model
loglike.prop<- -negll.t(par.est.t,form,df) 
loglike.prop

#Residual deviance
res.dev<-2*(loglike.t.sat - loglike.prop) 
res.dev

#Null deviance
res.null<-2*(loglike.t.sat - loglike.t.null) 
res.null

#AIC
AIC.t <- function(par,formula,data){
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
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- 2*p1-2*val1
  return(val2)
}

AIC.t(par.est.t,form,df) 

AICc.t <- function(par,formula,data){
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
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- (2*p1-2*val1)+((2*p1^(2)+2*p1)/(length(y)-p1-1))
  return(val2)
}

AICc.t(par.est.t,form,df) 

#BIC
BIC.t <- function(par,formula,data){
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
  val1<- sum(log(choose(size,y))+y*log(prob)+(size-y)*log(1-prob))
  val2<- 2*log(length(y))-2*val1
  return(val2)
}

BIC.t(par.est.t,form,df) 





