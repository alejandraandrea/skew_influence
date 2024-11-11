rm(list = ls())
library(optimx)
source('functions.R')
# Girls ------------------------------------------------------------------------
df<-read.table("./app/girls/girls.txt", sep="\t", header=TRUE)
X <- cbind(1,df$Ida)
y <- df$Mens
Ni <- df$Entr
pl <- ncol(X)
psl <- pl+1
EPSL<-matrix(NA, psl, 1)
form <- cbind(df$Mens,df$Entr-df$Mens)~df$Ida

## Logit  ----------------------------------------------------------------------
link <- 'logit'
fit.l <- glm(form, family = binomial(link = link), data = df)
EPL <- fit.l$coefficients
beta0 <- EPL
#lambda0 <- fbeta(par = beta0, y, X, Ni, link = link)
#theta0 <- c(lambda0, beta0)
theta0 <- c(0, beta0)
fit.sl <- optimx(par = theta0,  y = y, X = X, Ni = Ni, link = link,
                 fn = negll.sl.group,
                 method = c('BFGS'),
                 control = list(maxit = 500))

lambda <- coef(fit.sl)[1]
beta <- coef(fit.sl)[2:psl]


hessian.sl<- attributes(fit.sl)$details["BFGS", "nhatend"][[1]]
oim.sl<- hessian.sl
std.error<-sqrt(diag(solve(oim.sl)))
std.error

## Probit  ---------------------------------------------------------------------
link <- 'probit'
fit.l <- glm(form, family = binomial(link = link), data = df)
EPL <- fit.l$coefficients
beta0 <- EPL
#lambda0 <- fbeta(par = beta0, y, X, Ni, link = link)
#theta0 <- c(lambda0, beta0)
theta0 <- c(0, beta0)
fit.sl <- optimx(par = theta0,  y = y, X = X, Ni = Ni, link = link,
                 fn = negll.sl.group,
                 method = c('BFGS'),
                 control = list(maxit = 500))

lambda <- coef(fit.sl)[1]
beta <- coef(fit.sl)[2:psl]

hessian.sl<- attributes(fit.sl)$details["BFGS", "nhatend"][[1]]
oim.sl<- hessian.sl
std.error<-sqrt(diag(solve(oim.sl)))
std.error

## Skew-t  ---------------------------------------------------------------------
link <- 'skewt'
fit.l <- glm(form, family = binomial(link = 'logit'), data = df)
EPL <- fit.l$coefficients/1.7
beta0 <- EPL
#lambda0 <- fbeta(par = beta0, y, X, Ni, link = link)
#theta0 <- c(lambda0, beta0)
theta0 <- c(0, beta0)
fit.sl <- optimx(par = theta0,  y = y, X = X, Ni = Ni, link = link,
                 fn = negll.sl.group,
                 method = c('BFGS'),
                 control = list(maxit = 500))

lambda <- coef(fit.sl)[1]
beta <- coef(fit.sl)[2:psl]

hessian.sl<- attributes(fit.sl)$details["BFGS", "nhatend"][[1]]
oim.sl<- hessian.sl
std.error<-sqrt(diag(solve(oim.sl)))
std.error


