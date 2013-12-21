library(regr0)
# source('../R/regr.R')
options(digits=3)
# load("/u/stahel/data/regdata")
##- load("../data/d.blast.rda")
##- load("../data/d.surveyenvir.rda")
##- load("../data/d.rehab.rda")
##- load("../data/d.fossiles.rda")

options(verbose=1)
## ========================================================================
## --- lm
t.d <- data.frame(y = c(1,5,4,4,3,8,6,7),
                   f1 = factor(c('a', 'b', 'b', 'c', 'a', 'a', 'c','a')),
                   x1 = c(1,2,3,4,5,6,7,8), x2 = c(4,2,8,7,5,3,4,3))
r.r1 <- regr(formula = y~f1+x1+x2, data=t.d)
r.r2 <- regr(formula = y~f1+x1+x2+f1:x1, data=t.d)

t.d[1,"x2"] <- NaN
r.r3 <- regr(y~f1+x1+x2+f1:x2, data=t.d)

r.r4 <- regr(y~1, data=t.d)

r.r5 <- regr(formula = y~f1+x1+x2, data=t.d, weights=0:7)

##- plot(r.r2,plotselect=
##-      c( qq = NA, yfit=2, ta=3, tascale = NA, weights = NA, hat = 0),ask=c.ask)

data(d.blast)
r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)

with(r.blast,
     stopifnot(all.equal(
     testcoef[,"signif"],
                          c(13.58841618, 3.40649491, -12.0192639, 8.1947910))
                )
     )

r.ad <- add1(r.blast)
r.bl2 <- update(r.blast, ~.+location:log10(distance)+I(log10(charge)^2) )
r.bl3 <- step(r.bl2, trace=FALSE)

stopifnot(r.bl3$anova[2,"Step"]=="- I(log10(charge)^2)")

r.mt <- modelTable(list(large=r.bl2,reduced=r.bl3,original=r.blast))
r.mt2 <- r.mt[,-2]
stopifnot(c(round(r.mt2[4,1]$p,7))==0.2757776)

r.ct <- compareTerms(large=r.bl2,reduced=r.bl3,original=r.blast)

plresx(r.blast, vars=~distance,
       pch=d.blast$location, smooth.group=d.blast$location,
       smooth.col=c("blue","red","darkgreen","purple","brown",
         "orange","cyan","black"),
       smooth.legend="bottomright")
## unreproduzierbarer Fehler
## ========================================================================
## robust
data(d.blast)
r.rob <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast,
            robust=T)
## ========================================================================
## ordered regression (using polr)
## load('../data/d.surveyenvir.rda')
data(d.surveyenvir)

r.survey <- regr(disturbance~age+education+location, data=d.surveyenvir)

with( r.survey, 
     stopifnot(
          all.equal(unname(coefficients[1:2]),
                          c(-0.0034971305,  0.0697008869), tol=10e-6)
,
          all.equal(unname(zeta),
                          c(-0.15772572, 1.37302143, 2.82205070 ), tol=10e-4)
               )
     )

## ========================================================================
## multivariate regression
data(d.fossiles)
r.mregr <-
  regr(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+lChlorophyll+region+N,
                data=d.fossiles)
plot(r.mregr)
## ========================================================================
## Baby Survival
data(d.babysurv)
t.d <- d.babysurv
t.r <- regr(Survival~.,data=t.d)
t.rglm <- glm(Survival~.,data=t.d,family=binomial)
t.rs <- step(t.r, trace=F)  ## ???
t.r <- r.babysurv <- regr(Survival~Weight+Age+Apgar1,data=t.d,family=binomial)
plot(r.babysurv,xplot=~Weight,cex=0.7,symbol.size=NULL,res.lim=c(-5,5))
plot(r.babysurv,glm.restype="cond", ask=FALSE)

mframe(2,2)
plresx(r.babysurv,vars=~Age+Apgar1+Apgar5+pH,data=d.babysurv,
       weights=F,cex.lab=0.2)

plresx(t.r,data=t.d,vars=~.+Apgar5,sequence=T)

data(d.babysurv.w)
t.d <- d.babysurv.w
t.r <- regr(cbind(Survival.1,Survival.0)~Weight,data=t.d,family=binomial)
t.r <- glm(cbind(Survival.1,Survival.0)~Weight,data=t.d,family=binomial)
plot(t.r,cex=1.2)

## ===========================================================
## Gamma
data(d.transaction)
r.gam <- regr(Time~Type1+Type2,data=d.transaction,family=Gamma)
r.gami <- regr(Time~Type1+Type2,data=d.transaction,family=Gamma(link=identity)) 

## ===========================================================
##- multinom
data(d.surveyenvir)
require(nnet)
r.mnom <- multinom(responsibility~age+education+disturbance+sex,
          data=d.surveyenvir) 
r.mnom <- regr(responsibility~age+education+disturbance+sex,
          data=d.surveyenvir, family="multinomial", na.action=na.omit) 
r.mnom2 <- regr(responsibility~age+education, data=d.surveyenvir[1:100,])
drop1(r.mnom2) 

## ===================================================================
## nonlinear
## sysfile('external/d.reacten.rda', package='regr0')
data(d.reacten)
t.d <- d.reacten[300:1700,]
r.lin <- regr(log10(q)~time,data=t.d)
t.cl <- unname(coefficients(r.lin))
t.r <- regr(q~theta1*exp(-theta2*time), data=t.d, nonlinear=T,
             start=c(theta1=10^t.cl[1],theta2=t.cl[2]))

t.d <- subset(DNase, Run == 1)
t.r <- regr(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
            nonlinear = TRUE, 
            data = t.d, control=list(maxiter=3, warnOnly=T),
            start = list(Asym = 1, xmid = 0, scal = 1))
##- data(d.treated)
##- t.r <- regr(~weighted.MM(rate, conc, Vm, K), data = d.treated, nonlinear=T,
##-        start = list(Vm = 200, K = 0.1))  # d.treated is found in example(nls)
##- plot(t.r)
## ===================================================================
require(survival)
t.rs <- survreg(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, 
    dist = "weibull")
t.rss <- summary(t.rs)
##- t.rs <- survreg(formula = Surv(log(futime), fustat) ~ ecog.ps + rx,  
##-     data = ovarian, dist = "extreme") ## not the same
r.surv <- regr(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian,
            family="weibull")
plot(r.surv)
r.coxph <- regr(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian)

bladder1 <- bladder[bladder$enum < 5, ]
t.cph <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + 
               cluster(id), bladder1) 
t.r <- regr(Surv(stop, event) ~ rx + size + number, bladder1) 
r.sr <- regr(Surv(stop, event) ~ rx + size + number, bladder1, family="weibull") 

bladder1$sizej <- jitter(bladder1$size)
plresx(r.sr, vars="sizej")
## d.cmbscores <- t.d
## ===================================================================
## Tobit

##- t.d <- d.cmbscores
##- dim(na.omit(t.d))
##- t.d$y <- Tobit(t.d$EVAPOR, 10, log=T)
##- t.r <- regr(y~Temp+Time+lWindspeed, data=t.d)

require(survival)
