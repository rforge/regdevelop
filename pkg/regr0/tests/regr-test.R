library(regr0)## <- do *NOT* change 'lib' here! {it must work on R-forge, CRAN, ..}
## source('../R/regr.R')
## library(regr0, lib="/u/stahel/R/regdevelop/pkg/regr0_1.0-5.tar.gz")
## attach("../misc/regr0.rda")
## attach("../misc/data-regr0.rda")
options(digits=3)

options(verbose=1)
## ========================================================================
## --- lm
t.d <- data.frame(y = c(1,5,4,4,3,8,6,7),
                   f1 = factor(c('a', 'b', 'b', 'c', 'a', 'a', 'c','a')),
                   x1 = c(1,2,3,4,5,6,7,8), x2 = c(4,2,8,7,5,3,4,3))
r.r1 <- regr(formula = y~f1+x1+x2, data=t.d)
r.r2 <- regr(formula = y~f1+x1+x2+f1:x1, data=t.d)
r.r2o <- regr(formula = y~ordered(f1)+x1+x2+f1:x1, data=t.d)

t.d[1,"x2"] <- NaN
r.r3 <- regr(y~f1+x1+x2+f1:x2, data=t.d)

r.r4 <- regr(y~1, data=t.d)

## DB()
names(t.d)
r.r5 <- regr(formula = y~f1+x1+x2, data=t.d, weights=0:7)

t.d$f2 <- ordered(t.d$f1)

t.ctr <- contr.wsumpoly(t.d)

r.r6 <- regr(formula = y~f2+x1+x2, data=t.d)
##- plot(r.r2,plotselect=
##-      c( qq = NA, yfit=2, ta=3, tascale = NA, weights = NA, hat = 0),ask=c.ask)

data(d.blast)
r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)

with(r.blast,
     stopifnot(all.equal(
     termtable[-1,"signif"],
                          c(3.40649491, -12.0192639, 8.1947910) )
##                          c( 3.477157958, -12.133975476,   8.169091379) )
               )
     )
stopifnot(length(residuals(r.blast))==nrow(d.blast))

plot(r.blast, smooth.group=location)
t.r <- regr(tremor ~ distance + charge + location, data=d.blast)
plot(r.blast, smresid=TRUE, xplot=F, plsel=c(ta=0, leverage=0))

dd <- cbind( d.blast, indicator= factor(1==1:nrow(d.blast)) )
t.rr <- regr(tremor ~ distance + charge + location + indicator, data=dd)
plot(t.rr, smresid=TRUE, xplot=F, plsel=c(ta=0, leverage=2))

r.ad <- add1(r.blast)
r.bl2 <- update(r.blast, ~.+location:log10(distance)+I(log10(charge)^2) )
r.bl3 <- step(r.bl2, trace=FALSE, k=4)
r.bl4 <- step(r.blast, expand=TRUE)

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

data(d.hail)
r.hail <- regr(logst(EGR) ~ logst(E0) * factor(SEED) + TB + TB:logst(E0) + H0,
               data=d.hail)
t.pr <- predict(r.hail, newdata=d.hail[1:5,])
stopifnot( all(abs(t.pr-fitted(r.hail)[1:5])<0.0001)  )

## splines
require(splines)
data(d.pollZH16d, package="regr0")
rr <- regr(log10(NO2) ~ bs(temp, df=5) +daytype, data=d.pollZH16d[-c(2,3),])
add1(rr)
## ------------------------------------
## all subsets
data(d.birthrates)
t.formula <- fertility ~ catholic + single24 + single49 + eAgric + eIndustry +
eCommerce + eTransport + eAdmin + gradeHigh + gradeLow + educHigh +
bornLocal + bornForeign + altitude + language

r.ae <- regrAllEqns(t.formula, data=d.birthrates, nbest=100, really.big=TRUE)
print(r.ae)
plot(r.ae, mar=9, ncharhorizontal=0, main="birthrates example",
     legend="topright")
## extract the formula of the best model
( t.formula <- regrAllEqnsXtr(r.ae) )
## ... and fit it
regr(t.formula, data=d.birthrates)

r.ae0 <-
  regrAllEqns(fertility ~ catholic + bornLocal + altitude + language +
                I(single49>30) -1,
              data=d.birthrates, force.in=~language, force.out=~altitude,
              codes=c(catholic="c",bornLocal="b",altitude="a",
                      "I(single49 > 30)"="s",language="L",foo="f"))

print(r.ae0, printcrit=T)
plot(r.ae0, legend=c(5.5,3000), cex=1)
                     
## ========================================================================
## robust [MM: now works thanks to quote(robustbase::lmrob) hack]
data(d.blast)
r.rob <-
  regr(log10(tremor) ~ location+log10(distance)+log10(charge),
       data = d.blast, robust=TRUE)
r.rob
modelTable( list(ls=r.blast, robust=r.rob) )
## ========================================================================
## ordered regression (using polr)
## load('../data/d.surveyenvir.rda')
data(d.surveyenvir)

## require(MASS)  ## should not be needed
r.survey <- regr(disturbance ~ age+education+location, data=d.surveyenvir,
                 contrasts="contr.treatment")

##- with( r.survey,
##-      stopifnot(
##-           all.equal(unname(coefficients[1:2]),
##-                           c(-0.0034971305,  0.0697008869), tol=100e-6)
##-          ,
##-           all.equal(unname(zeta),
##-                           c(-0.15772572, 1.37302143, 2.82205070 ), tol=100e-4)
##-          ) )

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
t.rs <- step(t.r, trace=FALSE)  ## ???
t.r <- r.babysurv <- regr(Survival~Weight+Age+Apgar1,data=t.d,family=binomial)
plot(r.babysurv,xplot=~Weight,cex=0.7,symbol.size=NULL,res.lim=c(-5,5))
plot(r.babysurv,glm.restype="cond", ask=FALSE)

mframe(2,2)
plresx(r.babysurv,vars=~Age+Apgar1+Apgar5+pH,data=d.babysurv,
       weights=FALSE,cex.lab=0.2)

plresx(t.r,data=t.d,vars=~.+Apgar5,sequence=TRUE, addcomp=TRUE)

data(d.babysurvGr)
t.d <- d.babysurvGr
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
## require(nnet)
##- r.mnom <- multinom(responsibility~age+education+disturbance+sex,
##-           data=d.surveyenvir)
r.mnom <- regr(responsibility~age+education+disturbance+sex,
          data=d.surveyenvir, family="multinomial", na.action=na.omit)
r.mnom2 <- regr(responsibility~age+education, data=d.surveyenvir[1:100,])
drop1(r.mnom2)
## !!! wieso geht drop1 in der regr-fn nicht?
## ===================================================================
## nonlinear
## sysfile('external/d.reacten.rda', package='regr0')
data(d.reacten)
t.d <- d.reacten[300:1700,]
r.lin <- regr(log10(q)~time,data=t.d)
t.cl <- unname(coefficients(r.lin))
t.r <- regr(q~theta1*exp(-theta2*time), data=t.d, nonlinear=TRUE,
             start=c(theta1=10^t.cl[1],theta2=t.cl[2]))

t.d <- subset(DNase, Run == 1)
t.r <- regr(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
            nonlinear = TRUE,
            data = t.d, control=list(maxiter=3, warnOnly=TRUE),
            start = list(Asym = 1, xmid = 0, scal = 1))
##- data(d.treated)
##- t.r <- regr(~weighted.MM(rate, conc, Vm, K), data = d.treated, nonlinear=T,
##-        start = list(Vm = 200, K = 0.1))  # d.treated is found in example(nls)
##- plot(t.r)
## ===================================================================
data("ovarian", package="survival")
require(survival) ## for dataset ovarian and >3 functions
##- t.rs <- survreg(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian,
##-     dist = "weibull")
##- t.rss <- summary(t.rs)
##- t.rs <- survreg(formula = Surv(log(futime), fustat) ~ ecog.ps + rx,
##-     data = ovarian, dist = "extreme") ## not the same
r.surv <- regr(formula = Surv(futime, fustat) ~ ecog.ps + rx,
               data = ovarian, family="weibull")
plot(r.surv)
r.coxph <- regr(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian)
length(residuals(r.coxph, type="martingale"))

data(bladder, package="survival")
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
data("tobin", package="survival")
dd <- tobin
dd[1,"durable"] <- NA
dd[2,"age"] <- NA
t.r <- regr(Tobit(durable) ~ age + quant, data=dd)
plot(t.r)

t.rf <- fitted(t.r)
t.rp <- predict(t.r)
t.rr <- residuals(t.r)
t.ci <- confint(t.r)
stopifnot(round(t.ci["quant",1],3)==-0.223)

add1(t.r)
##- set.seed(0)
##- dd$random <- rnorm(nrow(dd))
##- t.r2 <- regr(Tobit(durable) ~ age + quant + random, data=dd)
t.rst <- step(t.r)
stopifnot(round(t.rst$anova[3,"Deviance"],3)==0.835)
## =======================================================================
## trouble with complicated formula
t.loc <- d.blast$location
t.ri <- regr(formula = log10(tremor) ~ factor(t.loc)*log10(distance)+
               ordered(charge),
           data = d.blast, contrasts="contr.treatment")
r.fc <- fitcomp(t.ri)
plresx(t.ri)

## -----------------------------------------------------
## last
df <- data.frame(X=c(2,5,3,8), F=LETTERS[1:4], G=c(TRUE,FALSE,FALSE,TRUE))
last(df,3,-2)
last(df, ncol=1)
t.x <- matrix(11:30,5)
last(t.x,1)
last(t.x,,1)
last(t.x,1,drop=FALSE)
## ----------------------------------------------------------
plmatrix(cbind(log10(Petal.Width),Petal.Length)~Species, data=iris,
         col=as.numeric(Species)+1, xaxmar=1)
plmatrix(log10(tremor)~as.factor(device), data=d.blast, col=location)

plmboxes(tremor~location, data=d.blast, ilim=c(0,20),ilimext=0.05)
plmboxes(tremor~location, data=d.blast, at=c(1,3,4,7,NA,10,11,12),
         width=c(1,0.5))
