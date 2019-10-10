## install.packages("regr", repos="http://R-forge.R-project.org")
## install.packages("~/R/regdevelop/pkg/regr_1.0-5.tar.gz", repos=NULL)
## library(regr0)## <- do *NOT* change 'lib' here! {it must work on R-forge, CRAN, ..}
## library(regr) ##, lib="~/R/regdevelop/pkg/regr.Rcheck")
## attach("../misc/regr0.rda")
## attach("../misc/data-regr0.rda")
options(digits=3)
require(regr)
## require(plgraphics)
## options(verbose=1)
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
r.r5 <- regr(formula = y~f1+x1+x2, data=t.d, weights=0:7)

t.d$f2 <- ordered(t.d$f1)

t.ctr <- contr.wsumpoly(t.d)

r.r6 <- regr(formula = y~f2+x1+x2, data=t.d)

data(d.blast, package="plgraphics")
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

## plot(r.blast, smooth.group=location)
## t.r <- regr(tremor ~ distance + charge + location, data=d.blast)
plot(r.blast, smresid=TRUE, xvars=F, plotselect=c(resfit=0, leverage=0))

dd <- cbind( d.blast, indicator= factor(1==1:nrow(d.blast)) )
t.rr <- regr(tremor ~ distance + charge + location + indicator, data=dd)
plot(t.rr, smresid=TRUE, xvars=F, plotselect=c(resfit=0, leverage=2))

r.ad <- add1(r.blast)
r.bl2 <- update(r.blast, ~.+location:log10(distance)+I(log10(charge)^2) )
r.bl3 <- step.regr(r.bl2, trace=FALSE, k=4)
r.bl4 <- step.regr(r.blast, expand=TRUE)

stopifnot(r.bl3$anova[2,"Step"]=="- I(log10(charge)^2)")

r.mt <- modelTable(list(large=r.bl2,reduced=r.bl3,original=r.blast))
r.mt2 <- r.mt[,-2]
stopifnot(c(round(r.mt2[4,1]$p,7))==0.2757776)

r.ct <- compareTerms(large=r.bl2,reduced=r.bl3,original=r.blast)

r.sim <- simresiduals(r.blast)

k <- d.blast
t.r <- regr(tremor~distance+charge+location, data=k)
add1(t.r) ## has produced an error due to wrong eval() environment 

plresx(r.blast, xvars=~distance,
       plab=d.blast$location, smooth.group=d.blast$location,
       colors.smgrp=c("blue","red","darkgreen","purple","brown",
         "orange","cyan","black"),
       smooth.legend="bottomright")

## splines
require(splines)
data(d.pollZH16d)
rr <- regr(log10(NO2) ~ bs(temp, df=5) +daytype, data=d.pollZH16d[-c(2,3),])
add1(rr)
## ------------------------------------
## all subsets
data(d.birthrates, package="plgraphics")
t.formula <- fertility ~ catholic + single24 + single49 + eAgric + eIndustry +
  eCommerce + eTransport + eAdmin + gradeHigh + gradeLow + educHigh +
  bornLocal + bornForeign + altitude + language

r.ae <- regrAllEqns(t.formula, data=d.birthrates, nbest=100, really.big=TRUE)
print(r.ae)
plot(r.ae, mar=9, ncharhorizontal=0, main="birthrates example",
     legend="topright")
## extract the formula of the best model
t.formula <- regrAllEqnsXtr(r.ae)
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
data(d.blast, package="plgraphics")
r.rob <-
  regr(log10(tremor) ~ location+log10(distance)+log10(charge),
       data = d.blast, robust=TRUE)
r.rob
modelTable( list(ls=r.blast, robust=r.rob) )
## ========================================================================
## ordered regression (using polr)
data(d.surveyenvir)

r.survey <- regr(disturbance ~ age+education+location, data=d.surveyenvir,
                 contrasts="contr.treatment")
plot(r.survey, plotselect=c(resfit=3, default=0), xvars=F, mf=c(1,1))

##- with( r.survey,
##-      stopifnot(
##-           all.equal(unname(coefficients[1:2]),
##-                           c(-0.0034971305,  0.0697008869), tol=100e-6)
##-          ,
##-           all.equal(unname(zeta),
##-                           c(-0.15772572, 1.37302143, 2.82205070 ), tol=100e-4)
##-          ) 

## ========================================================================
## multivariate regression
data(d.fossileSamples, package="plgraphics")
r.mregr <-
  regr(cbind(sAngle,lLength,rWidth)~SST+Salinity+lChlorophyll+Region+N,
                data=d.fossileSamples)
plot(r.mregr)
## ========================================================================
## Baby Survival
data(d.babysurvival, package="plgraphics")
t.d <- d.babysurvival
t.d$Age[2] <- NA
t.r <- regr(Survival~.,data=t.d)
t.rglm <- glm(Survival~.,data=t.d,family=binomial)
t.rs <- step.regr(t.r, trace=FALSE)  ## ???

t.r <- r.babysurvival <- regr(Survival~Weight+Age+Apgar1,data=t.d,family=binomial)
plot(r.babysurvival, xvars=~Weight,cex.plab=0.7, ylim=c(-5,5))
plot(r.babysurvival, glm.restype="condquant")

"is plgraphics in the path?"
search()
ls("package:plgraphics")
help(package="plgraphics")

plmframes(2,2)
plresx(r.babysurvival, xvars=~Age+Apgar1+Apgar5+pH, data=d.babysurvival,
       weight=FALSE,cex.plab=0.2)

plresx(t.r, data=t.d, xvars=~.+Apgar5,sequence=TRUE, addcomp=TRUE)

data(d.babysurvGr, package="plgraphics")
t.d <- d.babysurvGr
t.r <- regr(cbind(Survival.1,Survival.0)~Weight,data=t.d,family=binomial)
t.r <- glm(cbind(Survival.1,Survival.0)~Weight,data=t.d,family=binomial)
plot(t.r,cex=1.2)

## ===========================================================
## Gamma
data(d.transaction)
r.gam <- regr(Time~Type1+Type2, data=d.transaction, family=Gamma)
r.gami <- regr(Time~Type1+Type2, data=d.transaction,
               family=Gamma(link=identity))

##- data(wafer, package="faraway")
##- r.gammalog <-
##-   regr(resist ~ x1 + x2 + x3 + x4, family = Gamma(link = "log"), data=wafer)
##- stopifnot(r.gammalog[["family"]][["link"]]=="log")

## ===========================================================
##- multinom
data(d.surveyenvir)
## require(nnet)
##- r.mnom <- multinom(responsibility~age+education+disturbance+sex,
##-           data=d.surveyenvir)
r.mnom <- regr(responsibility~age+education+disturbance+sex,
          data=d.surveyenvir, family="multinomial")
r.mnom2 <- regr(responsibility~age+education, data=d.surveyenvir[1:100,])
drop1(r.mnom2)
## !!! wieso geht drop1 in der regr-fn nicht?
## ===================================================================
## nonlinear

data(d.reacten)
t.d <- d.reacten[300:1700,]
r.lin <- regr(log10(q)~time, data=t.d)
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
## require(survival) ## for dataset ovarian and >3 functions
##- t.rs <- survival::survreg(formula = survival::Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian,
##-     dist = "weibull")
##- t.rss <- summary(t.rs)
##- t.rs <- survival::survreg(formula = survival::Surv(log(futime), fustat) ~ ecog.ps + rx,
##-     data = ovarian, dist = "extreme") ## not the same
r.surv <- regr(formula = survival::Surv(futime, fustat) ~ ecog.ps + rx,
               data = ovarian, family="weibull")
## plot(r.surv)  # !!!
r.coxph <- regr(formula = survival::Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian)
length(residuals(r.coxph, type="martingale"))

data(bladder, package="survival")
bladder1 <- bladder[bladder$enum < 5, ]
t.cph <- survival::coxph(survival::Surv(stop, event) ~ (rx + size + number) * strata(enum) +
               cluster(id), bladder1)
t.r <- regr(survival::Surv(stop, event) ~ rx + size + number, bladder1)
r.sr <- regr(survival::Surv(stop, event) ~ rx + size + number, bladder1, family="weibull")

bladder1$sizej <- jitter(bladder1$size)
plresx(r.sr, xvars="size") ## !!! falsche limits
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
## stopifnot(round(t.ci["quant",1],3)==-0.223)

add1(t.r)
##- set.seed(0)
##- dd$random <- rnorm(nrow(dd))
##- t.r2 <- regr(Tobit(durable) ~ age + quant + random, data=dd)
t.rst <- step.regr(t.r)
stopifnot(round(t.rst$anova[3,"Deviance"],3)==1.371) ## 0.835

t.rr <- survival::survreg(Tobit(durable) ~ quant, data=tobin, dist="gaussian")
add1(t.rr, scope=~.+age)
step.regr(t.rr, scope=~.+age)

## =======================================================================
## trouble with complicated formula
t.loc <- d.blast$location
t.ri <- regr(formula = log10(tremor) ~ factor(t.loc)*log10(distance)+
               ordered(charge),
           data = d.blast, contrasts="contr.treatment")
r.fc <- fitcomp(t.ri)
plresx(t.ri)

## -----------------------------------------------------
