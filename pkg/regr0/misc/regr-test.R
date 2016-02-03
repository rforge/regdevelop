## attach("/u/stahel/regr0/misc/regr0.rda")
require(regr0,lib="/u/stahel/R/regdevelop/pkg/regr0.Rcheck")
##- save(list=ls(), file="data-regr0.rda")
attach("data-regr0.rda")
## ---------------------------------------------------------------
## --- regr.Rd examples
# data(d.blast)
options(verbose=1)

( r.blast <-
  regr(logst(tremor)~location+log10(distance)+log10(charge), data=d.blast) )
plot(r.blast)

plot(r.blast, plotselect=c(ta=2),xplot=F, seq=F, mf=1, cex.lab=1)

userOptions(project='regr0.demo', step='blast')
# u.pslatex('p-plotregr-ta')
##-   options(colors.ra =
##-           c("black","gray4","blue4","blue3","darkgreen","green",
##-             "burlywood4","burlywood3","burlywood4"))
par(lwd=2)
plot(r.blast, plotselect=c(leverage=2), xplot=F, pch=1, mf=1,
     lwd=c(2,1.5,2,1.5,1.5,1,1,1,1))
# ps.end()

t.pr <- predict(r.blast, interval="prediction")

t.r <- regr(tremor~distance+charge, data=d.blast)
t.xd <- xdistResdiff(t.r, nsim=1000)
plot(t.xd)

dd <- d.blast
dd$random <- c(NA,NA,rnorm(nrow(dd)-2))
t.r <- regr(log10(tremor)~log10(distance)+log10(charge)+random, data=dd)
r.st <- step(t.r)

## environment problem
f.rr <- function(fo) {
  ldt <- d.blast[1:20,]
  regr(fo, data=ldt)
}
f.rr(tremor~distance)

t.rlm <- lm(formula = log10(tremor) ~ location + log10(distance) +
                log10(charge) + I(log10(charge)^2) + location:log10(distance),
            data = d.blast)
t.r <- lm(formula = log10(tremor) ~ log10(charge) + I(log10(charge)^2)
          +location, data = d.blast)
t.rr <- lm(formula = log10(tremor) ~ location + log10(distance) +
                log10(charge) + location:log10(distance),
           data = d.blast) ## works

## ------------------------------------------------------------
t.x <- factor(c(1,1,1,2,3,3))
contr.wsum(t.x)
contr.wsum(data.frame(x=t.x))
t.x2 <- factor(c("a","b","c","a","c","b"))
t.d <- data.frame(fac=t.x, f2=t.x2, y=c(1,2,3,10,6,8))
t.r <- regr(y~fac, t.d, contrasts="contr.wsum")
t.r
sum(table(t.x)*t.r$allcoef$fac[,"estimate"])
model.matrix(y~fac, t.d)
model.matrix(y~fac, t.d, contrasts.arg=list(fac="contr.wsum"))

t.x <- factor(c(1,1,1,1,1,2,2,2,3,3,3,3))
t.x2 <- factor(c("a","b","c","a","a","a","b","c","a","b","c","a"))
t.d <- data.frame(fac=t.x, f2=t.x2, y=c(1:5,10:12,6:3))
model.matrix(y~fac*f2, t.d, contrasts.arg=list(fac="contr.wsum"))
t.r <- regr(y~fac*f2, t.d, contrasts=c("contr.wsum","contr.poly"))
t.ri <- regr(y~fac:f2, t.d, contrasts=c("contr.wsum","contr.poly"))
t.rlm <- lm(y~fac:f2, t.d, contrasts=c("contr.wsum","contr.poly"))
model.matrix(t.r)
predict(t.r, newdata=rbind(t.d,c(3,NA)))
fitted(t.r)
t.rlm <- lm(y~fac, t.d)
fitted(t.rlm)

t.d[,"fac"] <- "3"
t.r <- regr(y~as.numeric(fac)*f2, t.d, contrasts=c("contr.wsum","contr.poly"))


t.ri <- regr(formula = log10(tremor) ~ factor(location)*log10(distance)+
               factor(charge),
           data = d.blast)
t.ri <- regr(formula = log10(tremor) ~ factor(location)*log10(distance)+
               factor(charge),
           data = d.blast, contrasts="contr.treatment")

r.fc <- fitcomp(t.ri)
plresx(t.ri)

## ------------------------------------------------------------
## Anorexia
data(anorexia, package="MASS")
r.anorexia <- regr(Postwt ~ Treat + Prewt + offset(Prewt),
                   data = anorexia)
## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
     ## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
d.dob <- data.frame(group = gl(2,10,20, labels=c("Ctl","Trt")),
                    weight = c(ctl, trt))
(r.dob <- regr(weight ~ group, data=d.dob))

## multinomial regression
## data(d.surveyenvir)
dd <- d.surveyenvir
class(dd$disturbance) <- "factor"
t.rm <- multinom(disturbance~age+education+location, data=dd)
t.rm <- multinom(disturbance~age+education+location, data=na.omit(dd))
t.rm <- multinom(disturbance~age+education+location, data=dd, na.action=na.exclude)
t.rr <- regr(disturbance~age+education+location, data=dd)

## ordered regression
dd <- d.surveyenvir
t.rp <- polr(disturbance~age+education+location, data=dd)
t.rr <- regr(disturbance~age+education+location, data=dd)
##- t.r <- regr(Sat ~ Infl + Type + Cont, weights = housing$Freq, data = housing)
##- plot(t.r)

## multivariate regression
## data(d.fossiles)
r.mregr <-
  regr(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+lChlorophyll+region+N,
                data=d.fossiles)
plot(r.mregr)
## ------------------------------------------------
c.ask <- F
t.d <- data.frame(y = c(1,5,4,4,3,8,6,7),
                   f1 = factor(c('a', 'b', 'b', 'c', 'a', 'a', 'c','a')),
                   x1 = c(1,2,3,4,5,6,7,8), x2 = c(4,2,8,7,5,3,4,3))
t.r <- regr(formula = y~x1, data=t.d)
t.r <- regr(formula = y~f1+x1+x2, data=t.d)
t.r <- regr(formula = y~f1+x1+x2+f1:x1, data=t.d)
t.d[1,"x2"] <- NaN
t.r <- regr(y~f1+x1+x2, data=t.d)
plot(t.r,plotselect=
     c( qq = NA, yfit=2, ta=3, tascale = NA, weights = NA, leverage = 0),
     ask=c.ask)
t.r <- regr(y~f1+x1+x2+f1:x2, data=t.d)

data(anorexia, package="MASS")
r.anorex <- regr(Postwt ~ Prewt + Treat + offset(Prewt),
                 family = gaussian, data = anorexia)

r.savings <- regr(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
plot(r.savings, ask=c.ask)
## Spreng
t.d <- d.spreng
t.d <- t.d[t.d$stelle<=2,]
t.r <- regr(log10(ersch)~factor(stelle)*log10(dist)+log10(ladung),
            data=d.spreng)
t.form <- log10(ersch)~Stelle+log10(dist)+log10(ladung)
t.r <- regr(t.form, data=t.d)
t.rw <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=t.d,
             weights = 5+1:nrow(t.d))
plot(t.rw,ask=c.ask)
plresx(t.r,var=names(t.d),ask=c.ask)
plresx(t.r,var=,ask=c.ask, addcomp=T)

## adequate weights
t.y <- fitted(t.r)
t.e <- rnorm(length(t.y))/sqrt(4+1:length(t.y))
t.e <- t.e*t.r$sigma/mean(abs(t.e))
t.y <- t.y+t.e
#plot(t.y,log10(t.d$ersch))

t.rww <- regr(t.y~Stelle+log10(dist)+log10(ladung), data=t.d,
             weights = 5+1:nrow(t.d))
plot(t.rww,ask=c.ask)



t.d <- d.spreng
t.d$ersch[3:5] <- NA
t.d$ladung[5:7] <- NA
t.r <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=t.d,
             weights = 5+1:nrow(t.d))  
##- t.r <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=t.d,
##-              weights = 5+1:nrow(t.d), subset=sample(1:nrow(t.d), nrow(t.d)/2))
t.r <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=t.d,
             weights = 5+1:nrow(t.d), subset=stelle<=4)
## round(runif(nrow(t.d)))
plot(t.r,ask=c.ask)
drop1(t.r)
add1(t.r)
t.d$redist[15:20] <- NA
add1(t.r,scope=~redist)
t.r <- regr(log10(ersch)~1, data=t.d)
t.r1 <- step(t.r, scope=~.+Stelle+log10(dist)+log10(ladung),trace=F)
t.r2 <- update(t.r1, formula=.~.-log10(ladung))
t.r3 <- update(t.r2, formula=.~.-Stelle)
##- t.mt <- modelTable(c("t.r","t.r1"))
r.mt <- modelTable(list(large=t.r1,reduced=t.r2,small=t.r3))
r.mt[,-2]
compareTerms(large=t.r1,reduced=t.r2,small=t.r3)

t.r <- regr(Fertility ~ cut(Agriculture, breaks=4) + Infant.Mortality,
             data=swiss)
t.lm <- lm(Fertility ~ cut(Agriculture, breaks=4) + Infant.Mortality,
             data=swiss)
t.aov <- aov(Fertility ~ cut(Agriculture, breaks=4) + Infant.Mortality,
             data=swiss)
dummy.coef(t.r)
dummy.coef(t.lm)
dummy.coef.regr(t.lm)
t.r$allcoef
## -------------------------------------
## robust
t.r <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=d.spreng)
t.rr <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=d.spreng,
             robust=T)
modelTable(list(t.r,t.rr))
## -------------------------------------
plot(t.r)
## Baby Survival
##- t.d <- read.table("/u/stahel/data/babysurvival.dat",sep=",",header=T)
##- t.d <- t.d[t.d$Age<35,]
##- d.babysurv <- t.d
t.d <- d.babysurv
##- t.d$y <- factor(t.d$Survival)
##- t.r <- regr(y~.-Survival,data=t.d)
t.r <- regr(Survival~.,data=t.d) # family=binomial
t.rglm <- glm(Survival~.,data=t.d, family=binomial) 
t.rs <- step(t.r, trace=F)  ## ???
t.r <- r.babysurv <- regr(Survival~Weight+Age+Apgar1,data=t.d)
plot(r.babysurv,xplot=~Weight,cex=0.7,symbol.size=NULL,res.lim=c(-5,5))
plot(r.babysurv,glm.restype="cond", ask=c.ask)

mframe(2,2)
plresx(r.babysurv,vars=~Age+Apgar1+Apgar5+pH,data=d.babysurv,
       weights=F,cex.lab=0.2, ask=c.ask)

plresx(t.r,data=t.d,vars=~.+Apgar5,sequence=T, ask=c.ask)

t.d <- na.omit(t.d)
t.r <- regr(Survival~Weight+Age+Apgar1,data=t.d,family=binomial)
t.fit <- fitted(t.r)
t.rs <- step(t.r, k=5)
##- t.nsim <- 1000
##- t.sim <- rep(NA,t.nsim)
##- for (li in 1:t.nsim) {
##- t.d$Survival <- rbinom(nrow(t.d),1,t.fit)
##- t.r <- regr(Survival~Weight+Age+Apgar1,data=t.d,family=binomial)
##- t.sim[li] <- testBinFit(t.r, nclass=10)$p
##- }
##- hist(t.sim)
##- r.simtBFSurv10 <- t.sim
##- t.d <- read.table("/u/stahel/data/babysurvival-w.dat",sep=",",header=T)
##- t.d <- t.d[,-1]
##- d.babysurv.w <- t.d
t.d <- d.babysurv.w
t.r <- regr(cbind(Survival.1,Survival.0)~Weight,data=t.d,family=binomial)
t.r <- regr(cbind(Survival.1,Survival.0)~Weight,data=t.d,family=binomial,
            calcdisp=F)
#t.r <- glm(cbind(Survival.1,Survival.0)~Weight,data=t.d,family=binomial)
plot(t.r,cex=1.2, ask=c.ask)

## ===========================================================
# source("/u/stahel/data/umweltumf.dump")
t.d <- d.umweltumf
t.d$Beeintr2 <- as.numeric(t.d$Beeintr)>1
t.dg <- data.frame(table(t.d[,c("Ortsgroesse","Schule","Partei","Beeintr2")]))
t.nn <- nrow(t.dg)/2
t.db <- data.frame(t.dg[1:t.nn,c(1:3,5)],Beeintr.y=t.dg[t.nn+1:t.nn,"Freq"])
names(t.db)[4] <- "Beeintr.n"
t.r <- r.lgr <- regr(cbind(Beeintr.y,Beeintr.n)~.,data=t.db,family=binomial)
## !!! bug
plot(t.r, ask=c.ask)
t.r <- r.l <- regr(Beeintr2 ~ Alter+Geschlecht+Schule+Wohnlage+
                    Ortsgroesse+Partei,data=t.d,family=binomial)
## ===========================================================
##- t.d <- read.table('/u/stahel/data/transaction.dat',header=T)
##- d.transaction <- t.d
t.d <- d.transaction
t.rglm <- glm(Time~Type1+Type2,data=t.d,family=Gamma)
# summary(t.r)  ## !!! avoid deviance table
## !!! Gamma gibt warnings: untersuchen
t.r <- regr(Time~Type1+Type2,data=t.d,family=Gamma)
t.r <- regr(Time~Type1+Type2,data=t.d,family=Gamma(link=identity))  ## !!! warning

## ===========================================================
##- multinom
t.d <- d.umweltumf
##- r.rm <- multinom(Hauptv~Alter + Schule + Beeintr + Geschlecht, data=t.d)
# t.r <- multinom(Hauptv~Alter + Schule + Beeintr + Geschlecht, data=t.d)
t.r <- regr(Hauptv~Alter + Schule + Beeintr + Geschlecht, data=t.d,
             family="multinomial", na.action=na.omit) 
t.d <- d.umweltumf
t.d <- t.d[1:100,]
t.r <- regr(Hauptv~Alter + Schule + Beeintr + Geschlecht, data=t.d)
t.rm <- multinom(Hauptv~Alter + Schule + Beeintr + Geschlecht, data=t.d)
t.r <- regr(Hauptv~Alter + Schule + Beeintr + Geschlecht, data=na.omit(t.d))

t.r <- regr(Hauptv~Alter, data=t.d) ## ??? coeftable yields NA = drop1.multinom
t.dr <- drop1(t.r) 

##- t.ra <- t.rm
##- t.ra$resid <- t.ra$resid[,1]
##- plot.regr(t.ra)
## plfitpairs(t.r) ## !!! $y not available
## ===================================================================
##-   source("ftp://stat.ethz.ch/U/jaggi/NDK/R/umweltdump5.R")  # -> t.d
##-   r.g1 <- regr(cbind(ja,nein)~Schule, data=t.d,
##-                family="binomial")
##-   r.g1

  ##- r.g1 <- glm(cbind(ja,nein)~Schule, data=t.d,
  ##-             family="quasibinomial")
  ##- summary(r.g1)

##-   plot(r.g1)
##-   plot(r.g1,ask=F,xplot=T)
## ===================================================================
##-   source("ftp://stat.ethz.ch/U/jaggi/NDK/R/umweltdump5.R")
##-   source("ftp://stat.ethz.ch/U/jaggi/NDK/R/rg1-fkt.R")
##- 
##-   # glm
##-   r.gl <- glm(cbind(ja,nein)~Schule+Ortsgroesse, data=t.d,
##-                family="binomial")
##- 
##-   # regr
##-   r.fr <- regr(cbind(ja,nein)~Schule+Ortsgroesse, data=t.d,
##-                family="binomial")
##- 
##- 
##-   drop1(r.fr,test="Chisq")
## ===================================================================
t.d <- d.reacten[300:1700,]
r.lin <- regr(log10(q)~time,data=t.d, method="lm")
t.cl <- unname(coefficients(r.lin))
## !!! termtable
t.r <- regr(q~theta1*exp(-theta2*time), data=t.d, nonlinear=T,
             start=c(theta1=10^t.cl[1],theta2=t.cl[2]))
##- t.r <- nls(q~theta1*exp(-theta2*time), data=t.d, start=c(theta1=10^t.cl[1],theta2=t.cl[2]))
example(nls)
t.d <- Treated
t.r <- regr(~weighted.MM(rate, conc, Vm, K), data = t.d, nonlinear=T,
       start = list(Vm = 200, K = 0.1))  
plot.regr(t.r)
## ===================================================================
##- options(step="manova")
##- t.d <- read.table("/u/stahel/data/fossilien-samples.dat")
##- t.y <- as.matrix(t.d[,c("Length","rWidth","rCLength","Cratio")])
##- r.man <- manova(t.y~SST.Mean+Chlorophyll+Salinity,data=t.d)
##- summary(r.man, test="Wilks")
##- class(r.man)
##- t.ra <- aov(t.y~SST.Mean+Chlorophyll+Salinity,data=t.d)
##- # summary(t.ra)
##- drop1(t.ra)
##- r.lm <- lm(t.y~SST.Mean+Chlorophyll+Salinity,data=t.d)
##- # summary(r.lm)
##- plot.regr(r.lm)
##- t.url <- "http://stat.ethz.ch/Teaching/Datasets/NDK/fossilien2.dat"
##- d.fos <- read.table(t.url, header=T)
##- ##- names(d.fos) ; dim(d.fos)
##- d.foss <- d.fos
##- dump("d.foss","/u/stahel/data/d.foss.R")

source("/u/stahel/data/d.foss.R")

### Multivariate Regression
t.y <- d.fossiles[,c(35,36,32)]        # Zielvariablen
##- r.mlm <- lm(as.matrix(t.y)~SST.Mean+Salinity+lChlorophyll+region+N, data=d.foss)
##- r.man <- manova(as.matrix(t.y)~SST.Mean+Salinity+lChlorophyll+region+N,
##-                 data=d.foss)
##- summary.mreg(r.mlm)
r.mregr <- regr(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+region+N,
                data=d.fossiles)
r.mregr
plot(r.mregr, mfrow=c(3,3),ask=c.ask, lab=1)
drop1.mlm(r.mregr)
### Residuenanalyse
# plot.regr(r.mlm,mfrow=c(1,3),ask=c.ask)
# t.r <- regr(as.matrix(t.y)~1,data=d.foss) ## does not work !!!
## ===================================================================
require(MASS)
userOptions(step="polr")

data(housing, package="MASS")
t.r <- regr(Sat ~ Infl + Type + Cont, weights = housing$Freq,
            data = housing)
plot(t.r)

##- t.r <- polr(Beeintr~Alter+Geschlecht, data=d.umweltumf)
##- summary(t.r)
##- t.d <- d.umweltumf
##- t.r <- polr(Beeintr~Alter+Geschlecht, data=d.umweltumf) 
##- t.t <- ftable(Beeintr~Schule,data=d.umweltumf)
#t.r <- polr(t.t~x, data=data.frame(x=1:5)) # error!
t.r <- regr(Beeintr~Alter+Geschlecht, data=d.umweltumf) # multinom
t.d <- d.umweltumf
t.d$y <- ordered(t.d$Beeintr)
dd <- t.d[seq(1,nrow(t.d),5),]
dd$random <- c(rep(NA,20),rnorm(nrow(dd)-20))
## !!! attach in i.polr vermeiden! siehe i.mlm
t.fo <- y~Alter+Geschlecht+random
t.vars <- all.vars(t.fo)
t.rr <- regr(t.fo, data=dd)
t.rrd <- regr(t.fo, data=na.omit(dd[,t.vars]))
t.r <- polr(t.fo, data=dd)
t.rd <- polr(t.fo, data=na.omit(dd[,t.vars]))
t.rst <- step(t.rd)
t.r0 <- polr(y~Alter+Geschlecht, data=dd)
t.r0d <- polr(y~Alter+Geschlecht, data=na.omit(dd[,t.vars]))
add1(t.r0,scope="random")
add1(t.r0d,scope="random")
t.rr0 <- regr(y~Alter+Geschlecht, data=dd, subset=seq(1,nrow(dd),2))
t.rr0d <- regr(y~Alter+Geschlecht, data=na.omit(dd[,t.vars]))
add1(t.rr0,scope="random")
add1(t.rr0d,scope="random")

r.rp <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
##- summary(r.rp)
##- dd <- housing
##- dd$Cont
##- t.r <- regr(Sat ~ Infl + Type + Cont, weights = housing$Freq, data = housing) 

f <- function()
{
  t.dat <- housing
     r.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = t.dat)
     try(summary(r.plr))
  t.dat$Satnum <- as.numeric(t.dat$Sat)
     r.lm <- lm(Sat ~ Infl + Type + Cont, weights = Freq, data = t.dat)
     summary(r.lm)
}

t.d <- d.initiative <-
  read.table("http://stat.ethz.ch/Teaching/Datasets/NDK/initiat.dat",
                  sep="\t",header=T)
r.mod1 <- regr(y~sex+polit,data=t.d,family="gaussian")

drop1(r.mod1)
## ---------------------------------------
require(survival)
# d.ovarian <- ovarian
t.rs <- survreg(formula = Surv(futime, fustat) ~ ecog.ps + rx,
                data = ovarian, dist = "weibull")
t.rs <- survreg(formula = Surv(futime, fustat) ~ ecog.ps + rx,
                data = ovarian, dist = "weibull",scale=2)
t.rss <- summary(t.rs)
##- t.rs <- survreg(formula = Surv(log(futime), fustat) ~ ecog.ps + rx,  
##-     data = ovarian, dist = "extreme") ## not the same
t.r <- regr(formula = survival::Surv(futime, fustat) ~ ecog.ps + rx,
            data = ovarian, family="weibull")
plot(t.r)
t.rc <- coxph(formula = survival::Surv(futime, fustat) ~ ecog.ps + rx,
              data = ovarian)
t.r <- regr(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian)

str(d.kevlar)
d.kevlar$Spool <- factor(d.kevlar$Spool)
d.kevlar <- d.kevlar[d.kevlar$Stress>=24,]
rregr <- regr(Surv(Failure,rep(1,87))~Stress+Spool, data=d.kevlar, 
              family="weibull")
t.sr <- survreg(Surv(Failure,rep(1,87))~Stress+Spool, data=d.kevlar, 
              dist="weibull")

bladder1 <- bladder[bladder$enum < 5, ]
t.cph <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + 
               cluster(id), bladder1) 
t.r <- regr(Surv(stop, event) ~ rx + size + number, bladder1) 
t.r <- regr(Surv(time,status)~x, data=aml, family="weibull")
showd(residuals(t.r))

## 
t.d <- d.cheese
t.d$Anzahl[1] <- NA
t.d$Bakt[2] <- NA
dim(na.omit(t.d))
t.d$y <- Tobit(t.d$Anzahl, 10, log=T)
t.r <- regr(y~Temp+Bakt+Konz, data=t.d)
t.r <- regr(Tobit(Anzahl, 10, log=TRUE)~Temp+Bakt+Konz, data=t.d,
            family="gaussian")
plot(t.r)
t.d$y <- Tobit(t.d$Anzahl, 10, transform=sqrt)
t.rf <- fitted(t.r)
t.rp <- predict(t.r)
t.rr <- residuals(t.r)

add1(t.r)
step(t.r)
## margEffTobit(t.r)
confint(t.r)

t.rr <- residuals(t.r)
t.i <- which(t.rr[,"prob"]>0)
t.o <- order(t.rr[t.i,"fit"])
t.i <- t.i[t.o]
matplot(t.rr[t.i,"fit"],t.rr[t.i,c(1:3)],type="l")
t.rr[t.i,"fit"]
rr <- t.rr[t.i,]

example(survreg)
data("tobin", package="survival")
Tobit(tobin$durable, log=T)
t.r <- regr(Tobit(durable) ~ age + quant, data=tobin)
plot(t.r)
t.r <- regr(Tobit(durable, log=T) ~ age + quant, data=tobin)
plot(t.r)
## ----------------------------------------------------------
## quantile regr
t.r <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=d.spreng,
            method="rq", tau=0.7)
##- t.r <- rq(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=d.spreng,
##-           tau=0.7)
##- drop1(t.r, test="Chisq")
## -----------------------------------------------------------
showd(d.fossiles)
plmatrix(cbind(log10(Petal.Width),Petal.Length)~., data=iris)
plmatrix(~Sepal.Length+Species, ~log10(Petal.Width)+Petal.Length, data=iris)
plmatrix(~Sepal.Length+Species, ~log10(Petal.Width)+I(1:150), data=iris)
plmatrix(cbind(log10(Petal.Width),Petal.Length)~., data=iris,
         col=as.numeric(Species)+1)
plmatrix(cbind(log10(Petal.Width),Petal.Length)~Species, data=iris,
         col=as.numeric(Species)+1, xaxmar=1)
plmatrix(log10(tremor)~as.factor(device), data=d.blast, col=location)

mframe(1,2)
plmboxes(tremor~location, data=d.blast)
plmboxes(tremor~location, data=d.blast, ilimext=0.1)
plmboxes(tremor~location, data=d.blast, ilim=c(0,20),ilimext=0.05)
plmboxes(tremor~location, data=d.blast, ilim=c(0,15))
plmboxes(tremor~location+I(charge>5), data=d.blast)
plmboxes(tremor~1+location, data=d.blast)
plmboxes(tremor~location, data=d.blast, at=c(1,3,4,7,NA,10,11,12),width=c(1,0.5))

plmboxes(len~dose+supp, data=ToothGrowth, colors=c(med="red"),
         probs=c(0.1,0.25))
boxplot(len~dose+supp, data=ToothGrowth)

     utils::data(anorexia, package = "MASS")
plmboxes( I(Postwt-Prewt)~Treat , data=anorexia)
t.wtmed <- median(anorexia$Prewt)
d.anorexia <- anorexia
d.anorexia[sample(nrow(d.anorexia),5),"Postwt"] <- NA
plmboxes( I(Postwt-Prewt)~Treat + I(Prewt>t.wtmed) , data=d.anorexia,
         probs=c(0.1,0.25), na=T)


plot( I(Postwt-Prewt)~as.numeric(Treat), data=anorexia)

set.seed(2)
t.x <- sort(round(2*rchisq(20,2)))
table(t.x)
t.p <- ppoints(100)
plot(quinterpol(t.x,t.p),t.p, type="l")
lines(c(0,t.x), 0:length(t.x)/length(t.x) , type="s")


showd(wagnerGrowth)
t.rr <- regr(y~PA*Region, data=wagnerGrowth, contrast="contr.treatment")
t.rl <- lmrob(y~PA*Region, data=wagnerGrowth)


