source('/u/stahel/R/regdevelop/pkg/regr0/R/regr.R')
options(digits=3)
load("/u/stahel/data/regdata")
load("../data/d.blast.rda")
load("../data/d.surveyenvir.rda")
load("../data/d.rehab.rda")
load("../data/d.fossiles.rda")

##- source("/scratch/users/stahel/transfer/regr0/R/regr.R")
##- source("/scratch/users/stahel/transfer/regr0/R/drop1.R")
##- options(digits=3)
##- load("/scratch/users/stahel/data/regdata")

##- ldir <- "/u/stahel/R/Pkgs/regr0/man"
##- for (lnm in c("datatype","simresiduals", "d.rehab","d.fossiles"))
##-   prompt(lnm,paste(ldir,"/",lnm,".Rd",sep=""))
## ---------------------------------------------------------------
## --- regr.Rd examples
# data(d.blast)
options(verbose=1)

( r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast) )
plot(r.blast)

options(project='regr0.demo', step='blast')
u.pslatex('p-plotregr-ta')
  options(colors.ra =
          c("black","gray4","blue4","blue3","darkgreen","green",
            "burlywood4","burlywood3","burlywood4"))
par(lwd=2)
plot(r.blast, plotselect=c(ta=3), xplot=F, seq=F, pch=1, mf=1,
     lwd=c(2,1.5,2,1.5,1.5,1,1,1,1))
ps.end()
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
t.r <- regr(disturbance~age+education+location, data=d.surveyenvir)

## ordered regression
t.r <- regr(Sat ~ Infl + Type + Cont, weights = housing$Freq, data = housing)
plot(t.r)

## multivariate regression
data(d.fossiles)
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
     c( qq = NA, yfit=2, ta=3, tascale = NA, weights = NA, hat = 0),ask=c.ask)
t.r <- regr(y~f1+x1+x2+f1:x2, data=t.d)

data(anorexia, package="MASS")
r.anorex <- regr(Postwt ~ Prewt + Treat + offset(Prewt),
                 family = gaussian, data = anorexia)

r.savings <- regr(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
plot(r.savings, ask=c.ask)
## Spreng
t.d <- d.spreng
## t.d <- t.d[t.d$stelle<=2,]
t.r <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=d.spreng)
t.form <- log10(ersch)~Stelle+log10(dist)+log10(ladung)
t.r <- regr(t.form, data=t.d)
t.r <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=t.d,
             weights = 5+1:nrow(t.d))
plot(t.r,ask=c.ask)
plresx(t.r,var=names(t.d),ask=c.ask)

t.d$ersch[3:5] <- NA
t.r <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=t.d,
             weights = 5+1:nrow(t.d))  
plot(t.r,ask=c.ask)
drop1(t.r)
t.r <- regr(log10(ersch)~1, data=d.spreng)
t.r1 <- step(t.r, scope=~.+Stelle+log10(dist)+log10(ladung),trace=F)
t.r2 <- update(t.r1, formula=.~.-log10(ladung))
t.r3 <- update(t.r2, formula=.~.-Stelle)
##- t.mt <- modelTable(c("t.r","t.r1"))
t.mt <- modelTable(list(large=t.r1,reduced=t.r2,small=t.r3))
t.mt[,-2]
compareTerms(large=t.r1,reduced=t.r2,small=t.r3)
## -------------------------------------
## robust
t.r <- regr(log10(ersch)~Stelle+log10(dist)+log10(ladung), data=d.spreng,
            robust=T)
## -------------------------------------

## Baby Survival
##- t.d <- read.table("/u/stahel/data/babysurvival.dat",sep=",",header=T)
##- t.d <- t.d[t.d$Age<35,]
##- d.babysurv <- t.d
t.d <- d.babysurv
t.r <- regr(Survival~.,data=t.d)
t.rglm <- glm(Survival~.,data=t.d,family=binomial)
t.rs <- step(t.r, trace=F)  ## ???
t.r <- r.babysurv <- regr(Survival~Weight+Age+Apgar1,data=t.d,family=binomial)
plot(r.babysurv,xplot=~Weight,cex=0.7,symbol.size=NULL,res.lim=c(-5,5))
plot(r.babysurv,glm.restype="cond", ask=c.ask)

mframe(2,2)
plresx(r.babysurv,vars=~Age+Apgar1+Apgar5+pH,data=d.babysurv,
       weights=F,cex.lab=0.2, ask=c.ask)

plresx(t.r,data=t.d,vars=~.+Apgar5,sequence=T, ask=c.ask)

##- t.d <- read.table("/u/stahel/data/babysurvival-w.dat",sep=",",header=T)
##- t.d <- t.d[,-1]
##- d.babysurv.w <- t.d
t.d <- d.babysurv.w
t.r <- regr(cbind(Survival.1,Survival.0)~Weight,data=t.d,family=binomial)
t.r <- glm(cbind(Survival.1,Survival.0)~Weight,data=t.d,family=binomial)
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
##- t.r <- nls(q~theta1*exp(-theta2*time), data=t.d, 
##-              start=c(theta1=10^t.cl[1],theta2=t.cl[2]))
example(nls)
t.d <- Treated
t.r <- regr(~weighted.MM(rate, conc, Vm, K), data = t.d, nonlinear=T,
       start = list(Vm = 200, K = 0.1))  
plot(t.r)
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
plot(r.mregr, mfrow=c(3,3),ask=c.ask)
drop1.mlm(r.mregr)
### Residuenanalyse
# plot.regr(r.mlm,mfrow=c(1,3),ask=c.ask)
# t.r <- regr(as.matrix(t.y)~1,data=d.foss) ## does not work !!!
## ===================================================================
require(MASS)
options(step="polr")
##- t.r <- polr(Beeintr~Alter+Geschlecht, data=d.umweltumf)
##- summary(t.r)
##- t.d <- d.umweltumf
##- t.r <- polr(Beeintr~Alter+Geschlecht, data=d.umweltumf) 
##- t.t <- ftable(Beeintr~Schule,data=d.umweltumf)
#t.r <- polr(t.t~x, data=data.frame(x=1:5)) # error!
t.r <- regr(Beeintr~Alter+Geschlecht, data=d.umweltumf) # multinom
t.d <- d.umweltumf
t.d$y <- ordered(t.d$Beeintr)
## !!! attach in i.polr vermeiden! siehe i.mlm
t.r <- regr(y~Alter+Geschlecht, data=t.d) 

##- r.rp <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
##- summary(r.rp)
t.r <- regr(Sat ~ Infl + Type + Cont, weights = housing$Freq, data = housing) 

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
r.mod1 <- regr(y~sex+polit,data=t.d,family="gaussian",method="lm")

drop1(r.mod1)
## ---------------------------------------
require(survival)
# d.ovarian <- ovarian
t.rs <- survreg(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, 
    dist = "weibull")
t.rs <- survreg(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, 
    dist = "weibull",scale=2)
t.rss <- summary(t.rs)
##- t.rs <- survreg(formula = Surv(log(futime), fustat) ~ ecog.ps + rx,  
##-     data = ovarian, dist = "extreme") ## not the same
t.r <- regr(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = d.ovarian,
            family="weibull")
plot(t.r)
t.rc <- coxph(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = d.ovarian)
t.r <- regr(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = d.ovarian)

bladder1 <- bladder[bladder$enum < 5, ]
t.cph <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + 
               cluster(id), bladder1) 
t.r <- regr(Surv(stop, event) ~ rx + size + number, bladder1) 

## d.cmbscores <- t.d
t.d <- d.cmbscores
dim(na.omit(t.d))
t.d$y <- Tobit(t.d$EVAPOR, 10, log=T)
t.r <- regr(y~Temp+Time+lWindspeed, data=t.d)
plot(t.r)

