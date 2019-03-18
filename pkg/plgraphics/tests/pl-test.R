
require(plgraphics) ##, lib="/u/stahel/R/regdevelop/pkg/plgraphics.Rcheck")

plyx(Sepal.Width~Sepal.Length, data=iris)
plyx(iris[,c("Sepal.Width","Sepal.Length")]) ##!!! farben

## ploptions
ploptions("linewidth")
t.plo <- ploptions(linewidth=1.5)
ploptions("linewidth")
t.plo$linewidth
.ploptions$linewidth
ploptionsDefault$linewidth
t.plo <- ploptions(default="linewidth")
ploptions("linewidth")

##!!!
t.plo <- ploptions(col="magenta", smooth.col="darkgreen", assign=F)
attr(t.plo, "old")
ploptions("col")
t.plo$col

plmfg(2,2)
plyx(Sepal.Width~Sepal.Length, data=iris, ploptions=t.plo)
plyx(Sepal.Width~Sepal.Length, data=iris, 
     psize=Petal.Length^3, pcol=Species, pch=Species, cex=1.5)
plyx(Sepal.Width~Sepal.Length, data=iris, smooth=2, smooth.group=Species)
plyx(Sepal.Width~Sepal.Length, data=iris, smooth=TRUE, group=Species)
plmfg()
plyx(jitter(Sepal.Width) ~ jitter(Sepal.Length), data=iris, axp=7, plab=T)
plmfg(2,3, mar=c(NA, 0.5), oma=c(2,2,2,2)+2)
plyx(Petal.Length+Petal.Width~Sepal.Length+Sepal.Width, group=Species,
     data=iris)
plmfg(2,2)
plyx(Petal.Length ~ Sepal.Length+Sepal.Width, data=iris, smooth=TRUE,
     smooth.group=iris$Species, reflines=lm)
plyx(Sepal.Width~Sepal.Length, data=iris[1:50,], smooth=F, markextremes=0.1)
attr(iris$Sepal.Length, "ticksat") <- 
  structure(seq(4, 8, 0.5), small=seq(4,8,0.1))
iris$"(pcol)" <- as.numeric(iris$Species)
plyx(Sepal.Width~Sepal.Length, data=iris)

t.plargs <- pl.control(~Species+Petal.Length, ~Sepal.Width+Sepal.Length,
                       data=iris, smooth.group=Species, pcol=Species)
t.plargs$ploptions$group.col <- c("magenta","orange","cyan")
plpanel(iris$Petal.Length, iris$Petal.Width, plargs=t.plargs, frame=TRUE)

t.plo <- ploptions(col="blue")
plyx(Sepal.Width~Sepal.Length, data=iris, ploptions=t.plo)
plyx(Sepal.Width~Sepal.Length, data=iris)

ploptions(gridlines.col="lightblue")
t.plo <- ploptions(list=list(smooth.lty=4, smooth.lwd=5), assign=FALSE)
plyx(Sepal.Width~Sepal.Length, data=iris, ploptions=t.plo, gridlines=TRUE)

plyx(y=EuStockMarkets[1:40,], type="b") ## ??? 2 blaue linien
plyx(structure(1:40, varlabel="time"), EuStockMarkets[1:40,], type="b")

ff <- function(formula, data, smooth=T, pcol=1)
  plyx(formula, data=data, smooth=smooth, pcol=pcol)
ff(Sepal.Width~Sepal.Length, data=iris, pcol=I("gray"), smooth=T)

## plmatrix
plmatrix(iris, pch=as.numeric(Species))
plmatrix(~Sepal.Length+Sepal.Width, ~Petal.Length+Petal.Width, data=iris,
         smooth=TRUE, pch=as.numeric(iris[,"Species"]))

## plmboxes
plmboxes(Sepal.Width~Species, data=iris, labelsvert=1, main="iris")
plmboxes(Sepal.Length~Species, data=iris,
  widthfac=c(med=2), colors=c(med="red"), horizontal=TRUE)

## attributes of variables
data(d.blast)
dd <- genvarattributes(d.blast)
str(attributes(dd$tremor))
ddd <- varattributes(dd, list( tremor=list(axisat=seq(0,24,2),
  axislabels = seq(0,24,10)) ) )
str(attr(ddd$tremor, "axislabels"))

## =========================================
#require(regr0)
## attach("../div/pl-data.rda")
data(d.blast)
rr <- r.blast <-
  lm(logst(tremor)~location+log10(distance)+log10(charge), data=d.blast)
plregr(rr)
plregr(rr, xvar=FALSE)
plregr(rr, transformed=TRUE, reflinesband=TRUE)

dd <- d.blast[as.numeric(d.blast$location)<=3,]
dd[1,"distance"] <- 200
rr <- lm(log10(tremor)~log10(distance)+log10(charge)+location, data=dd)
plres2x(~ log10(distance) + log10(charge), reg=rr, transformed=F,
        pcol=location) ## ???

## utilities
showd(dd)
sumNA(dd)
tit(dd) <- "blasting"
plmatrix(dd, main="test plmatrix")

## functions generating elements
t.fc <- fitcomp(rr,se=TRUE)
t.fc$comp[1:10,]
t.fct <- fitcomp(rr, se=TRUE, transformed=TRUE)

rr <- lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
r.smooth <- gensmooth( fitted(rr), residuals(rr))
showd(r.smooth$y)
plot(fitted(rr), resid(rr), main="Tukey-Anscombe Plot")
abline(h=0)
lines(r.smooth$x,r.smooth$y, col="red")
## grouped data
t.plargs <- list(pldata=data.frame(d.blast$location), names="(smoothGroup)")
## residuals against regressor, without  plresx:
t.res <- naresid(structure(r.blast$na.action, class="exclude"), residuals(r.blast))
r.smx <- gensmooth( d.blast$dist, t.res, plargs=t.plargs)
plot(d.blast$dist, t.res, main="Residuals against Regressor")
abline(h=0)
plsmoothline(r.smx, d.blast$dist, t.res, plargs=t.plargs)

## --------------------------------------------------------
## multivariate regression
data(d.fossileSamples)
rr <-
  lm(cbind(sAngle,lLength,rWidth)~SST+Salinity+lChlorophyll+Region,
                data=d.fossileSamples)
plregr(rr)

## ================================================
## glm
data(d.babysurvival)
t.d <- d.babysurvival
t.d$Age[2] <- NA
rr <- glm(Survival~Weight+Age+Apgar1,data=t.d,family=binomial)
plregr(rr, xvar=~Weight, cex.plab=0.7, ylim=c(-5,5))
plregr(rr, condquant=FALSE)

## polr
data(housing, package="MASS")
rr <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
t.res <- residuals.regrpolr(rr)
showd(attr(t.res, "condquant"))
plregr(rr)
plregr(rr, factor.show="jitter")

## survreg
data(lung, package="survival")
lung$gender <- factor(c("m","f")[lung$sex])
r.sr <- survival::survreg(Surv(time, status) ~ age + gender + ph.karno, data=lung) 
plregr(r.sr, group=gender, pcol=gender, xvar=~age)
r.cox <- survival::coxph(Surv(time, status) ~ age + gender + ph.karno, data=lung) 
plregr(r.cox, group=gender, pcol=gender, xvar=~age)

## =================================================================
## gendate
rr <- gendate(year=2010, month=c("Jan","Apr"), day=c(3,30),
            hour=25, min=c(0,70), sec=c(0,300))
stopifnot(all(
  print(rr) == c("2010-01-04 01:00:00 CET", "2010-05-01 02:15:00 CEST")
  ))
td <- data.frame(datum=as.Date(c("2010-05-20","1968-05-01")),
                 tag=c(1.5, 3), min=c(30,70))
rr <- gendate(date=datum, day=tag, hour=6, data=td, min=min)
stopifnot(all(
  print(rr) == c( "2010-05-02 06:30:00 CEST",  "1968-05-03 07:10:00 CET")
  ))
gendate(day=tag, hour=4, data=td, min=min) 
gendate(day=tag, data=td, min=min, sec=8) 
