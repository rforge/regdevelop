##  regr.R  Functions that are useful for regression, W. Stahel,
## ==========================================================================
regr <- function(formula, data=NULL, tit=NULL, family=NULL, # dist=NULL,
                 calcdisp=NULL, suffmean=3,
                 nonlinear = FALSE, start=NULL,
                 robust = FALSE, method=NULL, init.reg="f.ltsreg",
                 subset=NULL, weights=NULL, na.action=nainf.exclude,
                 contrasts=getUserOption("regr.contrasts"),
                 model = FALSE, x = TRUE, termtable=TRUE, vif=TRUE, ...)
{
  ## !!! dispersion: allow to be set.
  ## Purpose:    fit all kinds of regression models
  ## -------------------------------------------------------------------------
  ## Arguments:
  ##   formula, data, ...  as with lm
  ##   tit        title (becomes tit attribute of result)
  ##   calcdisp   should dispersion be calculated for
  ##              family=binomial and family=poisson
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  Jan 02
  ## -------------------------------------------------------------------------
  ## --- preparations
  lcall <- match.call()
  ## b. --- convert character formula to formula
  lformula <- as.formula(formula)
  lcall$formula <- lformula
  ## d. --- data
  dataname <- deparse(substitute(data))
  data <- as.data.frame(data)
    lnobs <- nrow(data)
  if (lnobs==0&&ncol(data)>0) stop("!regr! no observations in data")
  ## f. --- names of all variables to be used (as in lm)
  ## get all vars
  mf <- match.call(expand.dots = FALSE)
  lj <- match(c("formula", "data", "weights", "offset"), # "subset",
             names(mf), 0L)
  mf <- mf[c(1L, lj)]
  mf$formula <- lformula
  mf[[1L]] <- as.name("get_all_vars")
  if (!(is.null(nonlinear)||as.character(nonlinear)=="FALSE")) {
    lcf <- setdiff(all.vars(formula),names(data))
    if(nonlinear=="lazy") {
      lcfundef <- setdiff(lcf, names(start))
      if (length(lcfundef)) {
        warning(":regr(nonlinear): check if  ",paste(lcfundef, collapse=", "),
                "  are coefficients. Starting values will be 0.")
        lst <- rep(0,length(lcfundef))
        names(lst) <- lcfundef
        start <- c(start,lst)
      }
    }
    lcall$start <- start <- start[names(start)%in%lcf]
    lfo <- setdiff(all.vars(formula),names(start))
    mf$formula <- as.formula(paste("~",paste(lfo,collapse="+")))
    lcall$nonlinear <- nonlinear <- TRUE
  } else {
    mf$formula <- lformula
  }
  lav <- try(eval(mf, parent.frame()))
  if (class(lav)=="try-error") stop("!regr! undefined variables in formula")
  ## convert character to factor
  for (lvn in 1:ncol(lav)) {
    lv <- lav[[lvn]]
    if (is.character(lv)|is.factor(lv)) lav[[lvn]] <- factor(lv)
  }
  ## -------------------------------------------
  ## h. --- response type
  if (length(lformula)==2) { # nonlinear called with formula of type ~...
    ly <- rep(0,NROW(lav))
    lytype <- "numeric"
  } else {
    lyf <- model.frame(lformula[1:2], lav)
    ltrm <- attr(lyf, "terms")
    lytype <- substring(attr(ltrm, "dataClasses"),1,5)
    lysimple <- lytype!="nmatr" ## not a matrix
    ly <- lyf[[1]]
    if (lysimple&&length(unique(ly))==2 &&
        (is.factor(ly[[1]]) || all(ly%in%0:1)))
        lytype <- "binary"
    if (inherits(ly,"Surv"))  {
        lytype <- "survival"
        require(survival)
    }
## strange variables
##-   l1v <- sapply(ldta, function(x) all(x==c(x[!is.na(x)],0)[1],na.rm=TRUE) )
##-                                 ## covers case of several or all NAs
##-   if (any(l1v)) {
##-     warning(paste(":regr: variable(s)", paste(lvars[l1v],collapse=", "),
##-                   "has (have) no distinct values")) #  -> dropped.
##-   }
  }
  ## k. --- family and fitting function
  lfam <- as.character(substitute(family))[1]
  if (is.null(lfam)||is.na(lfam)) 
    lfam <- switch(substring(lytype,1,5),
                   numer="normal", nmatr="normal", binar="binomial",
                   binco="binomial", order="cumlogit",
                   facto="multinomial", survi="ph", "unknown")
  if (lytype=="survival")
      lfam <- c( attr(ly,"distribution"), lfam)[1]
  else  if (substring(lfam,1,7)=="multinom") lfam <- "multinomial"
  ##
  lfitfun <-
      switch( lfam,
             gaussian="lm", normal="lm", binomial="glm", poisson="glm",
             Gamma="glm",
             cumlogit="polr", multinomial="multinomial",
             weibull="survreg", lognormal="survreg", loggaussian="survreg",
             loglogistic="survreg", extreme="survreg", ph="survreg",
             prop.hazard="survreg",
             "unknown")
  if (lfitfun=="unknown") stop("!regr! Fitting function not identified")
  ## additional checks
  if (lytype=="survival") {
    if (!inherits(ly,"Surv"))
      stop("!regr! bug: convert response to Surv object")
    ## !!! hier machen! lav[,1] ersetzen durch Surv davon
    lfitfun <- "survreg"
    if (is.null(family)) lfam <-attr(ly,"distribution")
  }
  else  if (lfitfun=="glm")
    lcall$control <- list(calcdisp=calcdisp, suffmean=suffmean,lcall$control)
  ## 
  lfitname <- paste("i",lfitfun,sep=".")
  if (!exists(lfitname)||!is.function(get(lfitname)))
    stop (paste("!regr! Fitting function",lfitname, "not found"))
  ## m. --- prepare call
  lcl <- lcall
  lcl[[1]] <- ## hack --> eval(.) works also when call is source()d ...
      switch(lfitname,
	     "i.lm" = quote(regr0::i.lm),
	     "i.glm" = quote(regr0::i.glm),
	     "i.multinomial" = quote(regr0::i.multinomial),
	     "i.polr" = quote(regr0::i.polr), 
	     "i.smooth" = quote(regr0::i.smooth), ## ??
	     "i.survreg" = quote(regr0::i.survreg),
	     ## default:
	     as.name(lfitname))
##  lcl[[1]] <- as.name(lfitname) ## sonst geht das debuggen nicht.
  lcl$fname <- lfam
  lcl$na.action <- substitute(na.action)
##  lcl$data <- as.name("lav") ## environment(formula)
  lcl$data <- eval(lcl$data, sys.parent())
  ## problem with environment if different for  data  and  formula
  old.opt <- NULL
  if(is.atomic(contrasts)&&length(contrasts)) {
      if(!is.character(contrasts))
          warning("!regr! invalid contrasts argument") else {
              old.opt <- options(contrasts=c(contrasts,"contr.poly")[1:2])
              lcl$contrasts <- NULL
          }
  }
## --------------------------------------------
  lreg <- eval(lcl, envir=environment(formula))
  ## --------------------------------------------
  if (length(old.opt)) options(old.opt)
  if (is.null(lreg$distrname)) lreg$distrname <- lfam
##  <<<<<<< .mine
#  lreg$Y <- data.frame(ly) # ly is a model.frame
#  if (ncol(lyy)>1) colnames(lyy) <- colnames(ly[[1]])  
#  lreg$Y <- lyy
## =======  !?!
  lyy <- as.matrix(ly) # ly is a model.frame
  if (ncol(lyy)>1) colnames(lyy) <- colnames(ly[[1]])
  lreg$Y <- lyy
  lreg$response <- ly
  if (nonlinear) lreg$r.squared <- 1-lreg$sigma^2/var(ly)
  ## >>>>>>> .r32
  lreg$allvars <- lav ## needed more than $model
             ## since $model contains transformed variables
  lreg$funcall <- lreg$call
  lcall$formula <- formula(lreg) # hope this never damages anything
  lreg$call <- lcall
  tit(lreg) <- if (length(tit)==0) attr(data,"tit") else tit
  doc(lreg) <- attr(data,"doc")
  if (model&&length(lreg$model)==0) {
      if (nonlinear) warning(":regr: no $model available for nonlinear regr.")
      else lreg$model <- lm(lformula, data, method="model.frame")
  }
  lterms <- if (nonlinear) NULL else terms(lreg)
  if ((!nonlinear) && is.null(attr(lterms, "predvars")))  ## needed for survreg
    attr(lreg$terms,"predvars") <- attr(attr(lreg$model,"terms"),"predvars")
  lresnm <- colnames(lreg$residuals)
  ## r. --- leverages, standardized res
  lhat <- lreg$leverage
  if (length(lhat)==0 && !nonlinear) {
    lmm <- lreg[["x"]]
    if (length(lmm)==0) lmm <- model.matrix(lterms,lreg$model)
    lreg$leverage <- lhat <- hat(lmm)
    if (length(lhat)==0) warning(":regr: no leverages available")
  }
  if (length(lhat)==NROW(lreg$stres))
    lreg$stres <- lreg$stres/ifelse(lhat>=1,1,sqrt(pmax(1e-10,1-lhat))) else
      if (length(lreg$stres))
        warning(":regr: bug: leverages and st.res. incompatible")
  if (class(lreg)[1]=="survreg") lreg$n.obs <- length(lreg$linear.predictor)
  ## misc
  if (is.null(x) || !x) lreg$x <- NULL
  class(lreg) <- if (class(lreg)[1]=="orig")  ##  nls shall not be regr
    class(lreg)[-1] else c("regr",class(lreg))
## result of regr
  lreg
}
## -----------------------------------------------------------------------
i.lm <- function(formula, data, family, fname="gaussian", nonlinear=FALSE,
                 robust=FALSE, method=NULL, control=NULL, 
                 vif=TRUE, termtable=TRUE, ...)
{
  ## Purpose:  internal: fit lm
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
  ## b. --- method
  lcall <- match.call()
  lmeth <- c(lcall$method,"")[1]
  lfn <- if (nonlinear) "nls" else {
    if (lmeth=="rlm") robust <- TRUE
    if (robust) {
      if (is.null(method)) method <- lcall$method
      if (is.null(method)) method <- c("lmrob","KS")
      method[1]
    } else "lm"}
  if (robust) {
    if (lfn=="lmrob") {
 ##     require(robustbase)  ## !?!
      if (length(method)>1) {
        if (substr(method[2],1,2)=="KS")
          method <- NULL
        lcall$setting <- "KS2011"
      }
      lcall$x <- TRUE
  }
    if (lfn=="rlm") {
##      require(MASS)    ##  !?!
      lcall$method <- c(method,"MM")
      lcall$x.ret <- TRUE
      lcall$robust <- NULL
    }
  } else  lcall$x <- TRUE
  if (lmeth=="rq"|lmeth=="quantreg") { # quantile regression
##    require(quantreg) ## !?!
    lfn <- "rq"
    lcall$method <- if(length(lcall$method)>1) lcall$method[-1] else NULL
    lcall$x <- NULL
  }
  ## d. --- call
  lcall$method <- if (length(method)>1) method[-1] else NULL
  ##  method[-1]  produces character(0) which is not NULL!
  mkFn <- function(fn) { ## hack --> eval(.) works also when call is source()d ...## ??? wieso function?
      switch(fn,
	     "lmrob" = quote(robustbase::lmrob),
	     "rlm" = quote(MASS::rlm),
             "rq" = quote(quantreg::rq),
	     ## default:
	     as.name(fn))
  }
  lcall[[1]] <- mkFn(lfn)
  lcall$fname <- lcall$family <- lcall$vif <- lcall$nonlinear <- NULL
  ## --------------------------
  lreg <- eval(lcall, envir=environment(formula))
  ## --------------------------
  ## f. --- collect results
  lreg$call$formula <- formula
  lreg$fitfun <- lfn
  lreg$distrname <- "gaussian"
  lttype <- switch(lfn,
                   rq="Chisq",
                   rlm="Chisq",
                   "F"
                   )
##-   ## leverage
##-   if (!nonlinear) {
##-     lhat <- pmax(0,hat(lreg$x))
##-     if (length(lhat)==0) warning(":regr/i.lm: no leverages")  ## else {
##- ##-       if (length(lhat)!=NROW(lreg$stres))
##- ##-         if (length(lreg[["w"]])==NROW(lreg$stres))
##- ##-           lhat <- u.merge(lreg$leverage, 0, lreg[["w"]]>0)
##- ##-     }
##-     lreg$leverage <- lhat
##-   }
  ## multivariate
  if (class(lreg)[1]=="mlm")
    return(i.mlmsum(lreg, termtable))
  ## 
  lreg1 <- summary(lreg)
  lsig <- lreg1$sigma
  if (is.null(lsig)) lsig <- lreg$scale ## lmrob
  if (is.null(lsig)) lsig <- sd(resid(lreg)) # !!! used for rq
  lreg$sigma <- lsig
  ## standardized residuals
  if (is.finite(lsig)&&lsig>0) {
    lreg$stres <- lreg$residuals/lsig
    if (length(lreg$weights)) lreg$stres <- lreg$stres*sqrt(lreg$weights)
  } ## standardization by lhat is done in regr
  if (class(lreg)=="lmrob") lreg1$cov.unscaled <- lreg$cov/lsig^2 ## !!!
  ## from summary
  lcomp <- c("r.squared","fstatistic","colregelation","aliased",
             "df","cov.unscaled")
  lreg[lcomp] <- lreg1[lcomp]
  if (nonlinear) {
    lreg$coefficients <- lreg1$coefficients[,1]
    lreg$testcoef <- lreg1$coefficients
#   lreg$r.squared <- 1-(lsig/lsdy)^2
  }
  ## degrees of freedom
  ldfr <- lreg$df.residual
  if (length(ldfr)==0||is.na(ldfr))
      lreg$df.residual <- ldfr <- lreg1$df[2] # rlm
  if (is.null(lreg$df)) # needed for rq
    lreg$df <- c(length(coef(lreg))-attr(terms(lreg),"intercept"),
                 length(resid(lreg))-length(coef(lreg)))
  lreg$adj.r.squared <- 1-(1-lreg$r.squared)*(length(lreg$residuals)-1)/ldfr
  ## cov of estimates
  lcov <- lreg$cov.unscaled*lsig^2
  lreg$covariance <- lcov
  lse <- sqrt(diag(lcov))
  lreg$correlation <- lcov/outer(lse, lse)
  ## --- table of terms
  if (!nonlinear) {
    if(termtable) {
      ly <- lreg$model[[1]]
      lsdy <- sqrt(var(ly))
      ltt <- i.termtable(lreg, lreg1$coef, data, lcov, lttype, lsdy=lsdy,
                         vif=vif, leverage=TRUE) 
      lcmpn <- c("testcoef","allcoef","leverage")
      lreg[lcmpn[lcmpn%in%names(ltt)]] <- ltt
    }  else  class(lreg) <- c("orig",class(lreg))
  }
  ## result of i.lm
  lreg
}
## -----------------------------------------------------------------------
i.mlmsum <- function(object, termtable=TRUE)
{
  ## Purpose:  internal: fit multivariate lm;  called from i.lm() only
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
##-   lreg <- lm(formula, data=data, weights=.Weights, model=model, ...)
##-   lreg$call$formula <- eval(lreg$call$formula) # patch
  lreg1 <- summary(object)
  lform <- formula(object)
  lts <- ltp <- NULL
  for (ly in 1:length(lreg1)) {
    lrg <- lreg1[[ly]]
    lts <- cbind(lts,c(lrg[["sigma"]],lrg[["r.squared"]],
                       lrg[["fstatistic"]]))
    ltp <- cbind(ltp,lrg[["coefficients"]][,4])
  }
  lmodel <- nrow(lts)>=5  # non-trivial model
  if (lmodel) {
    lts[4,] <- pf(lts[3,],lts[4,],lts[5,], lower.tail=FALSE)
    lts <- lts[1:4,]
  }
  dimnames(object$coefficients)[[2]] <- as.character(lform[[2]])[-1]
  dimnames(ltp) <- dimnames(object$coefficients)
  dimnames(lts) <- list(rep(c("sigma","r.squared","fstatistic","p-value"),
                            length=nrow(lts)), dimnames(ltp)[[2]])
  object$pvalues <- ltp
  object$stats <- lts
  object$sigma <- lsig <- lts["sigma",]
  lres <- residuals(object)
  if (all(lsig>0)) {
    object$stres <- sweep(lres,2,lsig,"/")
    if (length(object$weights))
      object$stres <- object$stres*sqrt(object$weights)
  }
  object$resmd <- mahalanobis(lres,0,var(lres))
  ldfr <- object$df.residual
  object$r.squared <- lr2 <- lts["r.squared",]
  object$adj.r.squared <- 1-(1-lr2)*(nrow(object$resid)-1)/ldfr
  lcomp <- c("aliased","df","cov.unscaled")
  object[lcomp] <- lreg1[[1]][lcomp]
  object$drop1 <- if (lmodel) drop1.mlm(object)
##  class(lreg) <- c("mregr","mlm","lm")
  object
} # i.mlmsum
## -----------------------------------------------------------------------
i.glm <- function(formula, data, family, fname,
                  control=NULL, vif=TRUE, termtable=TRUE, ...)
{
  ## Purpose:  internal: fit glm
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
  lfamily <- get(fname)
##-   environment(formula) <- environment()
  lcall <- match.call()
  lcall[[1]] <- as.name("glm")
  lcall$family <- lcall$fname
  lcall$fname <- lcall$control <- lcall$vif <- NULL
  lcall$x <- TRUE
  ## ---------------
  lreg <- eval(lcall, envir=environment(formula))
  ## ----------------
  lreg$leverage <- pmax(0,hat(lreg$x))
  lreg1 <- summary(lreg)
  lcoeftab <- lreg1$coef
  ly <- as.numeric(lreg$model[,1])
  ldisp <- lreg1$dispersion
  ## ---
  lfcount <- fname=="binomial"|fname=="poisson"
  lcalcdisp <- control$calcdisp
  lsuffmean <- TRUE
  if (lfcount) {
    lsuffmean <- mean(ly)>control$suffmean   # ,na.rm=TRUE
    lcd <- lcalcdisp
    if (length(lcd)==0) lcd <- lsuffmean
    if (lcd) {
      ldisp <- lreg1$deviance/lreg1$df.residual
      if (ldisp>1||length(lcalcdisp)>0) {
        lreg$distrname <- paste("quasi",fname,sep="")
        lcoeftab[,2] <- lcoeftab[,2]*sqrt(ldisp)
        lcoeftab[,3] <- lcoeftab[,3]/sqrt(ldisp)
##-         lcoeftab[,4] <- 2*pnorm(lcoeftab[,3],lower.tail=FALSE)
      }
      else ldisp <- 1
    }
  }  # else calcdisp <- FALSE
  attr(ldisp,"fixed") <- ldisp==1
  lreg$dispersion <- ldisp
  lreg$sigma <- sqrt(ldisp)
  ## ---
  if (ldisp>0) {
      lstr <- residuals(lreg,"pearson")/sqrt(ldisp)
      lnaa <- lreg$na.action
      if (class(lnaa)=="exclude") lstr <- lstr[-lnaa]
      lreg$stres <- lstr
  }
  ## bug? leverage not taken into account
  lcomp <- c("deviance","aic","df.residual","null.deviance", # "family",
    "df.null","iter","deviance.resid","aliased","df","cov.unscaled")
  lreg[lcomp] <- lreg1[lcomp]
  ## --- deviances
  ltesttype <- ifelse(ldisp==1,"Chisq","F")
  ldev <- unlist(lreg1[c("deviance", "null.deviance")])
  ldf <- lreg1$df[1:2]-c(attr(terms(lreg),"intercept"),0)
  ltbd <- cbind(deviance=c(diff(ldev),ldev), df=c(ldf,sum(ldf)),
                p.value=NA)
  dimnames(ltbd)[[1]] <- c("Model","Residual","Null")
  ltbd[1:2,3] <- pchisq(ltbd[1:2,1], ltbd[1:2,2], lower.tail=FALSE)
  if (!lsuffmean) ltbd[2,3] <- NA
  lreg$devtable <- ltbd
  ## ---
  ## cov of estimates
  ldisp <- lreg$dispersion
  if (is.null(ldisp)) ldisp <- 1
  lreg$covariance <- lcov <- lreg$cov.unscaled*ldisp
  lse <- sqrt(diag(lcov))
  lreg$correlation <- lcov/outer(lse, lse)
  ## ---
  lreg$fitfun <- "glm"
  if (termtable) {
    ltt <- i.termtable(lreg, lcoeftab, data, lcov, ltesttype, lsdy=1, vif=vif)
    lcmpn <- c("testcoef","allcoef","leverage")
    lreg[lcmpn[lcmpn%in%names(ltt)]] <- ltt
  }
  ## result of i.glm
  lreg
}
## -----------------------------------------------------------------------
i.multinomial <- function(formula, data, family, fname,
                          model=TRUE, vif=TRUE, termtable=TRUE, ...)
{
  ## Purpose:  internal: fit multinom
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
##  ltr <- control$trace
##  if (length(ltr)==0) ltr <- trace
##  require(nnet)   ## !?!
  lcall <- match.call()
  lcall[[1]] <- quote(regr0::i.multinomfit)
  lcall$fname <- lcall$family <- lcall$control <- lcall$vif <- NULL
  lcall$trace <- FALSE
  lreg <- eval(lcall, envir=environment(formula))
  ## ---------------
  if (length(lreg$na.action)) {
    lnaact <- attr(lreg$na.action,"class")
    attr(lreg$na.action,"class") <- "omit"
  } else  lnaact <- NULL ##  summary does not work with  exclude
  lreg$call$formula <- formula
  lreg1 <- summary(lreg)
  if (length(lnaact)) attr(lreg$na.action,"class") <- lnaact
  lreg$dispersion <- lreg$sigma <- 1
  lres <- lreg1$residuals
  lreg$residuals <- lres
  lcf <- lreg1$coefficients
  lreg$coefficients <- lcf
  lreg$aic <- lreg1$AIC
  ldfm <- lreg1$edf-nrow(lcf)
  lreg$df <- c(ldfm,prod(dim(lres)-1)-ldfm,ldfm)
##-   environment(lreg$call$formula) <- environment()
  lreg$fitfun <- "multinom"
  ldr1 <- if (u.debug()) drop1(lreg, test="Chisq", trace=FALSE) else 
              try(drop1(lreg, test="Chisq", trace=FALSE), silent=TRUE)
  if (class(ldr1)[1]=="try-error") {
    warning(paste(":regr/i.multinom: drop1 did not work.",
                  "I return the multinom object"))
    class(lreg) <- c("orig",class(lreg))
    return(lreg)
  } else {
  ldr1 <- ldr1[-1,]}
  names(ldr1) <- c("df", "AIC", "Chisq","p.value") #if(calcdisp) "F" else
  lreg$testcoef <- lreg$drop1 <- ldr1
## result of i.multinomial
  lreg
}
## -----------------------------------------------------------------------
i.polr <- function(formula, data, family, fname, weights = NULL, 
                   model=TRUE, vif=TRUE, termtable=TRUE, ...)
{
  ## Purpose:  internal: fit ordered y
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
##  require(MASS)   ## !?!
  lcall <- match.call()
  lcall[[1]] <- as.name("polr") ## quote(regr0::i.polrfit)
  lcall$fname <- lcall$control <- lcall$family <- lcall$vif <- NULL
  lcall$Hess <- TRUE
  lreg <- eval(lcall, envir=environment(formula))
##  lreg$call$formula <- formula
  lreg$w <- weights
  lreg$leverage <- hat(lreg[["x"]])
  lreg1 <- if (u.debug()) summary(lreg) else
           try(summary(lreg))
  if (class(lreg1)[1]=="try-error") {
    warning(paste(":regr/i.polr: summary did not work.",
                  "I return the polr object"))
##    lreg$call$data <- call$data
    class(lreg) <- c("orig","polr")
    return(lreg)
  }   ## ---
  ## model.matrix
##-   ldata <- eval(data, envir=environment(formula))
##  lreg$x <- model.matrix(formula, ldata)
  lcf <- lreg1$coefficients
  lreg$intercepts <- lcf[(lreg1$pc+1):nrow(lcf),1:2]
  lreg$stres <- NULL
  ## cov of estimates!
  lreg$covariance <- lcov <- vcov(lreg)
  lse <- sqrt(diag(lcov))
  lreg$correlation <- lcov/outer(lse, lse)
  ## --- deviances
  lreg$fitfun <- "polr"
  if (termtable) {
    ltt <- i.termtable(lreg, lreg1$coef, data, lcov, ltesttype="Chisq",
                       lsdy=1, vif=vif, leverage=TRUE)
    lcmpn <- c("testcoef","allcoef","leverage")
    lreg[lcmpn[lcmpn%in%names(ltt)]] <- ltt
  }
  lreg$dispersion <- 1
  ## result of i.polr
  lreg
}
## -----------------------------------------------------------------------
i.survreg <-
  function(formula, data, family, fname="ph", method, control,
           vif=TRUE, termtable=TRUE, ...)
{
  ## Purpose:  internal: fit ordered y
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
##  require(survival)  ## !?!
  lcall <- match.call()
  ## b. --- method
  if (fname=="ph") {
    lfitfun <- "coxph"
    lcall[[1]] <- quote(survival::coxph)
  } else {
    lfitfun <- "survreg"
    lcall[[1]] <- quote(survival::survreg)
    lcall$dist <- fname
    lcall$method <- lcall$control <- NULL
  }
  lcall$fname <- lcall$family <- lcall$vif <- NULL
  lcall$x <- TRUE
  ## ---
  lreg <- eval(lcall, envir=environment(formula))
  ## ---
##  lreg$call$formula <- formula
  lreg1 <- if (u.debug()) summary(lreg) else
           try(summary(lreg), silent=TRUE)
  if (class(lreg1)[1]=="try-error") {
    warning(paste(":regr/i.survreg: summary did not work.",
                  "I return the survreg object"))
##    lreg$call$data <- call$data
    class(lreg) <- c("orig",class(lreg))
    return(lreg)
  }   ## ---
  lcf <- lreg1$coefficients
  lreg$stres <- NULL
  ## --- deviances
  ## lreg$scale
  if (lfitfun=="survreg") {
    attr(lreg$scale,"fixed") <- length(lcall$scale)>0
  }
  if (lfitfun=="coxph") {
    lreg1$table <- lreg1$coefficients
    lreg$df.residual <- length(lreg$residual)-length(lreg$coefficients)
  }
  lreg$aic <- extractAIC(lreg)
  lreg$deviance <- -2*lreg$loglik
  lchi <- 2*diff(lreg1$loglik)
  ldf <- sum(lreg$df) - lreg$idf
  ltbd <- cbind(deviance=c(lchi,-2*lreg1$loglik[2]),
                df=c(ldf, lreg$df.residual+lreg$df),
                p.value=c(pchisq(lchi,ldf,lower.tail=FALSE),NA))
  dimnames(ltbd)[[1]] <- c("Model","Null")
  lreg$devtable <- ltbd
  lcov <- lreg$var
  lsd <- sqrt(diag(lcov))
  lreg$correlation <- lcov/outer(lsd,lsd)
  lreg$fitfun <- lfitfun
  if (termtable) {
    ltt <- i.termtable(lreg, lreg1$table, data, lcov, ltesttype="Chisq",
                       lsdy=1, vif=vif)
    ## log(scale): signif<-NA
    lcmpn <- c("testcoef","allcoef","leverage")
    lreg[lcmpn[lcmpn%in%names(ltt)]] <- ltt
  }
  lreg$distrname <- if (lfitfun=="coxph") "prop.hazard" else lreg$dist
  ## result of i.survreg
  lreg
}

## -----------------------------------------------------------------------
Tobit <- function(data, limit=0, transform=NULL, log=FALSE, ...)
{
  ## Purpose:   create a Surv object for tobit regression
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  1 Jan 2010, 21:49
##  require(survival)  ## !?!
  ltrs <- as.character(substitute(transform))
  data <- pmax(data,limit)
  if (log) {
      transform <- logst
      ltrs <- "logst"
  }
  if (!is.null(transform)) {
    if (is.character(transform)) transform <- get(transform)
    if (!is.function(transform))
      stop("!Tobit! argument 'transform' does not yield a function")
    ldt <- transform(c(limit,data), ...)
    data <- ldt[-1]
    limit <- ldt[1]
  }
  if (sum(data==limit,na.rm=TRUE)<=1)
    warning(":Tobit: <= 1 observation equal to `limit`")
  rr <- Surv(data, data!=limit, type="left")
  structure(rr, distribution="gaussian", transform=ltrs, limit=limit,
            class=c(class(rr), "matrix"))
}

## -----------------------------------------------------------------------
i.termtable <- function(lreg, lcoeftab, ldata, lcov, ltesttype="F",
                        lsdy, vif=TRUE, leverage=vif)
{
  ## Purpose:  generate term table for various models
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 15:37
  if(length(attr(terms(lreg),"term.labels"))==0)
    return(list(test=data.frame(coef=c(lreg$coef,NA)[1],stcoef=NA,signif=NA,
                  R2.x=NA,df=length(lreg$coef),p.value=NA)))
## drop1
  ldr1 <-
      if (class(lreg)[1]%in%c("lm","lmrob")) {
          if (u.debug()) 
              drop1Wald(lreg, test=ltesttype, scope=terms(lreg)) else
              try(drop1Wald(lreg, test=ltesttype, scope=terms(lreg)),
              silent=TRUE) } else {
          if (u.debug()) 
              drop1(lreg, test=ltesttype, scope=terms(lreg)) else
              try(drop1(lreg, test=ltesttype, scope=terms(lreg)),
                  silent=TRUE)
          }
  if (class(ldr1)[1]=="try-error") {
    warning(paste(":regr: drop1 did not work. I return the table produced by ",
                  lreg$fitfun))
##-     lsum <- summary(lreg)
##-     lcft <- lsum$coef
##-     if (length(lcft)==0) lcft <- lsum$parameters ## nls
##-     return(list(test=lcft)) # !!! noch reparieren
    return(list(testcoef=lcoeftab))
  }
  ldr1 <- ldr1[-1,]
  ldr1$RSS <- NULL # same ncol for lm and glm
  if (inherits(lreg,"rlm"))  ldr1[,4] <- ldr1[,2]/ldr1[,1] ## !!!
  if (inherits(lreg,"mlm")||inherits(lreg,"manova"))
    return(list(testcoef=ldr1))  ## !!! needs much more
## degrees of freedom
  ldfr <- lreg$df.residual
  if (length(ldfr)==0) ldfr <- lreg$df[2]
  ltstq <- if (ltesttype=="F") qf(0.95,ldr1[,1],ldfr) else {
    if (ltesttype=="Chisq") qchisq(0.95,ldr1[,1]) else NA }
## coefficients
  lcoef <- lreg$coefficients
##-   ldt <- lreg$model
##-   if (is.null(ldt)) ldt <- ldata
##  lmmt <- model.matrix(lreg$call$formula, data=ldt) ## or lreg$model ## needed for vif & la
  lmmt <- lreg[["x"]]
  if (length(lmmt)==0)
      lmmt <- model.matrix(lreg)
  lasg <- attr(lmmt,"assign")
  if (class(lreg)[1]%in%c("polr","coxph")) lasg <- lasg[-1]
## terms with 1 coef
##  lcont1 <- lcont <- lasg[!lasg%in%lasg[duplicated(lasg)]]
  lcont1 <- lcont <- lasg[!lasg%in%lasg[duplicated(lasg)]]
##  which(table(lasg)==1)-1
## vif --> R2.x
  lr2 <- NA
  if (vif) {
    lvift <-     ## lterms: n of levels for each term
        if (u.debug()) vif.regr(lreg, lcov, lmmt) else
        try(vif.regr(lreg, lcov, lmmt), silent=TRUE)
    if (class(lvift)[1]=="try-error" || length(lvift)==0) {
      warning(":regr/i.termtable: error in the calculation of R2.xs")
      lvif <- NA
    } else lvif <- lvift[,3]^2
    lr2 <- 1-1/lvif
  }
## prepare table
  lpvcol <- pmatch("Pr(",names(ldr1), nomatch=ncol(ldr1))
  ltb <- data.frame(coef=NA, stcoef=NA, ci25=NA, ci975=NA, signif=NA, R2.x=lr2,
                    df=ldr1[,1], p.value=ldr1[,lpvcol], testst=ldr1[,lpvcol-1])
  row.names(ltb) <- row.names(ldr1)
## intercept
  if ("(Intercept)"==names(lcoef)[1]) {
    ltstint <- # if(class(lreg)[1]%in%c("lm","nls","rlm"))
      lcoeftab[1,3]^2 # else lcoeftab[1,3]
    ltb <- rbind("(Intercept)"=c(NA,NA,NA,NA,NA,NA,1,NA,ltstint),ltb)
    lcont1 <- lcont+1  # row number in dr1
    ltstq <- c(qf(0.95,1,ldfr), ltstq)
  }
## signif
  ltb$signif <- lsg <- sqrt(ltb$testst/ltstq)
##-   lprcol <- pmatch("Pr(",names(ldr1))
##-   ldr1t <-  if (is.na(lprcol)) NA else ldr1[,lprcol]
##      -qnorm(ldr1[,lprcol]/2)/qnorm(0.975)
##-   lnt <- nrow(ldr1)
## significance
##-   ldr1t <- sqrt(pmax(0,ldr1[,4])/ltstq)
## coefficients for terms with 1 df
  if (length(lcont)) {
    ltlb <- dimnames(ltb)[[1]]
    lclb <- ltlb[lcont1]
    licf <- pmatch(lclb,dimnames(lcoeftab)[[1]])
    if (any(is.na(licf))) warning(":regr: bug: coefficients not found")
    ## !!! bug if some coefs have 0 df
    ljc <- match(lcont,lasg) # index of coefs for cont variables
    lcf <- lcoef[ljc]
    ## standardized coefficients
    lstcf <- lcf[lcont>0] *  # exclude intercept term
      sqrt(apply(lmmt[,names(lcf[lcont>0]),drop=FALSE],2,var)) / lsdy
    ## fill in
    ltb$coef[lcont1] <- lcf
    lci <- lcf*(1+outer(1/lsg[lcont1], c(-1,1)))
       ## confint(lreg,row.names(ltb)[lcont1]) does not always work...
    ltb[lcont1,c("ci25","ci975")] <- lci
    ltb$stcoef[lcont1[lcont>0]] <- lstcf
    ltb[lcont1,"signif"] <- sign(lcf)*ltb[lcont1,"signif"]
}
##-   if (row.names(lcoeftab)[nrow(lcoeftab)]=="Log(scale)") { # survreg
##-     ltsc <- lcoeftab[nrow(lcoeftab),]
##-     ltb <- rbind(ltb,"log(scale)"=c(ltsc[1],NA,NA,NA,ltsc[3]/qnorm(0.975),
##-                        NA,1,NA,NA))
##-  }
## --- dummy coef
  lrg <- lreg
  class(lrg) <- "lm"
  if (inherits(lreg,"polr")) lrg$coefficients <- c("(Intercept)" = NA, lcoef)
  lallcf <- dummy.coef.regr(lrg) 
##                try(dummy.coef(lrg), silent=TRUE)
##-   if (class(lallcf)=="try-error") {
##-     warning("dummy.coef did not work")
##-     lallcf <- NULL
##-   }
  rr <- list(testcoef=ltb, allcoef=lallcf)
  if (leverage) rr <- c(rr,list(leverage=hat(lmmt)))
  rr
}
## ==========================================================================
print.regr <- function (x, call=TRUE, correlation = FALSE,
    dummy.coef = getUserOption("show.dummy.coef"),
    testcoefcol = getUserOption("regr.testcoefcol"),
    digits = max(3, getUserOption("digits")-2), 
    symbolic.cor = p > 4, signif.stars = getOption("show.signif.stars"),
    residuals=FALSE, niterations=FALSE, ...)
{
##
  ## doc
  ldoc <- getUserOption("doc")
  if (length(ldoc)==0) ldoc <- 1
  if (ldoc>=1) if (length(tit(x)))
    cat("\n ",tit(x),"\n")
  if (ldoc>=2) if (length(doc(x)))
    cat(" ",paste(doc(x),"\n "))
  ## mlm
  if (inherits(x,"mlm")) return(invisible(print.mregr(x, ...)))
  ## preparation
##-   if (length(dummycoef)==0)
##-     dummycoef <- c(getUserOption("show.dummy.coef"),TRUE)[1]
  ## call, fitting fn, residuals
  if (call) {
    if(!is.null(x$call)) {
      cat("\nCall:\n")
      cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),"\n", sep = "")
    }
    cat("Fitting function: ",x$fitfun,"\n")
  }
  df <- x$df
    rdf <- c(x$df.resid,df[2])[1]
  if (residuals) {
    resid <- x$residuals
##-     cat(if (!is.null(x$w) && diff(range(x$w)))
##-         "Weighted ", "Residuals:\n", sep = "")
    cat("Residuals:\n")
    if (rdf > 5) {
      nam <- c("Min", "1Q", "Median", "3Q", "Max")
      rq <- if (length(dim(resid)) == 2)
        structure(apply(t(resid), 1, quantile),
                  dimnames = list(nam, dimnames(resid)[[2]]))
        else structure(quantile(resid), names = nam)
      print(rq, digits = digits, ...)
    } else {
      if (rdf > 0) print(resid, digits = digits, ...)
      else  cat("ALL", df[1],
                "residuals are 0: no residual degrees of freedom!\n")
    }
  }
  ## coefficients
    nsingular <- df[3] - df[1]
    if ((!is.na(nsingular))&&nsingular>0)
        cat("\nCoefficients: (", nsingular,
            " not defined because of singularities)\n", sep = "")
    # else {
    if (!is.null(x$sigma))
      if((!is.finite(x$sigma))||x$sigma<=0)
        cat("\n!!! Error variance is 0 !!!")
  ## coef table
  ltc <- x$testcoef
  if (length(ltc)>0) {
    lltc <- TRUE
    if(!is.null(testcoefcol)) {
      if (all(testcoefcol=="")) lltc <- FALSE else {
        ljp <- match(testcoefcol,colnames(ltc), nomatch=0)
        if (sum(ljp)==0)
          warning(":print.regr: no valid columns of  testcoef  selected") else
        ltc <- ltc[,ljp,drop=FALSE]
      }
    }
    if (lltc) {
      cat("\nTerms:\n")
      ljrp <- pmatch(c("R2","p."),colnames(ltc))
      lip <- nna(ljrp)
      if (length(lip)) ltc[,lip] <- round(as.matrix(ltc[,lip]),max(3,digits))
      ltcf <- format(ltc)
      lsigst <- signif.stars*( signif.stars && !is.na(ljrp[2]) &&
                              any((pv <- ltc[,ljrp[2]]) < 0.1, na.rm=TRUE) )
      if (lsigst) {
        lsignif <- symnum(pv, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                          symbols = c("***", "**", "*", ".", " "))
        ltcf <- cbind(ltcf,format(lsignif))
        names(ltcf)[ncol(ltcf)] <- " "
    }
      print(ltcf, quote=FALSE)
      if (lsigst > 1) 
        cat("---\nSignif. codes:  ", attr(lsignif, "legend"),"\n", sep = "")
    } ## end if(lltc)
  ## --- error block
  } else {
    if (length(x$coef)) {
      cat("\nCoefficients:\n")
      print(x$coef)
    }
  }
##-   if (length(x$binlevels)>0) {
##-     cat("\nFactor(s) with two levels converted to 0-1 variable(s):\n")
##-     print(as.matrix(data.frame(x$binlevels,row.names=0:1)))
##-   }
##    cat("\n")
  ## special for polr
  if (length(x$intercepts)) { # polr
    cat("Intercepts:\n")
    print(x$intercepts)
  }
  ## error
  if (length(x$sigma) && !u.true(attr(x$sigma,"fixed")))
    cat("\nSt.dev.error: ", formatC(x$sigma, digits = digits),
        "  on", rdf, "degrees of freedom\n")
  if (length(x$r.squared)&&!is.na(x$r.squared))
    cat("Multiple R^2: ", formatC(x$r.squared, digits = digits),
        "   Adjusted R-squared:",
        formatC(x$adj.r.squared, digits = digits),"\n"
        )
  if (length(x$fstatistic)>0) {
    cat("F-statistic:  ", formatC(x$fstatistic[1],
            digits = digits), "  on", x$fstatistic[2], "and", x$fstatistic[3],
            "d.f.,  p.value:", formatC(pf(x$fstatistic[1],
            x$fstatistic[2], x$fstatistic[3], lower.tail=FALSE),
                                       digits = digits),
         "\n")
    }
  ## deviances
  if (length(x$deviance)>0) {
    if (length(x$devtable)>0) {
  #    cat("\n")
      print(x$devtable)
    }
    cat("\nDistribution: ",x$distrname)
    if (length(x$dispersion))
      cat(".  Dispersion parameter: ",
          if ((!is.null(attr(x$dispersion,"fixed")))&&
              attr(x$dispersion,"fixed"))
             "fixed at ", format(x$dispersion))
    else if (length(x$scale))
      cat(".  Shape parameter (`scale`): ",
          if ((!is.null(attr(x$scale,"fixed")))&&
              attr(x$scale,"fixed"))
             "fixed at ", format(x$scale))
    cat("\nAIC: ", format(x$aic, digits = max(4, digits + 1)), "\n", sep = "")
    if (niterations&&length(x$iter)>0)
      cat("Number of iterations:", x$iter, "\n")
  }
  ## --- additional coefficients
  if (x$distrname=="multinomial") {
    cat("\nCoefficients:\n")
    print(t(x$coefficients))
  } else {
    if (length(ltc)&u.true(dummy.coef)) {        
      lidf <- match("df",colnames(x$testcoef))
      if (is.na(lidf)) {
        if (getOption("verbose"))
          warning(":print.regr: df of coef not available")
      } else { ## dummy coefficients
        mterms <-
          unique(c(row.names(x$testcoef)[x$testcoef[,"df"]>1],
                   names(attr(x$terms,"dataClasses")[-1]%in%
                         c("factor","ordered")) ))
        if (length(mterms)>0 & length(x$allcoef)>0) {
          imt <- mterms%in%names(x$allcoef)
          mt <- x$allcoef[mterms[imt]]
          if (length(mt)>0) {
            cat("\nCoefficients for factors:\n")
            print(unclass(mt)) } ## avoid  print.dummy_coef 
        } ## else  cat("\n")
      }}
  }
  ## ---- correlation
  correl <- x$correlation
  if (length(correl)>0 && correlation) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (symbolic.cor) {
        symbc <- symnum(correl, symbols=c(" ", ".", ",", "+", "*", "H"))
        symbl <- attr(symbc,"legend")
        attr(symbc,"legend") <- NULL
        print(symbc)
        cat("\nSymbols:  ",symbl,"\n")
      } else {
        correl[!lower.tri(correl)] <- NA
        print(correl[-1, -p, drop = FALSE], digits = digits, na = "")
      }
    }
  }
##  cat("\n")
  invisible(x)
}
## ==========================================================================
summary.regr <- function(object, ...)  object ## dispersion=NULL,
## ==========================================================================
dummy.coef.regr <- function (object, use.na = FALSE, ...) 
{
    if (length(object$xlevels)==0)  return(as.list(coef(object)))
    Terms <- terms(object)
    tl <- attr(Terms, "term.labels")
    int <- attr(Terms, "intercept")
    facs <- attr(Terms, "factors")[-1, , drop = FALSE]
    Terms <- delete.response(Terms)
    mf <- object$model
    if (is.null(mf)) mf <- model.frame(object)
    xtnm <- dimnames(facs)[[1]]  ## names
    xtlv <- lapply(mf[,xtnm, drop=FALSE],function(x) levels(x)) ## levels
    xtnl <- pmax(sapply(xtlv,length),1)  ## number of levels
    termnl <- apply(facs, 2L, function(x) prod(xtnl[x > 0]))
    nl <- sum(termnl)
    ## df.dummy: data frame of vars
    args <- setNames(vector("list", length(xtnm)), xtnm)
    for (i in xtnm)
        args[[i]] <- if (xtnl[[i]] == 1)  rep.int(1, nl)    else
          factor(rep.int(xtlv[[i]][1L], nl), levels = xtlv[[i]])
    df.dummy <- as.data.frame(args) # do.call("data.frame", args)
    names(df.dummy) <- xtnm
    ## rnn: names of rows
    pos <- 0
    rn <- rep.int(tl, termnl)
    rnn <- rep.int("", nl)
    ## fill df.dummy
    for (j in tl) {
        i <- unlist(xtnm[facs[, j] > 0])
        ifac <- i[xtnl[i] > 1]
        if (length(ifac) == 0L) {
            rnn[pos + 1] <- j
        }
        else if (length(ifac) == 1L) {
            df.dummy[pos + 1L:termnl[j], ifac] <- xtlv[[ifac]]
            rnn[pos + 1L:termnl[j]] <- as.character(xtlv[[ifac]])
        }
        else {
            tmp <- expand.grid(xtlv[ifac])
            df.dummy[pos + 1L:termnl[j], ifac] <- tmp
            rnn[pos + 1L:termnl[j]] <- apply(as.matrix(tmp), 
                1L, function(x) paste(x, collapse = ":"))
        }
        pos <- pos + termnl[j]
    }
    attr(df.dummy,"terms") <- attr(mf,"terms")
    lcontr <- object$contrasts
    lci <- sapply(df.dummy,is.factor)
    lcontr <- lcontr[names(lci)[lci]] ## factors with 1 level have disappeared (?) 
    mm <- model.matrix(Terms, df.dummy, lcontr)
    if (any(is.na(mm))) {
        warning("some terms will have NAs due to the limits of the method")
        mm[is.na(mm)] <- NA
    }
    coef <- object$coefficients
    if (!use.na) 
        coef[is.na(coef)] <- 0
    asgn <- attr(mm, "assign")
    res <- setNames(vector("list", length(tl)), tl)
    for (j in seq_along(tl)) {
        keep <- asgn == j
        ij <- rn == tl[j]
        res[[j]] <- setNames(drop(mm[ij, keep, drop = FALSE] %*% 
            coef[keep]), rnn[ij])
    }
    if (int > 0) {
        res <- c(list(`(Intercept)` = coef[int]), res)
    }
    class(res) <- "dummy_coef"
    res
}
## ====================================================================
confint.regr <- function(fitted, ...)
{
if (!inherits(fitted, c("glm","nls"))) {
    class(fitted) <- class(fitted)[-1]
    return(confint(fitted, ...))
}
if (inherits(fitted, "glm")) {
      ## confint needs $coefficients from object (a vector) as well as
      ## from its sumary (a matrix containing 'Std. Error"
  summary <- function(fitted) 
    list(coefficients = cbind(fitted$coefficients,
                                 "Std. Error"=sqrt(diag(fitted$covariance))) )
  class(fitted) <- class(fitted)[-1]
} else { ## workaround: call nls again, since  profile.nls  is difficult to adapt...
  call <- fitted$call
  call$start <- fitted$coefficients
  call$nonlinear <- NULL
  call[[1]] <- as.name("nls")
  fitted <- eval(call, parent.frame())
}
  confint(fitted, ...)
}
## ==========================================================================
drop1.regr <-
  function (object, scope=NULL, scale = 0, test = NULL, k = 2,
            sorted = FALSE, add=FALSE, ...)
{
  ## Purpose:    drop1/add1 for regr objects
  ## ----------------------------------------------------------------------
  lfam <- object$distrname
  lres <- object$residuals
  if (is.null(test)) test <- if (is.null(lfam)) "none" else {
    if ((lfam=="gaussian"&&as.character(object$fitfun)%in%c("lm","roblm"))|
        ((lfam=="binomial"|lfam=="poisson")&&object$dispersion>1)) {
          if (inherits(object,"mlm")) "Wilks" else "F" }
    else "Chisq"
  }
  if (length(scope)==0) {
    scope <- if (add) terms2order(object) else drop.scope(object)
  } else
##-     if (!is.character(scope))
##-             scope <- attr(terms(update.formula(object, scope)),
##-                 "term.labels")
    scope
  if (length(scope)==0) {
    warning(":drop1/add1.regr! no valid scope")
    ldr1 <- data.frame(Df = NA, "Sum of Sq" = NA, RSS =NA, AIC = NA,
                      row.names = "<none>")
    return(ldr1)
  }
  class(object) <- setdiff(class(object), "regr")
  dr1 <- if (add) {
    if (class(object)[1]=="lmrob")
        stop("!add1.regr! 'add1' not (yet) available for 'lmrob' objects")
    ldata <- eval(object$call$data)
    li <- row.names(ldata)%in%RNAMES(object$residuals)
##-     if (any(!li)) {
##-       ldata <- ldata[li,]
##-     }
    if (is.null(ldata[li,])) stop("!step.regr! no data found ")
    lvars <-unique(c(all.vars(formula(object)),
                     if (is.formula(scope)) all.vars(scope) else scope))
    lvars <- lvars[lvars%in%names(ldata)]
    linna <- li & !apply(is.na(ldata[,lvars]),1,any)
    lnobs <- sum(linna)
##-     ldt <- ldata[!lina,]
##-     ldt <- na.omit(ldata[,lvars,drop=FALSE])
    lnrd <- sum(li) # NROW(object$residuals)
    lfc <- object$funcall
    if (lnobs!= lnrd) {
      warning(":add1.regr: refitting object to ",lnobs," / ",lnrd,
              " observations due to missing values")
      lfc <- object$funcall
      if(!is.null(lsubs <- eval(lfc$subset))) {
        lnsubs <- rep(FALSE,length(linna))
        lnsubs[lsubs] <- TRUE
        linna <- linna &!lnsubs
      }
      lfc$subset <- linna
      object <- eval(lfc)
##-       object$call[[1]] <-
##-         if (is.null(lfc)) as.name(class(object)[1]) else
##-            lfc[[1]]
##-       object <- update(object, subset=linna)
#      environment(object$call$formula) <- environment()
##-       class(object) <- setdiff(class(object), "regr")
    }
    if (!all(linna)) { ## needed if NA's have been generated by transformations
      lfc$subset <- linna
      object <- eval(lfc)
    }
    add1(object, scope=scope, scale=scale, test=test, k=k, ...)
  } else {
    if (class(object)[1]%in%c("lmrob")) ## to be expanded
       drop1Wald(object, test="F", ...) else {
    ldata <- object$allvars # eval(object$call$data)
    if (is.null(ldata)) stop("!drop1.regr! no data found ")
    lina <- apply(is.na(ldata),1,any)
    if (any(lina)) ldata[lina,] <- NA
    object$call$data <- ldata
    drop1(object, scope=scope, scale=scale, test=test, k=k, ...)
  }}
##-   rnm <- row.names(dr1)
##-   row.names(dr1) <- paste(ifelse(substring(rnm,1,1)=="<","",
##-                                  if (add) "+ " else "- "),rnm,sep="")
  attr(dr1,"drop") <- !add
##-   if(add) attr(dr1,"ndropped") <- lndiff
  if (sorted) {
    lsrt <- nna(match(c("AIC","p.value"),colnames(dr1)))
    if (length(lsrt)) dr1 <- dr1[order(dr1[, lsrt[1]]), ]
  }
  dr1
}
## ==========================================================================
add1.regr <-
  function (object, scope=NULL, scale = 0, test = NULL, k = 2,
            sorted = FALSE, ...)
{
  ## Purpose:    add1 for regr objects
  ## ----------------------------------------------------------------------
  if (!is.null(scope)) {
    if (is.character(scope)) scope <- paste(scope,collapse="+")
    if (is.formula(scope)) scope <- last(as.character(scope))
    scope <- as.formula(paste("~ ",formula(object)[3],"+",scope))
  }
  drop1.regr(object, scope=scope, scale=scale, test=test, k=k,
             sorted=sorted, add=TRUE, ...)
}
## ==========================================================================
drop1Wald <-
  function (object, scope=NULL, scale = 0, test = c("none", "Chisq", "F"),
            k = 2) 
{
    x <- model.matrix(object)
    offset <- model.offset(model.frame(object))
    n <- nrow(x)
    asgn <- attr(x, "assign")
    tl <- attr(object$terms, "term.labels")
    if (is.null(scope)) 
        scope <- drop.scope(object)
    else {
        if (!is.character(scope)) 
            scope <- attr(terms(update.formula(object, scope)), 
                "term.labels")
        if (!all(match(scope, tl, 0L) > 0L)) 
            stop("scope is not a subset of term labels")
    }
    ndrop <- match(scope, tl)
    ns <- length(scope)
    rdf <- object$df.residual
    chisq <- object$sigma^2 * rdf
    ## sum(weighted.residuals(object)^2, na.rm = TRUE)
    ## deviance.lm(object)
    dfs <- numeric(ns)
    RSS <- numeric(ns)
    cov <- object$cov.unscaled
    if (is.null(cov)) cov <- object$covariance/object$sigma^2
    if (is.null(cov)) stop("!drop1Wald! no covariance matrix found")
    cf <- object$coefficients
##-     y <- object$residuals + predict(object)
    for (i in 1:ns) {
        ii <- seq_along(asgn)[asgn == ndrop[i]]
        RSS[i] <- if (length(ii)==1) cf[ii]^2/cov[ii,ii] else
          cf[ii]%*%solve(cov[ii,ii])%*%cf[ii]
        dfs[i] <- length(ii)
##-         if (all.cols) 
##-             jj <- setdiff(seq(ncol(x)), ii)
##-         else jj <- setdiff(na.coef, ii)
##-         z <- if (iswt) 
##-             lm.wfit(x[, jj, drop = FALSE], y, wt, offset = offset)
##-         else lm.fit(x[, jj, drop = FALSE], y, offset = offset)
##-         dfs[i] <- z$rank
##-         oldClass(z) <- "lm"
##-         RSS[i] <- deviance(z)
    }
    scope <- c("<none>", scope)
    dfs <- c(object$rank, dfs)
    RSS <- chisq + c(0, RSS)
    if (scale > 0) 
        aic <- RSS/scale - n + k * dfs
    else aic <- n * log(RSS/n) + k * dfs
##-     dfs <- dfs[1] - dfs
##-     dfs[1] <- NA
    aod <- data.frame(Df = dfs, "Sum of Sq" = c(NA, RSS[-1] - 
        RSS[1]), RSS = RSS, AIC = aic, row.names = scope, check.names = FALSE)
    if (scale > 0) 
        names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
    test <- match.arg(test)
    if (test == "Chisq") {
        dev <- aod$"Sum of Sq"
        if (scale == 0) {
            dev <- n * log(RSS/n)
            dev <- dev - dev[1]
            dev[1] <- NA
        }
        else dev <- dev/scale
        df <- aod$Df
        nas <- !is.na(df)
        dev[nas] <- pchisq(dev[nas], df[nas], lower.tail = FALSE)
        aod[, "Pr(Chi)"] <- dev
    }
    else if (test == "F") {
        dev <- aod$"Sum of Sq"
        dfs <- aod$Df
        rdf <- object$df.residual
        rms <- aod$RSS[1]/rdf
        Fs <- (dev/dfs)/rms
        Fs[dfs < 1e-04] <- NA
        P <- Fs
        nas <- !is.na(Fs)
        P[nas] <- pf(Fs[nas], dfs[nas], rdf, lower.tail = FALSE)
        aod[, c("F value", "Pr(F)")] <- list(Fs, P)
    }
    head <- c("Single term deletions (Wald test)", "\nModel:",
              deparse(as.vector(formula(object))), 
        if (scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

## ==========================================================================
step <- function(object, ...)
  UseMethod("step")
step.default <- get("step", pos="package:stats")

step.regr <- function (object, scope, scale = 0,
  direction = c("both", "backward","forward"), trace = 1, keep = NULL,
  steps = 1000, k = 2, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 15 May 2012, 07:58
    mydeviance <- function(x, ...) {
        dev <- deviance(x)
        if (!is.null(dev))
            dev
        else extractAIC(x, k = 0)[2L]
    }
    cut.string <- function(string) {
        if (length(string) > 1L)
            string[-1L] <- paste0("\n", string[-1L])
        string
    }
    re.arrange <- function(keep) {
        namr <- names(k1 <- keep[[1L]])
        namc <- names(keep)
        nc <- length(keep)
        nr <- length(k1)
        array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr,
            namc))
    }
    step.results <- function(models, fit, object, usingCp = FALSE) {
        change <- sapply(models, "[[", "change")
        rd <- sapply(models, "[[", "deviance")
        if(is.matrix(rd)) rd <- rd[2,]
        dd <- c(NA, abs(diff(rd)))
        rdf <- sapply(models, "[[", "df.resid")
        ddf <- c(NA, diff(rdf))
        AIC <- sapply(models, "[[", "AIC")
        heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
            "\nInitial Model:", deparse(formula(object)), "\nFinal Model:",
            deparse(formula(fit)), "\n")
        aod <- data.frame(Step = I(change), Df = ddf, Deviance = dd,
            `Resid. Df` = rdf, `Resid. Dev` = rd, AIC = AIC,
            check.names = FALSE)
        if (usingCp) {
            cn <- colnames(aod)
            cn[cn == "AIC"] <- "Cp"
            colnames(aod) <- cn
        }
        attr(aod, "heading") <- heading
        fit$anova <- aod
        fit
      }
    ## end step.results
    Terms <- terms(object)
    object$call$formula <- object$formula <- Terms
##-     ## !!! begin mod by WSt for regr objects
##-     ldata <- eval(object$call$data)
##-     if (is.null(ldata)) stop("!step.regr! no data found ")
##-     lvars <- all.vars(object$formula)
##-     lina <- apply(is.na(ldata[,lvars]),1,sum)>0
##-     lnna <- sum(lina)
##-     ldt <- ldata[!lina,]
##-     lcdata <- object$call$data
##-     object$call$data <- ldt
##-     ## !!! end
    md <- missing(direction)
    direction <- match.arg(direction)
    backward <- direction == "both" | direction == "backward"
    forward <- direction == "both" | direction == "forward"
    if (missing(scope)) {
        fdrop <- numeric()
        fadd <- attr(Terms, "factors")
        if (md)
            forward <- FALSE
    }
    else {
        if (is.list(scope)) {
            fdrop <- if (!is.null(fdrop <- scope$lower))
                attr(terms(update.formula(object, fdrop)), "factors")
            else numeric()
            fadd <- if (!is.null(fadd <- scope$upper))
                attr(terms(update.formula(object, fadd)), "factors")
        }
        else {
            fadd <- if (!is.null(fadd <- scope))
                attr(terms(update.formula(object, scope)), "factors")
            fdrop <- numeric()
        }
    }
    models <- vector("list", steps)
    if (!is.null(keep))
        keep.list <- vector("list", steps)
    n <- nobs(object, use.fallback = TRUE)
    fit <- object
    bAIC <- extractAIC(fit, scale, k = k, ...)
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]
    if (is.na(bAIC))
        stop("AIC is not defined for this model, so `step` cannot proceed")
    nm <- 1
    if (trace) {
        cat("Start:  AIC=", format(round(bAIC, 2)), "\n", cut.string(deparse(formula(fit))),
            "\n\n", sep = "")
        utils::flush.console()
    }
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n -
        edf, change = "", AIC = bAIC)
    if (!is.null(keep))
        keep.list[[nm]] <- keep(fit, bAIC)
    usingCp <- FALSE
    while (steps > 0) {
        steps <- steps - 1
        AIC <- bAIC
        ffac <- attr(Terms, "factors")
        scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
        aod <- NULL
        change <- NULL
        if (backward && length(scope$drop)) {
            aod <- drop1(fit, scope$drop, scale = scale, trace = trace,
                k = k, ...)
            rn <- row.names(aod)
            row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep = " "))
            if (any(aod$Df == 0, na.rm = TRUE)) {
                zdf <- aod$Df == 0 & !is.na(aod$Df)
                change <- rev(rownames(aod)[zdf])[1L]
            }
        }
        ## forward
        if (is.null(change)) {
          if (forward && length(scope$add)) {
            aodf <- add1(fit, scope$add, scale = scale, trace = trace,
                         k = k, ...)
            rn <- row.names(aodf)
            row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], sep = " "))
            aod <- if (is.null(aod)) aodf
            else {
              names(aodf) <- names(aod)
              rbind(aod, aodf[-1, , drop = FALSE])
            }
            }
          attr(aod, "heading") <- NULL
          nzdf <- if (!is.null(aod$Df))
            aod$Df != 0 | is.na(aod$Df)
          aod <- aod[nzdf, ]
          if (is.null(aod) || ncol(aod) == 0)
            break
          nc <- match(c("Cp", "AIC"), names(aod))
          nc <- nc[!is.na(nc)][1L]
          o <- order(aod[, nc])
          if (trace)
            print(aod[o, ])
          if (o[1L] == 1)
            break
          change <- rownames(aod)[o[1L]]
        }
        ## update
        usingCp <- match("Cp", names(aod), 0L) > 0L
        fit <- update(fit, paste("~ .", change), evaluate = FALSE)
        fit <- eval.parent(fit)
        nnew <- nobs(fit, use.fallback = TRUE)
        if (all(is.finite(c(n, nnew))) && nnew != n) {
          warning(":step.regr: number of rows in use has changed: \n  ",
                    nnew," observations instead of ", n)
          n <- nnew
        }
        Terms <- terms(fit)
        bAIC <- extractAIC(fit, scale, k = k, ...)
        edf <- bAIC[1L]
        bAIC <- bAIC[2L]
        ## output
        if (trace) {
            cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n",
                cut.string(deparse(formula(fit))), "\n\n", sep = "")
            utils::flush.console()
        }
        if (bAIC >= AIC + 1e-07)
            break
        nm <- nm + 1
        models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n -
            edf, change = change, AIC = bAIC)
        if (!is.null(keep))
            keep.list[[nm]] <- keep(fit, bAIC)
    }
    if (!is.null(keep))
        fit$keep <- re.arrange(keep.list[seq(nm)])
    ## !!!
##-     if (lnna) {
##-       lv <- all.vars(fit$call$formula)
##-       linan <- apply(is.na(ldata[,lv,drop=FALSE]),1,sum)>0
##-       if (any(lina!=linan))
##-           warning(":step.regr: ",# sum(lina),
##-               "observations deleted (or added) because of missing values.\n",
##-               "  You may want to refit the resulting model and/or",
##-               "call step again on it.")
##-     }
##-     fit$call$data <- lcdata
    ## !!!
    step.results(models = models[seq(nm)], fit, object, usingCp)
}
## ==========================================================================
terms2order <- function(object, squared = TRUE, interactions = TRUE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 15 Oct 2009, 16:01
  ltadd <- NULL
  if (squared) {
#    lvrs <- all.vars(formula(object)[-2])
    if (!is.list(object)||length(lterms <- object$terms)==0)
      stop("!terms2order! first argument has no 'terms' component")
    ltl <- attr(lterms, "term.labels")
    lts <- c(grep("\\^",ltl),grep(":",ltl))
    if (length(lts)) ltl <- ltl[-lts]
    ltn <- ltl[attr(object$terms,"dataClasses")[ltl]=="numeric"]
    if (length(ltn)) ltadd <- paste(paste("I(", ltn,"^2)",sep=""), collapse="+")
  }
  if (interactions) {
    ltint <- "(.)^2"
    if (is.null(ltadd)) ltadd <- ltint else ltadd <- paste(ltadd, ltint, sep=" + ")
  }
  if (is.null(ltadd)) {
    warning(":terms2order: nothing to add")
    return(formula(object))
  }
  ltadd <- paste("~.+", ltadd)
##-     attr(terms(update.formula(object, ltadd)), "term.labels")
  update.formula(object, ltadd)
}
## ==========================================================================
fitted.regr <-
  function (object, ...)
{
  class(object) <- setdiff(class(object), "regr")
  if (pmatch("fitted",names(object),nomatch=0))
    return(fitted(object, ...))
  lres <- residuals(object)
  if (inherits(lres, "condquant"))
    structure(lres[,"fit"], names=row.names(lres))  else predict(object, ...)
}
## ==========================================================================
predict.regr <-
function (object, newdata = NULL, scale = object$sigma, type = NULL, ...)
  ## bug: if used with NULL newdata, predictions will be produced
  ## for obs in model.matrix, which excludes those with missing values
{
##-   lglm <- inherits(object,"glm")
##-   lmeth <- object$call$method
##-   lnls <- length(lmeth)>0 && lmeth=="nls"
  if (length(type)==0)
    type <- if (inherits(object,"glm")) "link" else
  if (inherits(object, "polr")) "link" else "response"
  ## !!!
  if (object$fitfun=="rlm")
    if (!is.matrix(object[["x"]]))
      object$x <- model.matrix(formula(object), data=object$allvars) ## !!! was $model
##-   if (length(scale)==0) scale <- c(object$sigma,1)[1]
  class(object) <- class(object)[-1]
##-   if (missing(newdata) || length(newdata)==0) {
##-     lpred <- if (lglm)
##-       predict.glm(object, type=type, se.fit=se.fit,
##-                   dispersion=lscale^2, terms=terms, na.action = na.action )
##-     else {
##-       if (lnls) predict.nls(object, type=type, se.fit=se.fit,
##-                             na.action = na.action )  else
##-       predict.lm(object, type=type, se.fit=se.fit,
##-                     terms=terms, na.action = na.action )
##-     }
##-     lpred <- predict(object, type=type, se.fit=se.fit,
##-                   dispersion=lscale^2, terms=terms, na.action = na.action )
##-   } else {
    ldt <- newdata
    ## analyze variables
  if (!is.null(ldt))
    for (lvn in names(ldt)) {
      lv <- ldt[[lvn]]
      ## binary factors
      if (is.factor(lv)&&match(lvn,names(object$binlevels),nomatch=0)>0)
        ldt[[lvn]] <- match(lv,object$binlevels[[lvn]])-1
      else  if (is.factor(lv))   ldt[[lvn]] <- factor(lv)
    }
  if (is.null(ldt))
    predict(object, type=type, scale=object$sigma,
            dispersion=object$dispersion^2, ... )  else
    predict(object, newdata=ldt, type=type, scale=object$sigma,
            dispersion=object$dispersion^2, ... )
}
## ==========================================================================
##- extractAIC.regr <- function (fit, scale = 0, k = 2, ...)
##- {  ##- #  AIC, divided by n to allow for comparing models with different n
##- ##-   lres <- fit$residuals
##- ##-   if (is.null(lres)) {
##- ##-     lfit <- fit$fitted
##- ##-     fit$fitted.values <- nna(lfit)
##- ##-   } else  fit$residuals <- nna(lres)
##-   class(fit) <- setdiff(class(fit),"regr")
##-   extractAIC(fit, scale = scale, k = k, ...)
##- }
## ==========================================================================
vif.regr <- function(mod, cov, mmat)
{
  ## Purpose:   vif.lm  of library  car
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: modified by Werner Stahel, Date: 11 Mar 2005, 09:18
##-     v <- mod$sigma^2* mod$cov.unscaled
##-     cls <- dimnames(model.matrix(mod))[[2]]%in%dimnames(v)[[2]]
##-                                         # needed for singular cases
##-     assign <- attributes(model.matrix(mod))$assign[cls]
    cls <- dimnames(mmat)[[2]]%in%dimnames(cov)[[2]]
##-                                         # needed for singular cases
    assign <- attr(mmat,"assign")[cls]
    terms <- labels(terms(mod))
    n.terms <- length(terms)
    if (n.terms < 2) {
##-         stop("model contains fewer than 2 terms")
      return(matrix(1,1,3))
    }
    if (length(cov)==0) { # ||n.terms!=nrow(cov)|nrow(cov)!=ncol(cov)
      warning(":vif.regr: mod$cov.unscaled  is inappropriate. no vifs")
      return(matrix(NA,n.terms,3))
    }
    if (names(coefficients(mod)[1]) == "(Intercept)") {
        cov <- cov[-1, -1]
        assign <- assign[-1]
    }
    else if (mod$fitfun!="polr")
      warning("No intercept: vifs may not be sensible.")
    sd <- 1/sqrt(diag(cov))
    if (any(!is.finite(sd))) {
      warning(":vif.regr: zero variances of estimates. no R2x")
      return(NULL)
    }
    R <- cov/outer(sd,sd)
    result <- matrix(0, n.terms, 3)
    rownames(result) <- terms
    colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
    for (term in 1:n.terms) {
      subs <- which(assign == term)
      result[term, 1] <- det(as.matrix(R[subs, subs])) *
        det(as.matrix(R[-subs,-subs]))/det(R)
      result[term, 2] <- length(subs)
    }
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
    result
}
## =================================================================
drop1.multinom <-
  function (object, scope, test = c("Chisq","none"), ...)
{
    if (!inherits(object, "multinom"))
        stop("Not a multinom fit")
    if (missing(scope))
        scope <- drop.scope(object)
    else {
        if (!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)),
                "term.labels")
        if (!all(match(scope, attr(object$terms, "term.labels"),
            nomatch = FALSE)))
            stop("scope is not a subset of term labels")
    }
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>",
        scope), c("Df", "AIC")))
    ans[1, ] <- c(object$edf, object$AIC)
    if (test[1]=="Chisq") ans <- cbind(ans,Chisq=NA, p.value=NA)
    env <- environment(formula(object))
    for(i in seq(ns)) {
      tt <- scope[i]
##        cat("trying -", tt, "\n")
      nfit <- update(object, as.formula(paste("~ . -", tt)),
                     evaluate = FALSE)
      nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)
##-        nobject <- update(object, paste("~ . -", tt))
      if (nfit$edf == object$edf)
        nfit$AIC <- NA
      ans[i+1, ] <- c(nfit$edf, nfit$AIC,
                    if (test[1]=="Chisq") unlist(anova(object,nfit)[2,6:7]))
    }
    as.data.frame(ans)
}
## =================================================================
## FIXME: this is nowhere called, nor documented ...
summary.mregr <- function(object, ...)
{
  ## Purpose:   collect results for mlm object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 15 Feb 2005, 09:11
  lsum <- summary(object)
  lts <- ltp <- NULL
  for (ly in 1:length(lsum)) {
    lts <- cbind(lts,c(lsum[[ly]][["sigma"]],lsum[[ly]][["r.squared"]],
                       lsum[[ly]][["fstatistic"]]))
    ltp <- cbind(ltp,lsum[[ly]][["coefficients"]][,4])
  }
  lts[4,] <- pf(lts[3,],lts[4,],lts[5,], lower.tail=FALSE)
  lts <- lts[1:4,]
  dimnames(ltp) <- dimnames(object$coefficients)
  dimnames(lts) <- list(c("sigma","r.squared","fstatistic","p-value"),
                        dimnames(ltp)[[2]])
  list(coefficients=object$coefficients, pvalues=ltp, stats=lts)
}
## =================================================================

##- @title drop1() method for multivariate lm()`s
##- @param object
##- @param scope
##- @param test
##- @param total
##- @param add
##- @param ...
##- @note R`s  stats:::drop1.mlm() just has stop("no `drop1` method for \"mlm\" models")
drop1.mlm <- function (object, scope = NULL,
                       test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"),
                       total=TRUE, add=FALSE, ...)
{
  ## add=TRUE   do add instead of drop1
    Pillai <- function(eig, q, df.res) {
        test <- sum(eig/(1 + eig))
        p <- length(eig)
        s <- min(p, q)
        n <- 0.5 * (df.res - p - 1)
        m <- 0.5 * (abs(p - q) - 1)
        tmp1 <- 2 * m + s + 1
        tmp2 <- 2 * n + s + 1
        c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s *
            tmp2)
    }
    Wilks <- function(eig, q, df.res) {
        test <- prod(1/(1 + eig))
        p <- length(eig)
        tmp1 <- df.res - 0.5 * (p - q + 1)
        tmp2 <- (p * q - 2)/4
        tmp3 <- p^2 + q^2 - 5
        tmp3 <- if (tmp3 > 0)
            sqrt(((p * q)^2 - 4)/tmp3)
        else 1
        c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q,
            p * q, tmp1 * tmp3 - 2 * tmp2)
    }
    HL <- function(eig, q, df.res) {
        test <- sum(eig)
        p <- length(eig)
        m <- 0.5 * (abs(p - q) - 1)
        n <- 0.5 * (df.res - p - 1)
        s <- min(p, q)
        tmp1 <- 2 * m + s + 1
        tmp2 <- 2 * (s * n + 1)
        c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
    }
    Roy <- function(eig, q, df.res) {
        p <- length(eig)
        test <- max(eig)
        tmp1 <- max(p, q)
        tmp2 <- df.res - tmp1 + q
        c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
    }
##-     if (!is.null(object$drop1)) return (object$drop1)
    if (!(inherits(object, "maov") || inherits(object, "mlm")))
        stop("object must be of class \"maov\" or \"mlm\"")
    test <- match.arg(test)
    asgn <- object$assign[object$qr$pivot[1:object$rank]]
    tl <- attr(object$terms, "term.labels")
## scope
    if (is.null(scope))
        scope <- if (add) attr(terms(update.formula(object, ~(.)^2)),
                "term.labels") else
          drop.scope(object)
    else {
        if (!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)),
                "term.labels")
##-         if (!all(match(scope, tl, FALSE)))
        if (!(add||all(match(scope, tl, FALSE))))
            stop("!drop1.mlm! scope is not a subset of term labels")
    }
    ns <- length(scope)
    rdf <- object$df.residual
    res <- resid(object)
## ::: needed for finding the data later
    ldata <- eval(object$call$data,
                 envir=environment(formula(object)))
    ladd <- 1
    if (add) {
      ladd <- -1
      lna <- i.add1na(object, scope)
      if (!is.null(lna)) res[lna,1] <- NA
    }
    res <- nainf.exclude(res)
    lna <- attr(res,"na.action")
    if (!is.null(lna)) ldata <- ldata[-lna,]
## full model
    rss <- crossprod(as.matrix(res))
    rss.qr <- qr(rss)
    if (rss.qr$rank < NCOL(res))
      stop(paste("!drop1.mlm! residuals have rank", rss.qr$rank, "<",
                 ncol(res)))
    stats <- matrix(NA,length(scope),4)
    dimnames(stats) <- list(scope,c(test,"F.stat","dfnum","dfden"))
    tstfn <- switch(test, Pillai = Pillai, Wilks = Wilks,
                    "Hotelling-Lawley" = HL, HL = HL, Roy = Roy)
    object$call[[1]] <- as.name("lm")
## loop through scope
    for (lsc in scope) {
      lfo <- as.formula(paste(if (add) "~.+" else "~.-",lsc))
      lrg <- update(object, lfo, data=ldata, model=FALSE) # ,data=data
      dfj <- ladd * (lrg$df.residual - rdf)
      bss <- ladd * (crossprod(resid(lrg))-rss)
      eigs <- Re(eigen(qr.coef(rss.qr, bss),symmetric = FALSE)$values)
      stats[lsc,] <- tstfn(eigs, dfj, rdf)
    }
    ldf <- stats[1,3:4]
    names(ldf) <- c("numerator","denominator")
    if (total) {
      lpr <- predict(object)
      if (length(lna)) lpr <- lpr[-lna,] # drop rows with NA
      yy <- scale(res + lpr, scale=FALSE)
      bss <- crossprod(yy)-rss
      eigs <- Re(eigen(qr.coef(rss.qr, bss), symmetric = FALSE)$values)
      stats <- rbind(stats, "<total>"= tstfn(eigs, object$df[1], rdf))
    }
    data.frame(stats,
               p.value = pf(stats[,2],stats[,3],stats[,4], lower.tail = FALSE))
##    attr(stats,"df") <- ldf
##    stats
} ## {drop1.mlm}
## ==========================================================================
add1.mlm <-
  function (object, scope=NULL, test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"), ...)
{
  ## Purpose:    add1 for regr objects
  ## ----------------------------------------------------------------------
  drop1.mlm(object, scope=scope, test=test, total=FALSE, add=TRUE, ...)
}
## ===========================================================================
i.add1na <- function (object, scope)
{
  ##  determine rows with NA`s in model.frame for expanded model
    Terms <-
      terms(update.formula(object, paste("~.+",paste(scope, collapse = "+"))))
    fc <- object$call
    fc$formula <- Terms
    fob <- list(call = fc, terms = Terms)
    class(fob) <- oldClass(object)
    m <- model.frame(fob, xlev = object$xlevels)
    r <- cbind(resid(object))
    if (nrow(r)!=nrow(m)) {
      warning(gettextf("!add1! using the %d/%d rows from a combined fit",
                nrow(m), nrow(r)), domain = NA)
      lna <- !row.names(r)%in%row.names(m)
    }
    else lan <- NULL
  }
## ==========================================================================
## currently only called from print.regr():
print.mregr <- function(x, ...)
{
  ## Purpose:   collect results for mregr object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: Feb 2008
  f.prv <- function(x) paste(paste(names(x),x,sep=" = "),collapse=", ")
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),"\n", sep = "")
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nP-values:\n")
  print(round(x$pvalues,4))
  cat("\nStatistics for individual response models:\n")
  print(x$stats)
  cat("\nResidual degrees of freedom: ",x$df,"\n")
  ldr <- x$drop1
  if (!is.null(ldr)) {
    cat("\nMultivariate tests for all responses\n  Degrees of freedom: ",
        f.prv(attr(ldr,"df")),"\n")
    print(ldr[,])
  }
  invisible(x)
}

## ===================================================
fitted.polr <- function(object, type="link", na.action=object, ...) {
  if (pmatch(type,"link",nomatch=0)) {
    lfit <- object$linear.predictor
    if (length(lfit)==0) { # if called by original polr
      Terms <- delete.response(object$terms)
      environment(Terms) <- environment() ## ! WSt
      ## from predict.polr
      m <- object$model
      if (length(m)==0)
        m <- model.frame(Terms, object$data, na.action = function(x) x,
                         xlev = object$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, m)
      X <- model.matrix(Terms, m, contrasts = object$contrasts)
      xint <- match("(Intercept)", colnames(X), nomatch = 0)
      if (xint > 0)  X <- X[, -xint, drop = FALSE]
      lfit <- drop(X %*% object$coefficients)
    }
  } else
    object$fitted
  if (type=="class")
    lfit <- factor(max.col(lfit), levels = seq_along(object$lev),
            labels = object$lev)
  naresid(na.action,lfit)
}
## ===================================================
predict.polr <-
  function (object, newdata=NULL,
            type = c("class", "probs", "link"), ...)
  ## type link added by WSt, newdata=NULL
{
    if (!inherits(object, "polr"))
        stop("not a \"polr\" object")
    type <- match.arg(type)
    if (length(newdata)==0) {
      if (type=="link")
        eta <- fitted(object, type="link", na.action=NULL)
        Y <- object$fitted
      na.action <- object$na.action
    }
    else {
      na.action <- NULL
        newdata <- as.data.frame(newdata)
        Terms <- delete.response(object$terms)
        environment(Terms) <- environment() ## ! WSt
        m <- model.frame(Terms, newdata, na.action = function(x) x,
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(X), nomatch = 0)
    ## changed by WSt
##-         if (xint > 0)
##-             X <- X[, -xint, drop = FALSE]
        eta <- drop(X[,names(object$coefficients)] %*% object$coefficients)
        n <- nrow(X)
        q <- length(object$zeta)
##      pgumbel <- function(q) exp(pweibull(log(q))) # ???
        pfun <- switch(object$method, logistic = plogis, probit = pnorm,
            cloglog = prevgumbel, cauchit = pcauchy)
        cumpr <- matrix(pfun(matrix(object$zeta, n, q, byrow = TRUE) -
            eta), , q)
        Y <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
        dimnames(Y) <- list(rownames(X), object$lev)
    }
##-     if (newdata) && !is.null(object$na.action))
##-         Y <- napredict(object$na.action, Y)
    switch(type, class = {
        Y <- factor(max.col(Y), levels = seq_along(object$lev),
            labels = object$lev)
    }, probs = {
    }, link = { Y <- eta })
    Y <- napredict(na.action,Y)
    drop(Y)
}
## ===================================================================
condquant <- function(x, dist="normal", sig=1, randomrange=0.9)
{
  ## Purpose:   conditional quantiles and random numbers
  ## works only for centered scale families
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 24 Oct 2008, 09:30
##-   fp <- switch(dist, normal="pnorm(x,0,sig)", gaussian="pnorm(x,0,sig)",
##-                logis="plogis(x,0,sig)", revgumbel="prevgumbel(x,0,sig)",
##-                # lognormal="pnorm(log(x),0,sig)"
##-                ## x is log regsidual ?
##-                )
  if (length(x)==0) stop("!condquant! bug: no data")
  fp <- switch(dist, normal=pnorm, gaussian=pnorm, unif=function (x) x,
               logis=plogis, revgumbel=prevgumbel)
  if (is.null(fp)) stop(paste("!condquant! distribution ", dist, " not known"))
  fq <- switch(dist, normal=qnorm, gaussian=qnorm, unif=function (x) x,
               logis=qlogis, revgumbel=qrevgumbel)
##  if (NCOL(x)>=2) stop("!condquant! x must have 2 columns")
  lx <- t(apply(rbind(x),1,sort))
  lp <- fp(lx/sig)
  lpp <- rbind(rbind(lp)%*%rbind(c(0.5,0.75,0.25),c(0.5,0.25,0.75)))
  lprand <- lp[,1]+(lp[,2]-lp[,1])*
      runif(nrow(lp),(1-randomrange)/2,(1+randomrange)/2)
##  <<<<<<< .mine
  lr <- cbind(cbind(median=fq(lpp[,1]),lowq=fq(lpp[,2]),
              uppq=fq(lpp[,3]), random=fq(lprand))*sig,
              prob=lp[,2]-lp[,1])
  if (any(lp0 <- lp[,2]<=0)) lr[lp0,1:4] <- matrix(lx[lp0,2],sum(lp0),4)
  if (any(lp1 <- lp[,1]>=1)) lr[lp1,1:4] <- matrix(lx[lp1,1],sum(lp1),4)
##  dimnames(lr)[[1]] <- row.names(x)
  class(lr) <- c("condquant", "matrix")
  lr
##- =======  !!!???!!!
##-   structure(cbind(sig*cbind(median= fq(lpp[,1]), lowq  = fq(lpp[,2]),
##- 			    uppq  = fq(lpp[,3]), random= fq(lprand)),
##- 		  prob=abs(lp[,2]-lp[,1])),
##- 	    class = c("condquant", "matrix")) # "matrix", e.g., for head()
##  >>>>>>> .r32
}
## ===================================================================
residuals.coxph <- function(object, type="CoxSnellMod", ...)
{
  ## Purpose:    conditional quantiles and random numbers for censored obs
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: Aug 2010
  lres <- object$residuals
  ly <- object$y
  lst <- ly[,2]  # status
  li <- lst!=1
  ##  martingale --> coxsnell
  lres <- lst - lres
  lrs <- qnorm(exp(-lres)) # --> normalscores
  ## fill matrix with values adequate for non-censored obs
  lrr <- matrix(lrs,length(lrs),5)
  dimnames(lrr) <- list(row.names(ly),c("median","lowq","uppq","random","prob"))
  lrr[,"prob"] <- 0
  ## censoring
  if (any(li)) {
    llim <- cbind(lrs[li],Inf)
    lr <- condquant(llim, "normal")
    lrr[li,] <- lr
  }
  lrr <- cbind(lrr, fit=object$linear.predictors)
  class(lrr) <- "condquant"
  naresid(object$na.action, lrr)
}
## ===================================================================
residuals.survreg <- function(object, type="response", ...)
{
  ## Purpose:    conditional quantiles and random numbers for censored obs
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 24 Oct 2008, 10:16
  ly <- object$y  ## log for weibull
  lsig <- object$sigma
  ldist <- if (length(object$dist)>0) object$dist  else  "normal"
  li <- match(ldist, c("weibull","lognormal","loglogistic"))
  if (!is.na(li)) ldist <- c("revgumbel","normal","logistic")[li]
  ## for user-defined survreg.distributions with transformation,
  ##   this is not enough.
  if (length(lsig)==0) lsig <- summary(object)$scale
  if (length(lsig)==0) {
    warning("!residuals.survreg! no sigma found. Setting =1")
    lsig <- 1
  }
  lfit <- object$linear.predictors
  lres <- ly[,1]-lfit
  ## fill matrix with values adequate for non-censored obs
  lrr <- matrix(lres,length(lres),5)
  dimnames(lrr) <- list(row.names(ly),c("median","lowq","uppq","random","prob"))
  lrr[,"prob"] <- 0
  ## censoring
  lst <- ly[,2]  # status
##   ltt <- attr(object$response[[1]], "type")
  ltt <- attr(object$response, "type")
##-   if (length(ltt)>0&&ltt=="left") lst[lst==0] <- 2
  ltl <- length(ltt)>0&&ltt=="left"
  li <- lst!=1
  if (any(li)) {
##-     llim <- cbind(lres[li],c(Inf,NA,-Inf)[lst[li]+1]) #
      llim <- if(ltl) cbind(-Inf,lres[li]) else cbind(lres[li],Inf)
      lr <- condquant(llim, ldist, lsig)
      lrr[li,] <- lr
  }
  lrr <- cbind(lrr, fit=lfit)
  class(lrr) <- "condquant"
  structure( naresid(object$na.action, lrr), class=c("condquant", "matrix"))
}
## ===================================================================
residuals.polr <- function(object, na.action=object, ...)
{
  ## Purpose:   residuals for cumulative logit regression
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   object   result of polr
  ## Value:     list with components
  ##   median     "conditional median" residual = median of conditional distr.
  ##            of latent response variable, given observed response
  ##            will probably be replaced by mean in the near future
  ##            in order to allow for adequate smoothing
  ##   lowq, uppq  lower ans upper quartiles of this cond. distribution
  ## Remark:    experimental function !!!
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  2 Oct 2007, 11:31
  if (!((lpolr <- inherits(object, "polr"))|
        (lbin <- inherits(object, "glm")&&object$family=="binomial")))
    stop ("!residuals.polr! unsuitable first argument")
  lyr <- as.character(formula(object)[2])
  ldi <- max(match(c("model","allvars"), names(object), nomatch=0))
  if (ldi) ldt <- object[[ldi]]  else {
      ldt <- if (u.debug())  model.frame(object) else
      try(model.frame(object),silent=TRUE)
      if (class(ldt)=="try-error") 
          stop ("!residuals.polr! no data found")
  }          
  ly <- ldt[,lyr]
  if (length(ly)==0) stop ("!residuals.polr! bug: no response values found")
  if (length(dim(ly))) {
    warning(":residuals.polr: returning simple deviance residuals for non-binary (grouped) data")
    return(residuals(object, type="deviance"))
  }
  #if (lpolr)
  ly <- as.numeric(ordered(ly))
  lfit <- fitted(object, type="link", na.action=list(na.action=NULL))
  lthres <- c(-100, if (lpolr) object$zeta else 0, 100)
  llim <- cbind(lthres[ly],lthres[ly+1])-lfit
  lr <- cbind(condquant(llim,"logis"),fit=lfit,y=ly)
##-   lp <- cbind(plogis(),plogis(lthres[ly+1]-lfit))
##-   lpp <- lp%*%rbind(c(0.5,0.25,0.75),c(0.5,0.75,0.25))
##-   lprand <- lp[,1]+(lp[,2]-lp[,1])*runif(length(lfit))
##-   lr <- cbind(median=qlogis(lpp[,1]),lowq=qlogis(lpp[,2]),
##-               uppq=qlogis(lpp[,3]), random=qlogis(lprand), fit=lfit,y=ly)
  dimnames(lr)[[1]] <- row.names(ldt)
  class(lr) <- c("condquant", "matrix") 
  naresid(na.action$na.action, lr)
}
## ===========================================================================
fitcomp <- function(object, data=NULL, vars=NULL, se=FALSE,
                      xm=NULL, xfromdata=FALSE, nxcomp=51)
{
  ## Purpose:    components of a fit
  ## ----------------------------------------------------------------------
  ## !!! make nxcomp >= maximal number of factor levels !!!
  ## !!! why vars??? can possibly be much simpler!!!
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  6 Jan 2004, 22:24
##-   lf.predict <-
##-     function(object, newdata, type="linear.predictor", se.fit = FALSE) {
##-     if (type=="linear.predictor") {
##-       lmod <- model.matrix(formula(object)[-2], newdata)
##-       if (inherits(object,"polr")) lmod <- lmod[,colnames(lmod)!="(Intercept)"]
##-       rr <- lmod%*%object$coefficients[colnames(lmod)]
##-     } else rr <- predict(object, newdata, type=type, se.fit=se.fit)
##-     rr
##-   }
  if (length(data)==0) {
    data <- object$allvars
    if (length(data)==0)
      data <- eval(object$call$data)
##-     datanm <- as.character(object$call)[3]
##-     if (is.na(datanm)) stop("no data found")
##-     data <- eval(parse(text=datanm))
  }
##-   lmeth <- object$call$method
##-   lnls <- length(lmeth)>0 && lmeth=="nls"
  lnls <- c(eval(object$call$nonlinear),FALSE)[1]
##-   lvars <- unique(c(vars,all.vars(formula(object)[[3]])))
  lvars <- all.vars(formula(object)[[3]])
  if (lnls)
    lvars <- lvars[match(lvars,names(coefficients(object)),nomatch=0)==0]
  lvmiss <- match(lvars,names(data),nomatch=0)==0
  if (any(lvmiss))
    stop(paste("!fitcomp! variable(s)",
               paste(lvars[lvmiss],collapse=", "),"not in data"))
  if (length(lvars)==0)
    stop("!fitcomp! no variables selected")
  ldata <- data[,lvars,drop=FALSE]
  for (lj in 1:length(ldata))
    if (is.character(ldata[[lj]])|is.factor(ldata[[lj]]))
      ldata[,lj] <- factor(ldata[,lj])
  if (!lnls) {
    ltc <- attr(terms(object),"dataClasses")
    for (lt in lvars)
      if (lt%in%names(ltc)&&ltc[lt]=="factor") {
        lv <- all.vars(formula(paste("~",lt)))
        ldata[lv] <- lapply(ldata[lv],factor)
      }
  }
  ## generate means  xm  if needed
  if (length(xm)>0) {
    if ((!is.data.frame(xm))||any(names(xm)!=names(ldata))) {
      warning(":fitcomp: arg. xm  not suitable -> not used")
      xm <- NULL } else xm <- xm[1,]
  }
## median point and prediction for it
  if (length(xm)==0) {
    xm <- ldata[1,,drop=FALSE]
    for (lj in 1:length(ldata)) {
      lv <- ldata[,lj]
      if (is.character(lv)) lv <- factor(lv)
      xm[1,lj] <- if (is.factor(lv))
        levels(lv)[which.max(table(as.numeric(lv)))[1]]
      else   ## median(as.numeric(lv),na.rm=TRUE)
        sort(lv)[ceiling(length(lv)/2)]
        ## median should be attained in some cases
    }
  }
  if (is.null(attr(terms(object), "predvars"))) { # from  model.frame
    lterms <- attr(lm(formula(object), data=data, method="model.frame"),"terms")
    attr(object$terms,"predvars") <- attr(lterms,"predvars")
  }
  lprm <- c(predict(object, newdata=xm)) # lf.
  lny <- length(lprm)
##  expand to matrix
  if (xfromdata) {
    lx <- ldata
  } else {
    lnxc <- sapply(ldata, function(x) if (is.factor(x)) length(levels(x)) else 0)
    lnxc <- max(nxcomp, lnxc)
    lx <- ldata[1,,drop=FALSE][1:lnxc,,drop=FALSE]
    row.names(lx) <- 1:lnxc
  }
  lxm <- lx
  for (lv in names(lxm)) lxm[,lv] <- xm[,lv]
##  components
  lcomp <- array(dim=c(nrow(lx),length(lvars),lny)) # [,lvars,drop=FALSE]
  dimnames(lcomp) <- c(dimnames(lx),list(names(lprm)))
  lcse <- if (se) lcomp  else NULL
  for (lv in lvars) {
    if (xfromdata) {
      ld <- lxm
      ld[,lv] <- ldata[,lv]  # eliminate extra levels of factors
    } else { # +++
      ldv <- ldata[,lv]
      if (is.factor(ldv)) { # ---
        ldx <- factor(levels(ldv))
        lnl <- length(ldx)
        ld <- lxm[1:lnl,,drop=FALSE]
        ld[,lv] <- ldx
        lx[,lv] <- factor(c(1:lnl,rep(NA,lnxc-lnl)),labels=levels(ldv))
        ##
        lpr <- predict(object, newdata=ld, se.fit = se) # lf.
        if (se) {
          lc <- lpr$fit
          lcse[1:lnl,lv,] <- lpr$se.fit
        } else lc <- lpr
        lcomp[1:lnl,lv,] <- lc
        next # end for loop
      } else { # ---
        ld <- lxm
        lx[,lv] <- ld[,lv] <-
          seq(min(ldv,na.rm=TRUE),max(ldv,na.rm=TRUE),length=lnxc)
      } # ---
    } # +++
    ## continuous variable or xfromdata
    lpr <- predict(object, newdata=ld, se = se) # lf.
    if (se) {
      lcomp[,lv,] <- lpr$fit
      lcse[,lv,] <- lpr$se.fit
    } else lcomp[,lv,] <- lpr
  }
  if (lny==1) {
    dim(lcomp) <- dim(lcomp)[1:2]
    dimnames(lcomp) <- dimnames(lx)
    lcomp <- lcomp-lprm
    if (se) {
      dim(lcse) <- dim(lcse)[1:2]
      dimnames(lcse) <- dimnames(lx)
    }
  } else   lcomp <- sweep(lcomp,3,lprm)
  list(comp=lcomp, x=lx[,lvars,drop=FALSE], xm=xm[,lvars,drop=FALSE], se=lcse)
}
## ==========================================================================
predict.mlm <-
  function (object, newdata=NULL, se.fit = FALSE, scale = NULL, df = Inf,
    interval = c("none", "confidence", "prediction"), level = 0.95,
    type = c("response", "terms"), terms = NULL, na.action = na.pass,
##-     pred.var = res.var/wgts, wgts = 1, ...)
    pred.var = NULL, weights = 1, ...) ## ... to absorb unused args
  ## predict.lm, extended for mlm
{
    tt <- terms(object)
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        mmDone <- TRUE
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action,
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        offset <- if (!is.null(off.num <- attr(tt, "offset")))
            eval(attr(tt, "variables")[[off.num + 1]], newdata)
        else if (!is.null(object$offset))
            eval(object$call$offset, newdata)
        mmDone <- FALSE
    }
    r <- cbind(object$residuals)
    n <- nrow(r)
    m <- ncol(r)
    ynm <- colnames(r)
    p <- object$rank
    p1 <- seq_len(p)
    piv <- object$qr$pivot[p1]
    if (p < ncol(X) && !(missing(newdata) || is.null(newdata)))
        warning("prediction from a rank-deficient fit may be misleading")
##-     beta <- object$coefficients
    beta <- cbind(object$coefficients)
##-     predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
    predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv,,drop=FALSE])
    if (!is.null(offset))
        predictor <- predictor + offset
    interval <- match.arg(interval)
    if (interval == "prediction") {
        if (missing(newdata))
            warning("Predictions on current data refer to _future_ responses\n")
        if (missing(newdata) && missing(weights)) {
            w <- weights(object) ## .default
            if (!is.null(w)) {
                weights <- w
                warning("Assuming prediction variance inversely proportional to weights used for fitting\n")
            }
        }
        if (!missing(newdata) && missing(weights) && !is.null(object$weights) &&
##-             missing(pred.var)
            is.null(pred.var)
            )
            warning("Assuming constant prediction variance even though model fit is weighted\n")
        if (inherits(weights, "formula")) {
            if (length(weights) != 2L)
                stop("`weights` as formula should be one-sided")
            d <- if (missing(newdata) || is.null(newdata))
                model.frame(object)
            else newdata
            weights <- eval(weights[[2L]], d, environment(weights))
        }
    }
    type <- match.arg(type)
    if (se.fit || interval != "none") {
        res.var <- if (is.null(scale)) {
##            r <- object$residuals
            w <- object$weights
##-             rss <- sum(if (is.null(w)) r^2 else r^2 * w)
            rss <- apply( if (is.null(w)) r^2 else r^2 * w ,2,sum)
            df <- n - p
            rss/df
        }  else {
##-         scale^2
        if (length(scale)==m) scale^2 else
           stop("!predict.lm! argument scale has wrong length")
        }
        res.var.mx <- matrix(res.var, p, m, byrow=TRUE)
        if (type != "terms") {
            if (p > 0) {
                XRinv <- if (missing(newdata) && is.null(w))
                  qr.Q(object$qr)[, p1, drop = FALSE]
                else X[, piv] %*% qr.solve(qr.R(object$qr)[p1, p1])
##-                 ip <- drop(XRinv^2 %*% rep(res.var, p))
                ip <- drop(XRinv^2 %*% res.var.mx)
            }
            else ip <- rep(0, n)
        }
    }
    if (type == "terms") {
        if (!mmDone) {
            mm <- model.matrix(object)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        hasintercept <- attr(tt, "intercept") > 0L
        if (hasintercept)
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            if (!mmDone) {
                mm <- model.matrix(object)
                mmDone <- TRUE
            }
            avx <- colMeans(mm)
            termsconst <- sum(avx[piv] * beta[piv])
        }
        nterms <- length(asgn)
        if (nterms > 0) {
##-             predictor <- matrix(ncol = nterms, nrow = NROW(X))
            predictor <- array(dim=c(NROW(X),nterms,m))
##-             dimnames(predictor) <- list(rownames(X), names(asgn))
            dimnames(predictor) <- list(rownames(X), names(asgn), ynm)
            if (se.fit || interval != "none") {
##-                 ip <- matrix(ncol = nterms, nrow = NROW(X))
##-                 dimnames(ip) <- list(rownames(X), names(asgn))
                ip <- predictor
                Rinv <- qr.solve(qr.R(object$qr)[p1, p1])
            }
            if (hasintercept)
                X <- sweep(X, 2L, avx, check.margin = FALSE)
            unpiv <- rep.int(0L, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq.int(1L, nterms, length.out = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0L] <- 0L
##-                 predictor[, i] <- if (any(iipiv > 0L))
##-                   X[, iipiv, drop = FALSE] %*% beta[iipiv]
                predictor[, i,] <- if (any(iipiv > 0L))
                  X[, iipiv, drop = FALSE] %*% beta[iipiv,]
                else 0
                if (se.fit || interval != "none")
##-                   ip[, i] <- if (any(iipiv > 0L))
##-                     as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii,, drop = FALSE])^2 %*%
##-                       rep.int(res.var,p)
                  ip[, i,] <- if (any(iipiv > 0L))
                    as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii,, drop = FALSE])^2 %*%
                      res.var.mx
                  else 0
            }
            if (!is.null(terms)) {
                predictor <- predictor[, terms, drop = FALSE]
                if (se.fit)
                  ip <- ip[, terms, drop = FALSE]
            }
        }
        else {
            predictor <- ip <- matrix(0, n, 0)
        }
        attr(predictor, "constant") <- if (hasintercept)
            termsconst
        else 0
    }
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        if (is.null(pred.var)) pred.var <- res.var.mx/weights ## !!!
        hwid <- tfrac * switch(interval, confidence = sqrt(ip),
            prediction = sqrt(ip + pred.var))
        if (type != "terms") {
            if (m==1) {  ## changed
              predictor <- cbind(predictor, predictor + hwid %o% c(1, -1))
              colnames(predictor) <- c("fit", "lwr", "upr")
            } else {  ## changed
              predictor <- array(c(predictor, predictor - hwid, predictor + hwid),
                               dim=c(n,m,3))
              dimnames(predictor)[[3]] <- c("fit", "lwr", "upr")
            }
        }
        else {
            lwr <- predictor + hwid
            upr <- predictor - hwid
        }
    }
    if (se.fit || interval != "none")
        se <- sqrt(ip)
##-     if (missing(newdata) && !is.null(na.act <- object$na.action)) { ## !!! not yet extended
##-         predictor <- napredict(na.act, predictor)
##-         if (se.fit)
##-             se <- napredict(na.act, se)
##-     }
    if (m==1) { ## !!!
      if (length(dim(predictor))==3) predictor <- predictor[,,1]
      if (se.fit) if (length(dim(se))==3) se <- se[,,1]
    }
    if (type == "terms" && interval != "none") {
##-         if (missing(newdata) && !is.null(na.act)) { # !!! not yet extended
##-             lwr <- napredict(na.act, lwr)
##-             upr <- napredict(na.act, upr)
##-         }
        list(fit = predictor, se.fit = se, lwr = lwr, upr = upr,
            df = df, residual.scale = sqrt(res.var))
    }
    else if (se.fit)
        list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
    else predictor
}
## ==========================================================================
compareTerms <-
  function(..., list=NULL, seq=NULL)
{
  ## Purpose:   compare terms of several models
  ## -------------------------------------------------------------------------
  ## Arguments:
  ##   models   character vector of names of model objects
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  2010
  lcl <- match.call()[-1]
  if (!is.null(names(lcl))) lcl <- lcl[names(lcl)!="seq"]
  ## !!! bug: list does not work !!!
##-   if (length(list)) lcl <- c(as.list(lcl[names(lcl)!="list"]),
##-                              c(lcl["list"]))
  lnmod <- length(lcl)
  lmnm <- names(lcl)
  if (is.null(lmnm)) lmnm <- as.character(lcl) else
    lmnm[lmnm==""] <- as.character(lcl)[lmnm==""]
  lterms <- list()
  for (lmd in 1:lnmod) {
    lmd <- eval(lcl[[lmd]])
    if (is.character(lmd)) lmd <- get(lmd)
    ltr <- if(is.list(lmd)) attr(terms(lmd),"term.labels") else NULL
    if (is.null(ltr)) stop("!compareTerms! inadequate argument", lmnm[lmd])
    lterms <- c(lterms, list(ltr))
  }
  ltrm <- unique(unlist(lterms))
  rr <- sapply(lterms, function(x) ltrm%in%x )
  dimnames(rr) <- list(ltrm,lmnm)
  if (!is.null(seq)) rr <- rr[order(match(ltrm,seq,nomatch=length(seq)+1)),]
  rr
}
## ==========================================================================
modelTable <-
  function(models, data=NULL, seq=NULL)
{
  ## !!! Standardisierung und reestimation unterdrueckbar machen.
  ## Purpose:   collect several models into a table
  ## -------------------------------------------------------------------------
  ## Arguments:
  ##   models   character vector of names of model objects
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  9 Aug 2000, 11:50
##-   stopifnot(is.atomic(models),is.character(models))
  t.islist <- is.list(models)
  if ((!t.islist)&!(is.atomic(models)&is.character(models)))
    stop("!modelTable! Inadequate argument  model")
  t.nm <- length(models)
  t.ls <- t.trm <- t.cf <- vector("list",t.nm)
  t.nobs <- t.df <- t.sig <- rep(NA,t.nm)
  t.mnm <- if (t.islist) names(models) else models
  if (length(t.mnm)==0) t.mnm <- paste("model",1:t.nm,sep=".")
  names(models) <- t.mnm
  names(t.ls) <- names(t.sig) <- names(t.trm) <- names(t.cf) <-
    names(t.nobs) <- names(t.df) <- t.mnm
  t.data <- data
  t.vars <- names(data)
  t.trmc <- NULL
#  t.dt <- as.character(substitute(data))
  for (li in t.mnm) {
    t.r <- if (t.islist) models[[li]] else get(models[li],envir=parent.frame())
    lcls <- NULL
    if (inherits(t.r,"lm")) lcls <- "lm"
    if (inherits(t.r,"glm")) lcls <- "glm"
    if (inherits(t.r,"multinom")) lcls <- "multinom"
    if (inherits(t.r,"polr")) lcls <- "polr"
    if (inherits(t.r,"survreg")) lcls <- "survreg"
    if (is.null(lcls))
      stop(paste("!modelTable! Model ",li," is not an adequate model"))
##-     t.dat <- try(eval(t.r$call$data),silent=TRUE)
##-     if (class(t.dat)=="try-error")  if (!is.null(t.data))
##-       t.r <- update(t.r, data=t.data)  else
##-     stop("!modelTable! no data available for model", li,
##-          ". Use argument  data  to provide it.")
    t.testtype <- ifelse(inherits(t.r,"lm"),"F","Chisq")
    t.av <- t.r$allvars
    if (is.null(t.data)) {
      t.data <- eval(t.r$call$data)
      t.vars <- names(t.data)
    }
##-     if (nrow(t.data)!=nrow(t.av))
##-       stop("!modelTable! Data for model ",li,
##-            " has unsuitable number of observations")
    t.nobs[li] <- lnr <- NROW(t.r$fitted.values)
    t.df[li] <- ldf <- lnr-t.r$df.residual
    if (any(names(t.av)%nin%names(t.data))) {
      t.data <- cbind(t.data,t.av)
      t.vars <- c(t.vars, names(t.av))
    }
    t.dr <- drop1(t.r,test=t.testtype)
    t.tnm <- attr(t.r$terms, "term.labels")
##    t.tnm <- substring(t.nm,3,max(nchar(t.nm)))
    t.cf[[li]] <- t.r$testcoef[match(t.tnm,row.names(t.r$testcoef),nomatch=0),1]
    t.ls[[li]] <- t.dr
    t.trm[[li]] <- t.tnm
    t.trmc <- c(t.trmc, attr(terms(t.r),"dataClasses"))
    t.sig[li] <- c(summary(t.r)$sigma,NA)[1]
##-     t.dt <- c(t.dt, as.character(t.r$call[3]))
  }
  if (length(unique(t.nobs))>1)
    warning(":step.regr: models have different numbers of observations")
  t.data <- t.data[unique(t.vars)]
  t.tr <- unique(unlist(t.trm))
  # --- standard deviations
  t.sd <- rep(NA,length(t.tr))
  names(t.sd) <- t.tr
  ## --- standardize coefs
##-   t.dt <- unique(t.dt)
##-   if (length(t.dt)>1)
##-     warning(":modelTable: data arguments of model calls are different: "",
##-             paste(t.dt,collapse="", ""),"". Using the first one.")
##-   t.data <- get(t.dt[1])
  t.trmc <- t.trmc[t.tr]
  t.trmc <- t.trmc[!is.na(t.trmc)]
  t.trn <- names(t.trmc[t.trmc=="numeric"])
  if (length(t.trn)) {
    t.form <- as.formula(paste("~-1+",paste(t.trn,collapse="+")))
    t.x <- model.matrix(t.form, data=t.data)
    t.vr <- apply(t.x,2,var)
    t.sd[t.trn] <- sqrt(t.vr)
  }
  ## --- coefs and p values
  t.nt <- length(t.tr)
  t.coef <- matrix(NA,t.nt,t.nm)
  dimnames(t.coef) <- list(t.tr,t.mnm)
  t.pr <- t.coef
  for (li in t.mnm) {
    t.dr <- t.ls[[li]][-1,]
##-     t.t <- substring(row.names(t.dr),3,max(nchar(row.names(t.dr))))
    t.pvcol <- grep("Pr\\(",names(t.dr))
    if (is.null(t.pvcol))
      stop("!modelTable! no column with p-values found")
    t.t <- row.names(t.dr)
    if (length(t.t)) {
      t.pr[t.t,li] <- t.dr[,t.pvcol]
      t.coef[t.trm[[li]],li] <- 0
      t.c <- t.cf[[li]]
      t.sdc <- t.sd[names(t.c)]
      t.isdc <- !is.na(t.sdc)
      t.c[t.isdc] <- t.c[t.isdc]*t.sdc[t.isdc]
      t.coef[names(t.c),li] <- t.c
    }
  }
  ## reorder
  if (length(seq)>0) {
    t.i <- match(seq, t.tr)
    t.t <- c(t.tr[t.i[!is.na(t.i)]],t.tr[!t.tr%in%seq])
    if (any(!t.tr%in%t.t|!t.t%in%t.tr)) warning("bug. not all terms. ask wst")
    t.coef <- t.coef[t.t,]
    t.pr <- t.pr[t.t,]
    t.sd <- t.sd[t.t]
  }
#  attr(t.coef,"standardized") <- t.trn
  if (all(is.na(t.sig))) t.sig <- NULL
  t.r <- list(coef=t.coef, p=t.pr, sd.terms=t.sd, sigma=t.sig, nobs=t.nobs,
              df=t.df)
  class(t.r) <- "modelTable"
  t.r
}
## ==========================================================================
"[.modelTable" <- function(object,rows=NULL,cols=NULL, reduce=TRUE) {
  if (is.null(rows)) rows <- 1:nrow(object$coef)
  if (is.null(cols)) cols <- 1:ncol(object$coef)
  lp <- object$p[rows,cols,drop=FALSE]
  li <- if(reduce) !apply(is.na(lp),1,all)  else 1:nrow(lp)
  if (length(li)==0) stop("![.modelTable! no terms left")
  lsd <- if (length(object$sd.terms)) object$sd.terms[rows][li]
  lsig <- if (length(object$sigma)) object$sigma[rows][li]
  rr <- list(coef=object$coef[rows,cols,drop=FALSE][li,,drop=FALSE],
             p=lp[li,,drop=FALSE], sd.terms=lsd, sigma=lsig,
             nobs=object$nobs[cols], df=object$df[cols])
#  attr(rr$coef,"standardized") <- attr(object$coef,"standardized")
  class(rr) <- class(object)
  rr
}
## ==========================================================================
format.modelTable <-
  function(x, digits=getUserOption("digits"),
           stars = c("***","** ","*  ",":  ",".  "), sep="", ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   x        a modelTable object
  ##   tex      if TRUE, the output will be suitable for pasting into
  ##            (la)tex source
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 23 Dec 2008, 10:09
  t.sd <- x$sd.terms
  lsd <- length(t.sd)>0
  t.stars <- c(stars, rep(" ",5))[1:5]
  t.st <- stars[cut(x$p,c(-1,0.001,0.01,0.05,0.1,1))]
  t.st[is.na(t.st)] <- "   "
  dim(t.st) <- dim(x$p)
##-   t.stdd <- attr(x$coef,"standardized")
##-   if (is.null(t.stdd)||!t.stdd)
##-     warning(":: Coefficients are not standardized",call.=FALSE)
  t.cf <- x$coef
  if (length(x$sigma)) {
    t.cf <- rbind(t.cf, .sigma.=x$sigma)
    t.st <- rbind(t.st, "   ")
    if(lsd) t.sd <- c(t.sd, .sigma.=NA)
  }
  t.cf <- rbind(t.cf, .df.=x$df)
  t.st <- rbind(t.st, .df.="  ")
  if(lsd) t.sd <- c(t.sd, .df.=NA)
  if (lnobs <- length(unique(x$nobs))>1) {
    t.cf <- rbind(t.cf, .nobs.=x$nobs)
    t.st <- rbind(t.st, "   ")
    if(lsd) t.sd <- c(t.sd, .nobs.=NA)
  }
  t.cfo <- format(t.cf, digits=digits)
  t.cfo[".df.",] <- paste("  ",format(t.cf[".df.",]))
  if (lnobs) t.cfo[".nobs.",] <- paste("",format(t.cf[".nobs.",]))
  t.cfo[grep("NA",t.cfo)] <- paste(c(rep(" ",2+digits/2),"-"),collapse="")
  t.ii <- (!is.na(t.cf))&t.cf==0
  t.cfo[t.ii] <- {
    t.cf0 <- gsub("0"," ",t.cfo[t.ii])
#    t.cf0 <- sub("\\.","+",t.cf0)
  }
  t.o <- paste(sep, t.cfo, sep, t.st)
  t.sdo <- if (lsd) paste(sep, sub("NA","  ", format(t.sd, digits=digits))) else NULL
  t.out <- cbind(t.sdo, matrix(t.o, nrow=nrow(t.cf)))
  t.nm <- row.names(t.cf)
  if (lsd) t.nm <- paste(t.nm,ifelse(!is.na(t.sd),"@",""))
  dimnames(t.out) <- list(t.nm,c(if (lsd) "sd",dimnames(x$p)[[2]]))
  t.out
}
## ==========================================================================
print.modelTable <- function(x, tex = FALSE, transpose=FALSE,...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   x        a modelTable object
  ##   tex      if TRUE, the output will be suitable for pasting into
  ##            (la)tex source
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 23 Dec 2008, 10:09
  if (tex) {
    sep <- "&"
    t.out <- format.modelTable(x, sep=sep, ...)
    if (transpose) t.out <- t(t.out)
    mc1 <- "&\\mc{2}{c}{"
    mc2 <- "}"
    end <- "\\\\"
    headend <- "\\\\ \\hline"
    tabstart <- "  \\begin{tabular}" # \\providecommand{
    tabend <- "  \\end{tabular}"
    t.end <- c(rep(end,nrow(x$p)), headend, "")
    cat(tabstart, "{", rep("rl",ncol(t.out)),"}\n  ",
        paste(sep,mc1,colnames(t.out),mc2,sep=""),headend,"\n")
    for (li in 1:nrow(t.out)) cat("  ",t.out[li,],t.end[li],"\n")
    cat(tabend,"\n")
  } else {
    t.out <- format.modelTable(x)
    if (transpose) t.out <- t(t.out)
    colnames(t.out) <- paste("",colnames(t.out))
    print(t.out, quote=FALSE)
  }
  if (length(x$sd.terms) && any(!is.na(x$sd.terms))) cat("\n  @ : standardized coefficients\n")
  invisible(x)
}
## ===================================================================
leverage <- function(fit)
{
  ## Purpose:  extract leverages
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 10 Nov 2010, 09:31
  lh <- fit$leverage
  names(lh) <- names(fit$resid)
  lnaa <- fit$na.action
  if (length(lnaa)) naresid(lnaa, lh) else lh
}

## ==========================================================================
## pseudoreplicate variability
xdistResdiff <- function(object, perc=c(3,10,80), trim=0.1, nmax=100,
                         nsim=100, out="aggregate")
{
  ## Purpose:   distance in x space and absolute residual difference
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 13 Oct 2011, 09:40
  if (!inherits(object,"lm")) stop("only suitable for lm like objects")
  lres <- object$stres # resid
  lqr <- object$qr
  ln <- nrow(lqr$qr)
  lq <- qr.qy(lqr, diag(1, nrow = ln, ncol = lqr$rank)) # like in hat
  li <- which(!is.na(lres))
  if (ln>nmax) {
    li <- sample(1:ln, nmax) # [replace=FALSE]
    lq <- lq[li,]
    lres <- lres[li]
    ln <- nmax
  }
  ldist <- dist(cbind(lq))
  lrd <- abs(outer(lres,lres,"-"))
  if (nsim) {
    lrsim <- matrix(NA,length(ldist),nsim)
    for (ls in 1:nsim) {
      li <- sample(ln)
      lrsim[,ls] <- as.dist(lrd[li,li])
    }
  }
  lnm <- names(lres)
  lm <- diag(ln)
  lid <- cbind(id1=lnm[rep(1:(ln-1),(ln-1):1)],
               id2=lnm[row(lm)[row(lm)>col(lm)]])
  lrd <- as.dist(lrd)
  lio <- order(ldist)
  rr <- data.frame(lid[lio,], xdist=ldist[lio], resdiff=lrd[lio])
  if (nsim) rr <- cbind(rr, rdsim=lrsim[lio,])
  class(rr) <- c("xdistResdiff","data.frame")
  if (out=="aggregate")  xdistResscale(rr, perc=perc)  else rr
}
## ====================================================================
xdistResscale <- function(x, perc=c(3,10,90), trim=1/6)
{
  ## Purpose:  aggregate  xdistResdiff  data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 18 Oct 2011, 08:40
  if (!inherits(x,"xdistResdiff"))
    stop ("only programmed for  xdistResdiff  objects")
  llim <- c(-0.001*max(x$xdist), x$xdist[perc*nrow(x)/100], max(x$xdist))
  lgrp <- list(cut(x$xdist, llim))
  lrd <- aggregate(sqrt(x$resdiff), lgrp, mean, trim=trim)[,2]
  lmn <- mean(sqrt(x$resdiff), trim=trim)
  ljsim <- FALSE
  if (ncol(x)>5) {
    ljsim <- TRUE
    lrdsim <- aggregate(sqrt(as.matrix(x[,-(1:4)])), lgrp, mean, trim=trim)[,-1]
    lnsim <- ncol(lrdsim)
    lrdmn <- apply(lrdsim,1,mean)
    lrdse <- apply(lrdsim,1,sd)
    lrss <- apply(sweep(lrdsim,1,lrdse,"/")^2,2,sum)
    lrss0 <- sum((lrd/lrdse)^2)
##-     ltestall <- sum(((lrd-lmn)/lrdse)^2)
##-     lpv <- pchisq(ltestall, df=length(lrd)-1, lower=FALSE)
##-     lpv1 <- pnorm((lrd[1]-lmn)/lrdse[1]) # one-sided is good
    lpv <- apply(lrdsim<lrd,1,sum)/lnsim
    names(lpv) <- paste("pv",1:length(lpv),sep="")
    lpv <- c(lpv, pv.rssq=sum(lrss0<lrss)/lnsim)
  }
  rr <- cbind(xdist=aggregate(x$xdist, lgrp, mean)[,2], resd.mean=lrd)
  if (ljsim) rr <- cbind(rr, resd.simmean=lrdmn, resd.se=lrdse)
  attr(rr,"limits") <- llim
  attr(rr,"resdMean") <- lmn
  attr(rr,"trim") <- trim
  attr(rr,"perc") <- perc
  if (ljsim) attr(rr,"pvalues") <- lpv #c(shortdist=lpv1, overall=lpv)
  class(rr) <- c("xdistResscale", "matrix")
  rr
}
## =======================================================================
plot.xdistResscale <- function(x, xlab="distance in x space",
       ylab="average abs. residual difference", col.aux="grey30", ...)
{
  ## Purpose:   plot average residual difference^2 vs. x distance
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 14 Oct 2011, 08:46
  lxdist <- sqrt(x[,"xdist"])
  lrd <- x[,"resd.mean"]
  ly2 <- lrd+2*x[,"resd.se"]
  ljse <- length(ly2)>0
  llim <- attr(x,"limits")
  lymax <- if (ljse) max(ly2) else max(lrd)
  plot(lxdist, lrd, xlab=xlab, ylab=ylab,type="b",
       xaxs="i", xlim=c(0,sqrt(max(llim))), yaxs="i", ylim=c(0,1.02*lymax),
       axes=FALSE, ...)
  box()
  axis(2)
  lxat <- pretty(c(0,0.8*last(llim)))
  axis(1, at=sqrt(lxat), labels=format(lxat))
  if (ljse) {
    segments(lxdist,lrd-2*x[,"resd.se"],lxdist,ly2,
             col=col.aux)
    lysim <- x[,"resd.simmean"]
    if (length(lysim)) {
      lxd <- sqrt(max(llim))/(nrow(x)*10)
      segments(lxdist-lxd, lysim, lxdist+lxd, lysim, col=col.aux)
    }
  }
  axis(3,at=c(0,sqrt(llim[-1])),labels=rep("",length(llim)), col=col.aux)
  abline(h=attr(x,"resdMean"), lty=3, col=col.aux)
  "plot.xdistResscale done"
}
## ==========================================================================
## plotting functions
## ==========================================================================
plot.regr <-
function(x, data=NULL, plotselect = NULL, sequence=FALSE,
         xplot = TRUE, x.se=FALSE, addcomp=FALSE, glm.restype = "deviance",
         condprobrange=c(0.05,0.8), leverage.cooklim = 1:2, 
         weights = NULL, wsymbols=NULL, symbol.size=NULL,
         markprop=NULL, lab=NULL, cex.lab=0.7, mbox = FALSE, jitterbinary=TRUE,
         smooth = TRUE, x.smooth = smooth, smooth.par=NA, smooth.iter=NULL,
         smooth.sim=19, nxsmooth=51, 
         multnrows = 0, multncols=0,
         lty = c(1,2,5,3,4,6,1,1), lwd = c(1,1,2,1,1.5,1,1,1),
         colors = getUserOption("colors.ra"), pch=NULL, col=NULL,
         main = NULL, cex.title = NULL,
         reslim = TRUE, reslimfac=4.0, reslimext=0.1, 
         resaxp=NULL, stresaxp=NULL,
         ylim=TRUE, ylimfac=3.0, ylimext=0.1, yaxp=NULL, 
         mf = NULL, mfcol=FALSE, mar=c(3,3,2,1), mgp=c(2,0.7,0),
         oma = 2*(prod(mf)>1), cex=par("cex"), ask = NULL,
         ...)
{
## Purpose:  more plots for residual analysis
## -------------------------------------------------------------------------
## Arguments:
##  x     result of regr, lm or glm (or nls?)
##  data
##  markprop      how many extreme positive and negative residuals should be
##             labeled?
##  lab        labels for points
##  cex.lab     character expansion for labels
##  mf         multiple frames: number of rows and columns
##    oma arguments passed to par
##  ask        graphics pararmeter: should R ask before erasing the screen?
##  lty, colors, lwd  to be used for different elements of the plots:
##     [1]         observations
##     [2]         reference lines
##     [3]         smooth
##     [4]         simulated smooths
##     [5]         terms in plresx
##     [6]         for confidence bands of terms
##     [7]         for conditional distribution of residuals: simulated points
##     [8]         ... interquartile ranges
##     [9]         ... medians
##  main       main title, defaults to the formula of  x
##             if it starts by : it will be appended to the formula
##  cex.title  character expansion for the main title
##  wsymbols   plot points by circles according to weights (of x)
##  symbol.size  ... to be used for these symbols
##  plotselect which plots should be shown?
##             0 : do not show,  1: show without smooth,  2: show with smooth
##             The default is
##             c( yfit=0, ta=3, tascale = NA, weights = NA, qq = NA, leverage = 2)
##             modify this vector to change the selection and the sequence in
##             which the plots appear
##    yfit     response against fitted
##    ta       Tukey-Anscombe plot (residuals against fit)
##    tascale  absolute standardized residuals against fit
##             (smooth is calculated from sqrt(abs(stres)) )
##    weights  standardized residuals against weights to be given in
##             argument weights or as x$weights
##    qq       normal plot of standardized residuals
##    leverage    residuals against leverage
##    resmatrix  scatterplot matrix of residuals for multivariate regression
##    qqmult   QQ plot of Mahalanobis lengths for multivariate regression
##  weights    weights to be used for the plot of residuals against weights
##             (Note: x$weights are used for symbol sizes in any case)
##  sequence   plot residuals on sequence?
##             if NA, this is done unless another variable represents the seq
##  xplot      if TRUE, show residuals against terms in the model
##             if a formula or a character vector, it contains the terms
##             against which the residuals are plotted
##             may contain variables that are not in the model
##  reslim    plotting limits for residuals
##  x.se       draw standard error bands in term plots
##  x.smooth   if TRUE, show smooth in term plots
##  x.ylim     y range(s) for terms plots
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date:  7 May 93 / 2002
  lnaaction <- x$na.action
  x$na.action <- NULL
## family
  if (inherits(x,"mulltinom"))
    stop("!plot.regr! I do not know how to plot results of a mulitnomial regresseion")
  lfam <- x$distrname
  if (length(lfam)==0) lfam <- x$family$family
  if (is.null(lfam) || lfam=="" || is.na(lfam)) lfam <- "gaussian"
  lfgauss <- lfam%in%c("gaussian","Gaussian")
  lglm <- inherits(x, "glm")
##  lpolr <- inherits(x,"polr")
  lnnls <- !inherits(x, "nls")
  lfcount <- lfam == "binomial" | lfam == "poisson" | lfam == "multinomial" |
             inherits(x,"polr")
  lform <- formula(x)
## weights
  lweights <- x$weights
  lwgt <- length(lweights)>1 && any(lweights!=lweights[1],na.rm=TRUE) # are there weights?
## weights used for plotting ?
  lwsymbols <- is.logical(wsymbols)&&wsymbols
  if (lwsymbols&!lwgt) {
    warning(":plot.regr: I do not find weights in model object.")
    lwsymbols <- FALSE }
  if (length(wsymbols)==0||is.na(wsymbols)) lwsymbols <- lwgt
## plot selection
  lplsel <- c( yfit=0, ta=3, tascale = NA, weights = NA, qq = NA,
              leverage = 2, resmatrix = 1, qqmult = 3)
  if (length(plotselect)>0) {
    lpls <- TRUE
    lplnm <- names(plotselect)
    if (length(lplnm)==0) {
      if (length(plotselect)==length(lplsel))
        lplnm <- names(lplsel)
      else {
        warning(":plot.regr: Inadequate argument plotselect")
        lpls <- FALSE}
    }
    if (lpls) {
      if ("default"%in%lplnm) {
        lplnm <- setdiff(lplnm, "default")
        lplsel[] <- if (plotselect["default"]==0) 0 else
        pmin(lplsel,plotselect["default"])
      }
      lina <- is.na(match(lplnm,names(lplsel)))
      if (any(lina)) {
        warning(":plot.regr: Inadequate elements in plotselect: ",
              paste(names(plotselect)[lina], collapse=", "))
        lplnm <- lplnm[!lina] }
      lplsel[lplnm] <- plotselect[lplnm]
    }
  }
  lseq <- NA
## -----------------------------------
## prepare objects needed for plotting
  rtype <- "response"
  if (lglm) rtype <- glm.restype
  if (inherits(x,"coxph")) rtype <- "martingale"
##-   if (lsurv) rtype <- "condquant"
##-   lcondq <- !is.na(pmatch(rtype,"condquant"))
##-   if(!lcondq)
##-     lres <- residuals(x, type=rtype) else {
##-       lres <- if (lpolr) residuals.polr(x) else residuals.survreg(x)
##-       lcondq <- NCOL(lres)>1
##-     }
  lres <- if (substring(rtype,1,4)=="cond")  residuals.polr(x)
    else residuals(x, type=rtype)
  lcondq <- inherits(lres,"condquant")
  if (NCOL(lres)==1) lres <- cbind(lres)  ## do not destroy class attr
  if (var(lres[,1],na.rm=TRUE)==0)
    stop("!plot.regr! all residuals are equal -> no residual plots")
  lmres <- if (lcondq) 1 else ncol(lres)
  lmult <- lmres>1
  lyname <- deparse(lform[[2]])
    if (nchar(lyname)>10) lyname <- "Y"
  if (lmult) {
    if (is.null(lcn <- colnames(lres))) lcn <- 1:ncol(lres)
    colnames(lres) <- paste("res",lcn,sep=".")
  lyname <- colnames(lres)
  }
##-   if (length(lyname)==0) lyname <- paste("Y",1:lmres,sep="")
  IR(any(is.na(lres)))
  ln <- nrow(lres)
  lrname <- paste("res(", lyname, ")")
##-   ldfres <- x$df.residual
  ldfres <- x$df
  if (is.null(ldfres)) ldfres <- x$df.residual
  if (is.null(ldfres)) {
    warning(":plot.regr: bug: no df of residuals. setting n-1")
    ldfres <- ln-1
  }
  if(length(ldfres)>1) ldfres <- ldfres[2]
  ## fit
  lf <- x$linear.predictors  # predict(x)  does not work for nls
  lfname <- "Linear Predictor"
  if(is.null(lf)) {
    lf <- fitted(x)
    lfname <- "Fitted Values"
  }
  lf <- cbind(lf)
  lsigma <- x$sigma
  if (length(lsigma)==0) lsigma <- x$scale
  if (length(lsigma)==0)
    lsigma <- if (lfcount) 0 else sqrt(apply(lres^2,2,sum)/ldfres)
##  lsigma <- x$sigma
## weights
  if (lwgt) {
    if (length(lweights)!=ln) {
      warning(c(":plot.regr: there is something wrong with the weights.",
                "They are not used."))
      lwgt <- lwsymbols <- FALSE }
  } else lweights <- NULL
## hat
  lhat <- x$leverage
  if(length(lhat)==0&&length(x$qr)>0) lhat <- hat(x$qr)
##-   if (lwgt) lhat <- lhat/lweights use almost-unweighted h for plotting
## standardized residuals
  lstres <- x$stres
  if (length(lstres)==0) {
    lstres <- if (all(lsigma>0) & !inherits(lres,"condquant")) {
      if (length(lhat)==0||any(is.na(lhat))) sweep(lres,2,lsigma,"/") else
      lres/outer(ifelse(lhat>=1,1,sqrt(pmax(0,1-lhat))),lsigma)
    } else  lres
  }
  lstres <- cbind(lstres)
##-     if (lwgt) lstres <- lstres*sqrt(lweights) ## !!! check
  lrabs <- abs(lstres)
  lstrname <- paste("st", lrname, sep = ".")
  if (lmult) {
    lresmd <- x$resmd
    if (is.null(lresmd)) lresmd <- mahalanobis(lres,0,var(lres))
  }
## smooth
  if (is.null(smooth.iter)||is.na(smooth.iter)) smooth.iter <- {
    ldn <- x$distrname
    50 * if (is.null(ldn)) 1 else !ldn%in%c("binomial","multinomial")
  }
  if (is.logical(smooth)) smooth <- i.smooth
  lsmpar <- if (is.na(smooth.par)) 5*ln^log10(1/2)*(1+lglm) else
               smooth.par# 2
# simulated residuals
  lnsims <- smooth.sim
#  if (lcondq) lnsims <- 0
  if (length(lnsims)==0) lnsims <- 0
  lnsims <- if (is.logical(lnsims)&&lnsims) 19 else as.numeric(lnsims)
  if (!lnnls) lnsims <- 0
  lfitfun <- x$fitfun
  if ((!is.null(lfitfun)&&lfitfun=="survreg")) lnsims <- 0
  if (lnsims>0 & !class(x)[1]%in%c("regr","lm","glm")) {
    warning(":plot.regr/simresiduals: ",
            "I can simulate only for `regr`, `lm` and `glm` objects")
    lnsims <- 0
  }
  if (lmult) lnsims <- 0 # not yet programmed for mlm
  lsimres <- NULL
  if (lnsims>0) {
      lsimr <- if(u.debug())  simresiduals(x, lnsims)  else
      try(simresiduals(x, lnsims),silent=TRUE)
    if (class(lsimr)=="try-error") 
      warning(":plot.regr/simresiduals: simresiduals did not work. ",
              "No simulated smooths")
    else {
      lsimres <- lsimr$simres
      lsimstres <- lsimr$simstres
    }
    if (length(lsimres)==0) lnsims <- 0
  }
## plot range
  if (is.logical(reslim)&&!is.na(reslim))
    if (reslim) {
      reslim <- if(lcondq) robrange(c(lres[,1:4]), fac=reslimfac) else
      apply(lres,2,robrange, fac=reslimfac)
      streslim <- reslim*mean(abs(lstres))/mean(abs(lres))
    } else  reslim <- streslim <- NULL
  if (!is.null(reslim))
    if (any(dim(cbind(reslim))!=c(2,lmres))) {
      warning(":plot.regr: unsuitable argument  reslim ")
      reslim <- NULL
    }
## color vector
  if (length(col)>ln) if(length(lnaaction)>0 && max(lnaaction)<=length(col))
    col <- col[-lnaaction]
## labels
  ## priorities:  lab , pch , row.names
  lpch <- pch
  if (length(lpch)==ln) { ## labels given in pch
    if (length(lab)==0) lab <- lpch
    lpch <- NULL }
  if (length(lpch)==0)
    lpch <- if (lcondq) ifelse(lres[,"prob"]==0,15,3) else
      3 # ifelse(ln>200,"+",3)
  lrown <- row.names(lres)
  if (length(lrown)==0) lrown <- as.character(1:ln)
  
  if (length(lab)>1&length(lab)!=ln) {
    warning(":plot.regr: argument  lab  has unsuitable length")
    lab <- NULL
  }
  llabels <- lrown
  if (length(lab)>0) {
    llabels <- if (is.character(lab)) lab else as.numeric(lab) # factors
    if (length(llabels)>ln) {
      if(length(lnaaction)>0 && max(lnaaction)<=length(llabels))
        llabels <- llabels[-lnaaction]
    }
    else
      llabels <- rep(llabels, length=ln)
}
  ## now, llabels always useful
  if (length(markprop)==0 || is.na(markprop))
      markprop <- if (length(lab)>1|lcondq) 0 else ceiling(sqrt(ln)/2)/ln
  if (markprop==0) {
    llab <- rep(lpch, length=ln)
  } else  {
    llab <- llabels
    if (markprop<1) {
      li <- if (lmult)  order(lresmd)[1:(ln*(1-markprop))]  else
      order(abs(lstres[,1]))[1:(ln*(1-markprop))]  # [,1] for condq
      llab[li] <- NA
    }
}
  ## llabna will contain NA where weights should determine symbol
  llabna <- if (lwgt&markprop==0) rep(NA, length=ln) else llab 
## weights to be used for symbols
## plot symbols
  if (length(symbol.size)==0||is.na(symbol.size)) symbol.size <- 3/ln^0.3
  if (lwsymbols) {
    liwgt <- is.na(llabna) ## labels precede weights: which weights are used?
    if (any(liwgt)) {
      lweights <- lweights/mean(lweights,na.rm=TRUE)
      lsyinches <- 0.02*symbol.size*par("pin")[1] # *max(lweights,na.rm=TRUE)
      llab[liwgt] <- ifelse(is.character(llab),"",0) ## could be NA?
    } else lwgt <- FALSE
  }
  ltxt <- is.character(llab) # &&any(nchar(llab)>1)
  lpty <- ifelse(ltxt,"n","p")
## -----------------------------------
## prepare graphical elements
  lty <- rep(c(lty,1:6),length=6)
  lwd <- rep(c(lwd,1),length=6)
  if (length(ask)==0) {
      ask <- getOption("ask")
      if(length(ask)==0)
          ask <- !last(c("",names(dev.list())))%nin%c("postscript","pdf")
  }
  lnewplot <- FALSE
  if (is.null(mf)) mf <- if (lmult) {
    if (lmres<=4) c(lmres, min(3,lmres)) else {
      lnewplot <- TRUE
      lmr1 <- ceiling(sqrt(lmres))
      c((lmres-1)%/%lmr1+1, lmr1)
    }} else c(2,2)
  loma <- if (length(oma)==4) oma else c(0,0,oma[c(1,1)])
  loldpar <- NULL
  if (length(mf)==2) {
    loldpar <- if (mfcol)
      par(ask = ask, mfcol=mf, oma = loma, cex=cex, mar=mar, mgp=mgp)  else
      par(ask = ask, mfrow=mf, oma = loma, cex=cex, mar=mar, mgp=mgp)
  } else par(ask = ask, cex=cex, mar=mar, mgp=mgp)
  colors <- if (length(colors)==0)
          c("black","gray","blue","cyan","darkgreen","green",
            "burlywood4","burlywood3","burlywood4")
  else
    rep(colors,length=9)
  if (is.na(colors[4])) colors[4] <- colorpale(colors[4])
  if (is.na(colors[6])) colors[4] <- colorpale(colors[6])
  lftext <- paste(as.character(lform)[c(2,1,3)],collapse="")
  if (length(main)==0) main <- lftext
  if (is.logical(main)) main <- if (main) lftext else ""
  main <- if (is.character(main) && substring(main,1,1)==":")
    paste(lftext,substring(main,2,30)) else as.character(main)
  if (length(main))  tit(main) <- tit(x)
  if (is.null(cex.title)) cex.title <- max(0.5, min(1.2,
      par("mfg")[4]*par("pin")[1]/(par("cin")[1]*nchar(main))))
  outer.margin <- par("oma")[3]>0
## --------------------------------------------------------------------------
  ## start plots
  if (!lmult) lplsel[c("resmatrix","qqmult")] <- 0
  lplsel <- lplsel[is.na(lplsel)|lplsel>0]
  if (length(lplsel))
for (liplot in 1:length(lplsel)) {
  lpllevel <- lplsel[liplot]
  lpls <- names(lpllevel)
  if (lnewplot) par(ask = ask, mfrow=mf, oma = loma, cex=cex, mar=mar, mgp=mgp)
##  y on fit
  if(lpls=="yfit") {
    if (is.na(lpllevel)) lpllevel <- 3*(lplsel["ta"]==0)
    if (is.na(lpllevel)) lpllevel <- 0
    if (lpllevel>0) {
      if (lfam=="multinomial")
        warning(":plot.regr: I cannot plot categorical Y on fit")
      else {
      ly <- x[["y"]]
      if (length(ly)==0) ly <- cbind(y=lf + lres)  else
        if (nrow(ly)>length(lf))
          if(length(lnaaction))
          ly <- ly[-lnaaction,,drop=FALSE]
      lsimr <- c(lf)+lsimres
      if (is.logical(ylim)&&!is.na(ylim)) ylim <-
        if (ylim) apply(ly,2,robrange, fac=ylimfac) else  NULL
      i.plotlws(lf,ly, lfname,lyname,main,outer.margin,cex.title,
                colors,lty,lwd,lpch, col,
                lpty, llabna,cex.lab,ltxt, lwsymbols,lwgt,lweights,liwgt,
                lsyinches, lpllevel>1,smooth,lsmpar,smooth.iter,
                ylim=ylim, ylimext=ylimext, yaxp=yaxp,
                reflinex=0,refliney=1,lnsims=lnsims, simres=lsimr,
                condprobrange=condprobrange)
  } } }
## ---
  if(lpls=="ta") {
    if (is.na(lpllevel)) lpllevel <- 3*(lplsel["yfit"]==0)
    if (is.na(lpllevel)) lpllevel <- 3
    if (lpllevel>0) {
      lrx <- rbind(apply(lf,2,mean,na.rm=TRUE))
      lry <- rbind(rep(-1,lmres))
      i.plotlws(lf,lres, lfname,lrname,main,outer.margin,cex.title,
                colors,lty,lwd,lpch, col, lpty, llabna,cex.lab,ltxt,
                lwsymbols,lwgt,lweights,liwgt,lsyinches,
                lpllevel>1,smooth,lsmpar,smooth.iter,
                ylim=reslim, ylimext=reslimext,
                reflinex=lrx,refliney=lry,
                lnsims=lnsims, simres=lsimres,
                condprobrange=condprobrange)
    }
  }
## ---
  if(lpls=="tascale"&&!lcondq) {
    if (is.na(lpllevel)) lpllevel <- 3*(!lfcount)
    if (lpllevel>0)
      i.plotlws(lf,lrabs, lfname, paste("|",lstrname,"|"),
        main, outer.margin, cex.title,
        colors,lty,lwd,lpch, col, lpty, llabna,cex.lab,ltxt,
        FALSE,FALSE, # lplwgt=F,wgt=F,
        lweights,liwgt,lsyinches, lpllevel>1, smooth,lsmpar,smooth.iter,
        smooth.power=0.5,
        ylim=c(0,max(abs(streslim))), yaxs="i",
                condprobrange=c(0.01,0))
    if (lpllevel>2&lnsims>0&&length(lsimstres)>0) {
      lio <- order(lf)
      for (lr in 1:lnsims) {
        ys <- smooth(lf, sqrt(abs(lsimstres[,lr])),
                     weights=lweights, par=lsmpar, iter=smooth.iter)
        if (length(ys))
          lines(lf[lio], ys[lio]^2, lty=lty[4], lwd=lwd[4], col=colors[4])
      }
    }
  }
## --- plot abs. res vs. weights
  if(lpls=="weights") {
  if (is.na(lpllevel)) lpllevel <- 3*(lfgauss&lwgt)
##  lplweights <- (is.logical(weights) && weights) | (length(weights)>1)
  if (lpllevel>0) { # plot on weights
    lwgts <- if (length(weights)>1) {
      if (length(weights)==ln) weights  else {
        if (length(weights)>ln)
          if(length(lnaaction)>0 && max(lnaaction)<=length(weights))
            weights[-lnaaction]
      }
    } else x$weights
    if (length(lwgts)!=ln)
      warning(":plot.regr: no weights found. cannot plot absres on weights")
    else {
      i.plotlws(lwgts,lrabs, "weight (log)",
        paste("|",lstrname,"|"),main,outer.margin, cex.title,
        colors,lty,lwd,lpch, col, lpty, llab,cex.lab,ltxt, FALSE,FALSE, # lplwgt=F,wgt=F,
        lweights,liwgt,lsyinches,
        lpllevel>1,smooth,lsmpar,smooth.iter, smooth.power=0.5,
        log="x", ylim=c(0,1.05*max(lrabs,na.rm=TRUE)),yaxs="i",
                condprobrange=condprobrange)
  } } }
## --- normal plot
  if(lpls=="qq") {
    if (is.na(lpllevel)) lpllevel <- 3*(lfgauss) # how about gamma? !!!
    for (lj in 1:lmres)
    if (lpllevel>0) {
##-       qqnorm(lstres[,lj], ylab = lstrname[lj], main="", col=colors[1])
      llr <- if(lcondq)
        lstres[,"random"]/lsigma  else  lstres[,lj]
      lxy <- qqnorm(llr, ylab = lstrname[lj], main="", type="n")
      abline(0,1,lty = lty[2], col=colors[2])
      if (lpllevel>2&lnsims>0) {
        lxx <- qnorm(ppoints(ln))
        lsimstr <- simresiduals(x, nrep=lnsims, resgen=rnorm)$simstres
        for (lr in 1:lnsims) {
          lines(lxx,sort(lsimstr[,lr]),lty=lty[4], lwd=lwd[4], col=colors[4])
        }
##        lines(lxx, sort(lstres[,lj]), col=colors[1])
      }
      li <- order(lxy$x)
      lxx <- lxy$x[li]
      lyy <- lxy$y[li]
      lines(lxx,lyy, col=colors[1])
      if (is.character(llab)) text(lxx,lyy,llab[li], col=colors[1]) else
        points(lxx,lyy, pch=if (length(lpch)==length(li)) lpch[li] else lpch,
               col=colors[1])
##-     lquart <- quantile(lstresa,c(0.25,0.75))
##-     abline(0, diff(lquart)/(2*qnorm(0.75)), lty = lty[2], col=colors[2])
      i.main(main, outer.margin=outer.margin)
      stamp(sure=FALSE) }
    if(lcondq) legend("bottomright",pch=c(15,3,NA),
                      legend=c("uncensored","simulated","for censored"))
  }
## --- leverage plot. If weight are present, use "almost unweighted h"
  if(lpls=="leverage")
  if ((!is.na(lpllevel))&&lpllevel>0 && lnnls) {
    if (diff(range(lhat,na.rm=TRUE))<0.001)
      warning(":plot.regr: all leverage elements equal, no leverage plot")
    else {
      llabh <- llab
      if (markprop>0 & markprop<1) {
        li <- order(lhat, decreasing=TRUE)[1:(ln*markprop/2)]
        llabh[li] <- llabels[li]
      }
      lhattit <- paste("leverages", if(lwgt) "(unweighted)")
      for (lj in 1:lmres) {
        i.plotlws(lhat, lstres[,lj], lhattit, lstrname[lj],
              main,outer.margin,cex.title,
              colors,lty,lwd,lpch, col, lpty, llabh,cex.lab,ltxt,
              lwsymbols,lwgt, lweights,liwgt,lsyinches,
              FALSE, ylim=streslim, ylimext=reslimext, yaxp=stresaxp,
                condprobrange=condprobrange)
  ##  line with constant Cook distance
      if (length(leverage.cooklim)>0) {
        llx <- seq(min(c(leverage.cooklim,4),na.rm=TRUE)^2*(1-ldfres/ln)/6,
                   max(lhat,na.rm=TRUE),length=50)
        llr <- outer(sqrt((1-llx)^2*(ln-ldfres)/(llx*ldfres)),
                     c(leverage.cooklim,-leverage.cooklim))
        matlines(llx, llr, lty=lty[2], lwd=lwd[2], col=colors[2])
      }
      }
    }
  }
  if(lpls=="resmatrix") { ## residual matrix for multivariate regr
    if (lmult) plmatrix(lres, pch=lab, main=main)
  }
  if(lpls=="qqmult") { ## qq plot of Mahalanobis lenghts for multivariate regr
  if ((!is.na(lpllevel))&&lpllevel>0) {
    lxx <- sqrt(qchisq(ppoints(lresmd),ncol(lres)))
    lor <- order(lresmd)
    lyy <- sqrt(lresmd[lor])
    lop <- par(mfrow=c(1,1))
    plot(lxx,lyy, xlab="sqrt(Chisq.quantiles)",type="n",
             ylab = "Mahal.oulyingness", main="", col=colors[1])
    lines(lxx,lyy)
    if (ltxt) text(lxx,lyy,llab[lor]) # else points(lxx,lyy,pch=llab[lor])
    axis(1)
    axis(2)
    abline(0,1,lty = lty[2], col=colors[2])
    i.main(main)
    stamp(sure=FALSE)
    par(lop)
  }}
} ## end lplsel
## ----------------------------------------------------------------
## plot residuals vs. explanatory variables by calling plresx
  lseq <- sequence
  if (length(lseq)==0||is.na(lseq)) lseq <- FALSE
  xvars <- FALSE
  if (is.logical(xplot)) xplot <- if(xplot) lform else NULL
  if(lxpl <- length(xplot) > 0) {
    if ((!is.formula(xplot)) && all(is.na(xplot))) lxpl <- FALSE  else {
      if (is.character(xplot))
        xplot <- formula(paste("y ~",paste(xplot,collapse="+")))
      xplot <- update(lform, xplot)
      xvars <- all.vars(xplot[-2])
    }
  } else lxpl <- FALSE
  lnnls <- !inherits(x, "nls")
  if (!lnnls)  xvars <- setdiff(xvars, names(eval(x$call$start)))
  ##  plotting limits
  lylim <- if (addcomp) ylim else reslim
  lylimfac <- if (addcomp) ylimfac else reslimfac
  lylimext <- if (addcomp) ylimext else reslimext
  if (lxpl || is.na(lseq) || lseq) {
    plresx(x, data=data, resid=lres, vars=xvars, sequence=lseq,
           se=x.se, partial.resid=TRUE, addcomp = addcomp,
           glm.restype = glm.restype, # weights = lwgts,
           wsymbols=lwsymbols, symbol.size=symbol.size,
           markprop=NULL, lab=llabna, cex.lab=cex.lab,
           mbox=mbox, jitterbinary = jitterbinary,
           smooth=x.smooth, smooth.par=lsmpar, smooth.iter=smooth.iter,
           smooth.sim=lsimres, lty=lty, lwd=lwd, colors = colors,
           pch=lpch, col=col,
           main=main, cex.title=cex.title,
           ylim = lylim, ylimfac = lylimfac, ylimext = lylimext, yaxp=resaxp,
           ...)
  }
## --- end
  par(loldpar)
  "plot.regr done"
}
## ==========================================================================
i.plotlws <- function(x,y, xlab="",ylab="",main="", outer.margin=FALSE,
  cex.title=1.2, colors=1:9, lty=1:6, lwd=1, pch=1, col=1,
  pty="p", lab="+", cex.lab=1, # cex.lab[2] is cex for points
  txt=FALSE,  plwgt=FALSE, wgt=1, weights=NULL, iwgt=NULL, syinches=0.1,
  do.smooth=TRUE, smooth=i.smooth, smooth.par=NULL, smooth.iter=10,
  smooth.power=1, # smooth.groups = NULL,
  ylim=NULL, ylimfac=3.0, ylimext=0.1, yaxp=NULL,
  reflinex=NULL, refliney=NULL, reflineyw=NULL, lnsims=0, simres=NULL,
  smooth.group = NULL, smooth.col=colors[3:4], smooth.pale=0.2,
  smooth.legend = TRUE,
  condprobrange=c(0.05,0.8), new=TRUE, axes=1:2, ...)
{
  ## Purpose:   panel for residual plots (labels, weights, smooth)
  ##   produces more than one panel for multivariate y
  ##   shows vertical bars if y is a matrix of class condquant
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## !!! simres only works with ncol(y)==1
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  5 May 2004, 06:46
  lx <- cbind(x)
  ly <- cbind(y)
  lcq <- inherits(y,"condquant")
#  if (lcq) lnsims <- 0
  lnc <- if (lcq) 1 else ncol(ly)
  lnr <- nrow(lx)
  cex.pch <- if (length(cex.lab)>1) cex.lab[2] else
     0.7*cex.lab*max(min(1,log(50)/log(lnr))^2,0.3) ## adjust symbol size to n
  if (all(pch==".")) cex.pch <- 4*cex.pch
  if(is.factor(pch)) pch <- as.character(pch)
  lxlab <- rep(xlab,length=lnc)
  lylab <- rep(ylab,length=lnc)
  if (is.null(weights)) weights <- rep(1,length(x))
  ldj <- ncol(lx)==lnc
  ## plot range
##-   if (is.logical(ylim)&&ylim) ylim <- apply(ly,2,robrange, fac=ylimfac)
##-   llimy <- length(ylim)>0 && !any(is.na(ylim))
##-   lylim <- if (llimy) {
##-     if (length(ylim)==1)  ylim <- c(-ylim,ylim)
##-     matrix(ylim,2,lnc)
##-   }  else  NULL
  lylim <- if (length(ylim)==2|length(ylim)==2*lnc)
    matrix(ylim,2,lnc) else NULL
  llimy <- length(lylim)>0
  lrf <- !is.null(reflinex)
  reflinex <- cbind(reflinex)
  refliney <- cbind(refliney)
  reflineyw <- cbind(reflineyw) # refline: width of confidence interval
  lrfyw <- length(reflineyw)>0
  ldjrflx <- ncol(reflinex)==lnc
  ldjrfly <- ncol(refliney)==lnc
## do plot(s)
  ljx <- ljrflx <- ljrfly <- 1
  for (lj in 1:lnc) {
    lxj <- lx[,ljx]
    if (lcq) { # condquant: only 1 y possible
      lyj <- y[,1]
      lcqint <- y[,"prob"]!=0
      lypl <- y
      lyplj <- if (llimy) {
        lylimj <- lylim
        plcoord(y[,1:4], range=lylimj, limext=ylimext)
      }  else  y[,1:4]
    } else { # usual case
      lyj <- ly[,lj]
      if (llimy) {
        lylimj <- lylim[,lj]
        lyplj <- plcoord(lyj, range=lylimj, limext=ylimext)
      } else lyplj <- lyj
    }
    ## plotting frame
    if (new) {
      lyrgj <- range(lyplj,na.rm=TRUE)
      plot(range(lxj,na.rm=TRUE), lyrgj, xlab = lxlab[lj], ylab = lylab[lj],
           type="n", bty="n", axes=FALSE, ...)
      if (1%in%axes) axis(1)
      ## box
      if (llimy & length(attr(lyplj,"nmod"))) {
        lusr <- par("usr")
        box(lty=3)
        ## inner box
        lines(lusr[c(1,1,2,2,1)],
              c(max(lylimj[1],lusr[3]),min(lylimj[2],lusr[4]))[c(1,2,2,1,1)],
              xpd=TRUE)
        if (2%in%axes) { ## axis labels only in inner range
          if (is.null(yaxp)) {
            lat <- pretty(lylimj, n=6, min.n=5)
            lat <- lat[lat>=lylimj[1]&lat<=lylimj[2]]
            axis(2,at=lat)
          } else axis(2,yaxp=yaxp)
        }
      } else {
        box()
        if (2%in%axes) axis(2)
      }
    }
    abline(0, 0, lty = lty[2], col=colors[2])
    lusr <- par("usr")
    ## conditional quantiles
    if(lcq) {
      ## vertical
      li <- lcqint & lypl[,"prob"]>=condprobrange[1] &
                     lypl[,"prob"]<=condprobrange[2]
      if (any(li,na.rm=TRUE))
      segments(lxj[li],lypl[li,"lowq"],lxj[li],lypl[li,"uppq"],col=colors[9])
      ## horizontal
      ldx <- diff(lusr[1:2])*0.02
      segments(lxj[lcqint]-ldx,lypl[lcqint,"median"],
               lxj[lcqint]+ldx,lypl[lcqint,"median"], col=colors[8])
    }
    if (lrf) {  ## prep reference line
      lrfxj <- reflinex[,ljrflx]
      lrfyj <- refliney[,ljrfly]
      if(length(lrfxj)==1) { # given as abline: intercept, slope
        lrfxj <- lusr[1:2]
        lrfyj <- lrfxj+lrfyj*lusr[1:2]
      }
      if (lrfyw) {
        lrfylj <- lrfyj-reflineyw[,ljrfly]*lrfxj
        lrfyuj <- lrfyj+reflineyw[,ljrfly]*lrfxj
      }
      if (llimy) { ## adjust ref lines to limited y range
        if (length(lrfxj)==2) {
          lrfxj <- seq(lusr[1],lusr[2],length=50)
          lrfyj <- reflinex[,ljrflx]+refliney[,ljrfly]*lrfxj
          if (lrfyw) {
            lrfylj <- lrfyj-reflineyw[,ljrfly]*lrfxj
            lrfyuj <- lrfyj+reflineyw[,ljrfly]*lrfxj
          }
        }
        lrfyjp <- plcoord(lrfyj, range=lylimj, limext=ylimext)
        lrfyl <- lrfyj<lrfyjp
        lrfyh <- lrfyj>lrfyjp
        lrfyr <- lrfyj==lrfyjp
        lines(lrfxj[lrfyr],lrfyj[lrfyr],lty=lty[5],col=colors[5],lwd=lwd[5])
        ## draw reference lines
        if (any(lrfyl,na.rm=TRUE))
          lines(lrfxj[lrfyl],lrfyjp[lrfyl],lty=lty[5],col=colors[5],lwd=lwd[5]/2)
        if (any(lrfyh,na.rm=TRUE))
          lines(lrfxj[lrfyh],lrfyjp[lrfyh],lty=lty[5],col=colors[5],lwd=lwd[5]/2)
        if (lrfyw) {
          lrfylj[lrfylj<lylimj[1]|lrfylj>lylimj[2]] <- NA
          lrfyuj[lrfyuj<lylimj[1]|lrfyuj>lylimj[2]] <- NA
        }
      }  else  lines(lrfxj,lrfyj,lty=lty[5],col=colors[5],lwd=lwd[5])
      if (lrfyw) {
        lines(lrfxj,lrfylj,lty=lty[6],col=colors[6],lwd=lwd[6])
        lines(lrfxj,lrfyuj,lty=lty[6],col=colors[6],lwd=lwd[6])
      }
    }
    ## smooth
    if(do.smooth) { 
      lna <- (is.na(lxj)|is.na(lyj))
      lxj[lna] <- NA
      lio <- order(lxj)[1:sum(!lna)]
      lxjo <- lxj[lio] # sorted without NA
      lwgo <- weights[lio]
      lyjo <- lyj[lio]
      lsmgrp <- rep(1,length(lio))
      lngrp <- 1
      if (length(smooth.col)==0) smooth.col <- colors[3:4] 
      lsmcol1 <- smooth.col[1]
      lsmcol2 <- smooth.col[length(smooth.col)]
      if (length(smooth.group)) { # groups
        lsmgrp <- factor(smooth.group[lio])
        lsmglab <- levels(lsmgrp)
        lsmgrp <- as.numeric(lsmgrp)
        lngrp <- length(lsmglab)
        if (NROW(smooth.col)<lngrp) smooth.col <- colors
        lsmcol1 <- rep(smooth.col, length=lngrp)
        lsmcol2 <- if(NCOL(smooth.col)>1)
                       rep(smooth.col[,2], length=lngrp)  else
        colorpale(lsmcol1, pale=smooth.pale)
      }
      if (lnsims>0) {
        simres <- simres[lio,]
        for (lr in 1:lnsims) {
          for (lgrp in 1:lngrp) {  ## smooth within groups (if >1)
            lig <- lsmgrp==lgrp
            lsms <- try(smooth(lxjo[lig], simres[lig,lr]^smooth.power,
                               weights=lwgo[lig], par=smooth.par,
                               iter=smooth.iter)^(1/smooth.power) )
            if (class(lsms)!="try-error" & length(lsms)) {
              if (llimy)  lsms[lsms<lylimj[1]|lsms>lylimj[2]] <- NA
              lines(lxjo[lig], lsms, lty=lty[4], col=lsmcol2[lgrp],lwd=lwd[4])
            }
          }
        }
      }
      for (lgrp in 1:lngrp) {  ## smooth within groups (if >1)
        lig <- lsmgrp==lgrp
        lsm <- try(smooth(lxjo[lig], lyjo[lig]^smooth.power, weights=lwgo[lig],
                      par=smooth.par, iter=smooth.iter)^(1/smooth.power))
        if (class(lsm)!="try-error" & length(lsm)) {
          if (llimy)  # lsm[lsm<lylimj[1]|lsm>lylimj[2]] <- NA
            lsm <- plcoord(lsm, range=lylimj, limext=ylimext)
          lines(lxjo[lig], lsm, lty=lty[3], col=lsmcol1[lgrp],lwd=lwd[3])
        }
    }
      if (lngrp>1) { ## legend for smooths
        if (length(smooth.legend)) {
          if (is.logical(smooth.legend))
            smooth.legend <- if (smooth.legend) "topright" else NULL
          if (length(smooth.legend))
            legend(smooth.legend, legend=lsmglab, lty=rep(lty[4],lngrp),
                   col=lsmcol1, lwd=3, bg="white")
        }
      }
    }
    ##- points
    lpclr <- rep(if (length(col)) col else colors[1], length=lnr)
    if (lcq) {
      lyplj <- lypl[,"random"]
      lpclr <- if (length(col)) rep(col,length=lnr) else
        ifelse(lcqint,colors[7],colors[1])
    }
    if (length(lyj)) {
    if (txt) {
      text(lxj,lyplj, lab, cex=cex.lab, col=lpclr)
      if (plwgt) {
        if (is.null(iwgt)) iwgt <- rep(TRUE,length(x))
        symbols(lxj[iwgt], lyplj[iwgt], circles=sqrt(weights[iwgt]),
                     inches=syinches, fg=lpclr[iwgt], add=TRUE)
      }
      else  if (any(li <- is.na(lab)|lab==""))
        points(lxj[li], lyplj[li], pch=pch, cex=cex.pch, col=lpclr[li])
    } else   points(lxj, lyplj, pch=pch, cex=cex.pch, col=lpclr)
##-       if (any(li <- is.na(lab)|lab==""))
##-         points(lxj[li], lyplj[li], pch=pch, cex=cex.pch, col=lpclr)
  }
    ljx <- ljx+ldj
    ljrflx <- ljrflx+ldjrflx
    ljrfly <- ljrfly+ldjrfly
  }
    i.main(main, cex = cex.title, outer.margin=outer.margin)
    if (new) stamp(sure=FALSE)
  NULL  ## return last yrange
}
## ==========================================================================
i.smooth <- function(x,y,weights,par=5*length(x)^log10(1/2), iter=50)
{
  if (length(x)<8) return(NULL)
  lsm <- ## ----------------------------------------------------------------
formTrsf <- function(formula, xtrsf) {
    ## recode formula to avoid evaluation of transformed variables
    ## xtrsf:  labels of transformed variables as used in data.frame
    ## formula:  --> all occurencies of elements of  xtrsf  will be
    ##   replaced by  .Xt1, .Xt2, ...
    lform <- as.character(formula)
    lf3 <- last(lform)
    lxtr <- xtrsf[io <- order(nchar(xtrsf),decreasing=TRUE)]
    lnew <- paste(".Xt",1:length(lxtr),sep="")
    lno <- lnew[io]
    for (li in seq_along(lxtr))
        lf3 <- gsub(lxtr[li], lno[li], lf3, fixed=TRUE)
    structure( as.formula(paste(if(length(lform)==3) lform[2], "~", lf3)),
              newvars=lnew)
}
  lsm <- if (u.debug())
             loess(y~x, weights=weights, span=par, iter=iter,
                   family=if (iter>0) "symmetric" else "gaussian",
                   na.action=na.exclude)  else
  try(loess(y~x, weights=weights, span=par, iter=iter,
            family=if (iter>0) "symmetric" else "gaussian",
            na.action=na.exclude), silent=TRUE)
  if (class(lsm)=="try-error") {
    lsm <- loess(y~x, weights=weights, span=0.99, iter=iter,
               family=if (iter>0) "symmetric" else "gaussian",
                 na.action=na.exclude)
    warning(":i.smooth: span was too small. Using 0.99")
  }
  fitted(lsm)
}
## i.smoothrob <- function(x,y,weights,par=3*length(x)^log10(1/2),iter=50)
## ==========================================================================
simresiduals <- function(object, nrep, resgen=NULL, glm.restype="deviance")
{
  ## Purpose:   simulate residuals according to regression model
  ##            by permuting residuals of actual model or by random numbers
  ## ----------------------------------------------------------------------
  ## Arguments:  resgen: how are residuals generated?
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 10 Aug 2008, 07:55
  if (!class(object)[1]%in%c("regr","lm","glm")) {
    warning(":simresiduals: ",
            "I can simulate only for `regr`, `lm`, or `glm` objects")
    return(NULL)
  }
  lcall <- object$call
  if ("weights"%in%names(lcall)) {
    warning(":simresiduals: I cannot simulate for weighted regression (yet)")
    ## get_all_vars  contains  weights  without parentheses -> danger!
    return(NULL)
  }
  lerror <- NULL
  ldata <- object$allvars
  if (is.null(ldata)) {
      ldata <- if (u.debug())  eval(lcall$data)  else
      try(eval(lcall$data))
    if (class(ldata)=="try-error"||is.null(dim(ldata))) {
      warning(":simresiduals: data not found -> No simulated residuals")
      return(NULL)
    }
  }
##-   if (length(object$na.action)) {
##-     attr(ldata,"na.action") <- object$na.action
##-     ldata <- na.action(ldata)
##-   }
  lres <- object$stres * object$sigma
  if (lnostres <- (length(lres)==0)||all(!is.finite(lres)))
    lres <- residuals(object)
  if (inherits(lres, "condquant"))
    lres <- structure(lres[,"random"], names=row.names(lres))
  if (length(lres)==0||all(lres==lres[1])) {
    warning(":simresiduals: no (distinct) residuals found -> No simulated residuals")
    return(NULL)
  }
  lfit <- object$fitted
  if (is.null(lfit)) lfit <- object$linear.predictors
  if (is.null(lfit)) lfit <- 0

  if (nrow(ldata)!=length(lres)) {
    li <- match(names(lres),row.names(ldata))
    if (any(is.na(li))) {
      warning(":simresiduals: data not suitable -> No simulated residuals")
      return(NULL)
    }
    ldata <- ldata[li,]
  }
        ##!!! weights
  lina <- is.na(lres)
  if (any(lina)) {
    lres <- lres[!lina]
    ldata <- ldata[!lina,]
    lfit <- rep(lfit,length=length(lres))[!lina]
  }
  if (nrow(ldata)<=2) {
    warning(":simresiduals: <=2 residuals found -> No simulated residuals")
    return(NULL)
  }
## ---
  if (length(resgen)) {
    if (!is.function(resgen)) resgen <- rnorm
    lsig <- object$sigma
    if (length(lsig)!=1) lsig <- 1  ## only standardized res useful!
    lrgen <- TRUE
  } else  lrgen <- FALSE
  ## -------
  lcall$data <- as.name("ldata")
  lform <- formula(object)
  lynm <- all.vars(lform[[2]])
  environment(lform) <- environment()
  lcall$formula <- lform
  lcall$model <- NULL
  lcall$termtable <- NULL
  lnrow <- nrow(ldata)
  lsimres <- matrix(NA,lnrow,nrep)
  lfam <- object$distrname
  if (length(lfam)==0) lfam <- object$family$family
  ## ---
  ## weibull not yet implemented
  if (lfam%in%c("gaussian","Gaussian")) {
    lcall$formula <- update(lform, paste(lynm,"~.")) ## needed for transformed y
    lsimstres <- if (lnostres) NULL else lsimres
    for (lr in 1:nrep) {
      ldata[,lynm] <- lfit + if (lrgen) resgen(lnrow)*lsig else sample(lres)
      lrs <- eval(lcall) ## update(x, formula=lfo, data=ldata)
      lrsr <- residuals(lrs)
      if (inherits(lrsr, "condquant")) lrsr <- lrsr[,"random"]
      lsimres[,lr] <- lrsr
      if (!lnostres) lsimstres[,lr] <- lrs$stres
    }
  ## ---
  }  else  {  ## glm
    ly <- object$response
    if (is.null(ly))  ly <- object$y
    lfam <- object$family$family
    if (is.null(lfam)) lfam <- ""
    switch(lfam,
           binomial = {
             if (NCOL(ly)==1 && length(unique(ly))!=2) {
               warning(":simresiduals: binomial distribution with unsuitable response.\n",
                       "No residuals simulated")
               return(list(simres=numeric(0)))
             }
             resgen <- function(x) {
               if(NCOL(ly)==1) rbinom(x, 1, lfit)   else {
                 lnbin <- ly[,1]+ly[,2]
                 ly1 <- rbinom(x, lnbin, lfit)
                 cbind(ly1,lnbin-ly1)
               }
             }
           },
           poisson = {
             resgen <- function(x) rpois(x, lfit)
           },
           {
             warning(":simresiduals: not (yet) available for this type of model.\n",
                     "No residuals simulated")
             return(list(simres=numeric(0)))
           }
           )
    lrgen <- TRUE
    lsimstres <- NULL
    for (lr in 1:nrep) {
      ldata[,lynm] <- if (lrgen) resgen(lnrow)  else  sample(lres)
      lrs <- eval(lcall) ## update(x, formula=lfo, data=ldata)
      lsimres[,lr] <- residuals(lrs, type=glm.restype)
    }
  }
  list(simres=lsimres, simstres=lsimstres)
}
## ==========================================================================
plresx <-
  function (x, data = NULL, resid=NULL, vars = NULL, sequence=FALSE,
            se = FALSE, partial.resid = TRUE, addcomp = FALSE, 
            glm.restype = "deviance", condprobrange=c(0.05,0.8),            
            weights = NULL, wsymbols=NULL, symbol.size=NULL,
            markprop=NULL, lab = NULL, cex.lab = 0.7,
            reflines = TRUE, mbox = FALSE, jitter=NULL, jitterbinary=TRUE,
            smooth = TRUE, smooth.par=NA, smooth.iter=NULL,
            smooth.sim=19, nxsmooth=51, smooth.group=NULL,
            smooth.col=NULL, smooth.pale=0.2, smooth.legend=TRUE, 
            multnrows = 0, multncols = 0, 
            lty = c(1,2,5,3,6,4,1,1), lwd=c(1,1,2,1,1.5,1,1,1),
            colors = getUserOption("colors.ra"), pch=NULL, col=NULL,
            rug = FALSE, 
            main = NULL, cex.title = NULL,  xlabs = NULL, ylabs = NULL, 
            ylim=TRUE, ylimfac=4.0, ylimext=0.1, yaxp=NULL, 
            cex = par("cex"), ask = NULL,
            ...)
{
## Purpose:        draw residuals against x variables,
##                 with component of the fit
## -------------------------------------------------------------------------
## Arguments:
##   glm.restype   type of residuals to be used for  glm  models
##   weights       weights to be used by smoother and for symbol sizes
##   lab           labels for plotting points. weights will only be
##                 used for plotting by circles if  lab  is null or
##                 for those points for which  lab  is NA.
##   cex.lab       character size to be used for plotting points
##   vars          variables to be used as horizontal axes
##   sequence      if TRUE, the sequence of the observations will be used
##                 as a variable
##   se            should confidence bands be plotted?
##   addcomp       should the  component be added to the residual?
##   smooth        a smoothing function, or TRUE
##   smooth.sim    how many simulated smoothers should be shown?
##                 can be a matrix of simulated residuals
##   data
##   ylim          plot ranges for y axes.
##                 columns correspond to variables in data
##                 if  TRUE , robust plot range will be determined using:
##     ylimfac       factor determining robust plot range
##     ylimext       extension of plotting area for showing outliers
##   colors, lty   colors and line types used for the different
##                 elements of the plot
##     [1]         observations
##     [2]         horizontal reference line
##     [3]         smooth
##     [4]         simulated smooths
##     [5]         reference lines = terms
##     [6]         confidence bands of terms
##   xlabs         labels of x variables to be used instead of the names
##                 of the data.frame data
##   wsymbols      plot points by circles according to weights (of x)
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: Jan 2002 / Jan 2004
  lnaaction <- x$na.action
  x$na.action <- NULL
## family
  lfam <- x$distrname
  if (is.null(lfam)) lfam <- x$family$family
  if (is.null(lfam) || lfam=="" || is.na(lfam)) lfam <- "gaussian"
  lfgauss <- lfam == "gaussian"
  lglm <- !lfgauss
  lpolr <- inherits(x,"polr")
  lnnls <- !inherits(x, "nls")
## residuals
  rtype <- "response"
  if (lglm) rtype <- glm.restype
  lres <- if (length(resid)) resid else residuals(x, type=rtype)
##-   if(is.na(pmatch(rtype,"condquant")))
##-     lres <- residuals(x, type=rtype) else {
##-       lres <- residuals.polr(x)
##-       lpolr <- NCOL(lres)>1
##-     }
  lcondq <- inherits(lres,"condquant")
  if (NCOL(lres)==1) lres <- cbind(lres)  ## do not destroy class attr
  lmres <- if (lcondq) 1 else ncol(lres)
  lmult <- lmres>1
  ln <- nrow(lres)
  llab <- if (length(lab)==0) rep(NA,length=ln)  else
    rep(lab,length=ln)
  ltxt <- is.character(llab) # &&any(nchar(llab)>1)
  lpty <- ifelse(ltxt,"n","p")
  liwgt <- is.na(llab)
  if (is.character(llab)) {
    liwgt <- liwgt| llab=="NA"|llab==""
    lpch <- if (length(pch)) pch else 3 # was ifelse(ln>200,".","+")
  } else lpch <- if (length(pch)) pch else 3
## weights
  lwgt <- if (length(weights)==1) as.logical(weights) else NA
  lweights <- if (length(weights)<=1) x$weights else {
    if (length(weights)>ln) {
      if(length(lnaaction)>0 && max(lnaaction)<=length(weights))
        weights[-lnaaction] } else weights
  }
  if (is.na(lwgt)) lwgt <- length(lweights)>1
  if (lwgt&& length(lweights)!=nrow(lres)) {
      warning(":plresx: no suitable weights found")
      lwgt <- FALSE
  }
  lwsymbols <- lwgt&any(liwgt)
  if (!is.null(wsymbols)) if ((lws <- wsymbols[1])&!lwsymbols)
      warning(":plresx: no weights available") else lwsymbols <- lws
## color vector
  if (length(col)>ln) if(length(lnaaction)>0 && max(lnaaction)<=length(col))
    col <- col[-lnaaction]
## data
  ldata <- if (length(data)==0) x$allvars else data
  if (is.null(ldata)) ldata <- try(eval(x$call$data), silent=TRUE)
  if (length(dim(ldata))!=2) stop("!plresx! unsuitable argument  data")
  if (nrow(ldata)!=nrow(lres)) {
    li <- match(row.names(lres),row.names(ldata))
    if (length(li)==0||any(is.na(li)))
      stop("!plresx! cannot match residuals and x`s")
    else {
      ldata <- ldata[li,]
    }
  }
## Prepare vars
  lvmod <- all.vars(formula(x)[[3]])
  if (!lnnls) lvmod <- setdiff(lvmod, names(eval(x$call$start)))
  if (match(".",lvmod,nomatch=0)>0) {
    lvmod <- names(x$allvars)[-1]
  }
  if ((is.logical(vars)&&!vars)) vars <- NULL else {
    if (length(vars)==0) {
      vars <- lvmod
      if (length(data)==0) data <- x$allvars
    } else {
      if (is.formula(vars)) vars <- all.vars(vars)
      if (!is.character(vars)) {
        warning(":plresx: invalid argument  vars ")
        return() }
      if ((lj <- match(".",vars,nomatch=0))>0) vars <- c(lvmod,vars[-lj])
      vars <- unique(vars)
    }
  }
  dvars <- names(ldata)
  if (length(vars)) {
    mv <- match(vars,dvars,nomatch=0)
    if (any(mv==0)) {
      lvm <- vars[mv==0]
      ldx <- eval(x$call$data)
      lj <- match(lvm,names(ldx),nomatch=0)
      if (any(lj>0)) {
        ldx <- ldx[,lj[lj>0],drop=FALSE]
        if (NROW(ldx)!=nrow(ldata)) {
          li <- match(row.names(ldata),row.names(ldx))
          if (any(is.na(li))) {
            warning(":plresx: Incompatible data.frames. Specify argument `data`.")
            lj[] <- 0
          } else  ldx <- ldx[li,,drop=FALSE]
        }
        if (any(lj>0)) ldata <- cbind(ldata, ldx)
      }
      if (any(lj==0)){
        warning (paste(":plresx: Variable(s) ",
                       paste(lvm[lj==0],collapse=", ")," not found"))
        vars <- setdiff(vars,lvm[lj==0])
      }
    }
    ldata <- ldata[,vars,drop=FALSE]
## convert character and variables with 2 values to factor
    for (lvn in 1:ncol(ldata)) {
      lv <- ldata[[lvn]]
      if (is.character(lv)||(jitterbinary&&length(unique(lv))<=2))
        ldata[[lvn]] <- factor(lv)
    }
  } # if (length(vars))
  sequence <- length(sequence) && (!is.na(sequence)) && sequence
##-        ( is.na(sequence) || (is.logical(sequence)&&sequence) )
  if (sequence) {
    ## is the seqence represented by any other variable?
    lseqvar <- if (length(vars)>0)
      sapply(ldata[,vars,drop=FALSE],function(x) {
           if (is.factor(x)||is.character(x)) FALSE else {
             ld <- diff(x)
             sum(ld==0)<0.1*length(x) && (all(ld<=0) | all(ld>=0)) }
           } ) else FALSE
    sequence <- !any(lseqvar)
    if (any(lseqvar)) warning(paste(":plresx: sequence represented by",
                                    paste(vars[lseqvar],collapse=", ")))
    ldata$.sequence <- 1:nrow(ldata)
    vars <- c(".sequence",vars)
    dvars <- c(".sequence",dvars)
  }
  ncd <- length(dvars)
  nvars <- length(vars)
  if (nvars==0) {
    warning(":plresx: I did not find any x variables")
    return() }
  terminmodel <- match(vars,lvmod,nomatch=0)>0
##  tin <- vars[terminmodel]
##  reference lines
  if (is.null(reflines)) reflines <- !inherits(x,"coxph")
## weight symbols
  if (lwgt) {
    lweights <- lweights/mean(lweights)
    if (length(symbol.size)==0||is.na(symbol.size)) symbol.size <- 3/ln^0.3
    lsyinches <- 0.02*symbol.size*par("pin")[1] # *max(lweights,na.rm=TRUE)
    if (lwsymbols) llab[liwgt] <- NA # ifelse(is.character(llab),"",0)
  } else  {
    lweights <- NULL
#    llab[liwgt] <- lpch
  }
  lipts <- !(lwsymbols&liwgt) ## points to be shown by  lab
## smooth
  if (is.null(smooth.iter)||is.na(smooth.iter)) smooth.iter <- {
    ldn <- x$distrname
    50 * if (is.null(ldn)) 1 else !ldn%in%c("binomial","multinomial")
  }
  if (is.logical(smooth)) smooth <- if (smooth)  i.smooth  else NULL
##-     function(x,y,weights,par,iter)
##-       fitted(loess(y~x, weights=weights, span=par, iter=iter,
##-                    na.action=na.exclude))
  ldosm <- length(smooth)>0
  lsmpar <- if (is.na(smooth.par)) 3*ln^log10(1/2)*(1+lglm) else
               smooth.par# 2
## simulated residuals
  lnsims <- smooth.sim
  if (length(lnsims)==0) lnsims <- 0
  lnsims <- if (is.logical(lnsims)&&lnsims) 19 else as.numeric(lnsims)
#  if (lcondq) lnsims <- 0
  if (lmult) lnsims <- 0
  if (length(lnsims)>1) {
    lnsims <- ncol(smooth.sim)
    if (length(lnsims)==0) {
      warning(":plresx: unsuitable argument smooth.sim")
      lnsims <- 0 } else lsimres <- smooth.sim
  } else  if (lnsims)  {
    lsimres <- simresiduals(x, lnsims)$simres
    if (length(lsimres)==0) lnsims <- 0
  }
## graphical elements
  lty <- rep(c(lty,1:6),length=6)
  lwd <- rep(c(lwd,1),length=6)
  colors <- if (length(colors)==0)
    c("black","gray","blue","cyan","red","magenta","darkgreen",
      "green","lightgreen")  else
    rep(colors,length=9)
  if (length(jitter)==0) jitter <- 0.3*(1-10^(-0.01*pmax(0,ln-10)))
  if (is.na(smooth.par)) smooth.par <- min(0.99, 5*ln^log10(1/2)*(1+lglm)) # 2
## ask
  if (length(ask)==0)
    ask <- interactive() && (prod(par("mfcol")) < nvars) &&
           (.Device != "postscript")
## type
  addcomp <- as.logical(addcomp)
  if (is.na(addcomp)) addcomp <- FALSE
  lty.term <- lty[3]
  lwd.term <- lwd[3]
  if (is.na(lty.term)) lty.term <- if (addcomp) 1 else 3
## x and y axes
  if (length(xlabs)>0)
    if (length(xlabs)!=length(dvars)) {
      warning("argument  xlabs  has wrong length")
      xlabs <- NULL }
  if (length(xlabs)==0) xlabs <- dvars
  xlabs[xlabs==".sequence"] <- "sequence"
  names(xlabs) <- dvars
  if (length(ylabs)==0)
      ylabs <- if (addcomp) {
        ifelse(terminmodel, paste("Partial for", dvars), "Residuals")
      } else "Residuals"
  ylabs <- rep(ylabs, length=ncd)
  names(ylabs) <- dvars
## plot range preparation
  if (is.logical(ylim)&&!is.na(ylim)) ylim <- if (ylim) {
      if(lcondq) robrange(c(lres[,1:3]), fac=ylimfac) else
      apply(lres,2,robrange, fac=ylimfac)
    } else  NULL
  lylim <- if (length(ylim)) {
    if (any(dim(cbind(ylim))!=c(2,lmres))) {
      warning(":plot.regr: unsuitable argument  ylim ")
      NULL
    } else cbind(ylim)
  } else NULL
##-   llimy <- (!is.logical(ylim))&&length(ylim)>0
##-   lylim <- NULL
##-   if (llimy) {
##-     if (length(ylim)==1)  ylim <- c(-ylim,ylim)
##-     lylim <- matrix(ylim,2,lmres) ## not same as plmatrix
##-   } else if(llimyrob & !addcomp)
##-     lylim <- apply(lres, 2, robrange, fac=ylimfac)
##-   if (length(lylim))
##-     dimnames(lylim) <- list(c("min","max"),colnames(lres))
##-   if (lmres==1) lylim <- c(lylim)
  llimy <- length(lylim)>0
## plot title
  lmain <- if (length(main)==0||(is.logical(main)&&main))
    paste(x$call["formula"],collapse="") else ""
  if (is.character(main)) lmain <- main
  lonetitle <- length(lmain)==1
  if (length(lmain))  tit(lmain) <- tit(x)
  loma <- par("oma")
  outer.margin <- lonetitle&&loma[3]>0
  if (is.null(cex.title)) cex.title <- max(0.5, min(1.2,
      par("mfg")[4]*par("pin")[1]/(par("cin")[1]*nchar(main))))
## machinery
##  ldfres <- x$df.residual
  ldfres <- x$df
  if(length(ldfres)>1) ldfres <- ldfres[2]
  qnt <- if (length(ldfres)==0) 2 else qt(0.975,ldfres)
  is.fac <- sapply(ldata[,vars,drop=FALSE], is.factor)
##-   if (sum(terminmodel)>0) {
##-     comps <- fitcomp(x, ldt, vars=tin, xfromdata=addcomp, se=se)
##-     if (addcomp) {
##-       if (length(ldata)==0&&(ln!=nrow(ldata))) # length(res)
##-         stop("plresx: BUG: number of residuals != nrow(comps)")
##-       xcomp <- ldata[,tin]
##-       lcmp <- comps$comp
##-     } else {
##-       xcomp <- comps$x
##-       lcmp <- -comps$comp
##-     }
##-   }
##-   environment(x$call$formula) <- environment()
##-   x$call$data <- as.name("ldata")
  if (sum(terminmodel)>0 && reflines &&lnnls) {
##-     lcmp <- try(fitcomp(x, xfromdata=FALSE, se=se))
    lcmp <- fitcomp(x, xfromdata=FALSE, se=se)
##-     if (class(lcmp)=="try-error") {
##-       warning(":plresx: fitcomp did not work. no reference lines")
##-       terminmodel[] <- FALSE
##-     } else {
      lcompx <- lcmp$x
      lcompy <- if (addcomp) lcmp$comp else -lcmp$comp
      lcompse <- lcmp$se
      if (addcomp) {
##-         if (length(ldata)==0&&(ln!=nrow(ldata))) # length(res)
##-           stop("!plresx! BUG: number of residuals != nrow(comps)")
        lcompdt <- fitcomp(x, ldata, vars=vars[terminmodel], xfromdata=TRUE)$comp
      }
##-     }
  }
## ------------------------------------------------------------------
  if (inherits(x,"mlm")) {
##-     ff <- function(v2, v1, j2, j1, pch, clr, clrsmooth) {
##-             points(v2,v1, col=clr, pch=pch)
##-             lines(lowess(v2,v1), col=clrsmooth) }
    lpanel <- function(xx, yy, jx, jy, pch, clr) {
      lcmpx <- lcmpy <- NULL
      ltin <- terminmodel[jx]
      lvx <- vars[jx]
      lcnt <- !is.fac[lvx]
      if (ltin) {
        lcmpy <- lcompy[,lvx,jy]
        if (lcnt) lcmpx <- lcompx[,lvx]
      }
      i.plotlws(xx,yy, "","","", TRUE, cex.title, colors, lty, lwd, lpch, col,
                lpty, llab,cex.lab,ltxt, lwsymbols,lwgt,lweights,liwgt,lsyinches,
                ldosm&lcnt,smooth,lsmpar,smooth.iter,
##-                 ylim=lylim,
                yaxp=yaxp,
                lnsims=lnsims, simres=lsimres, new=FALSE,
                reflinex=lcmpx, refliney=lcmpy,
                smooth.group=smooth.group, smooth.col=smooth.col,
                smooth.pale=smooth.pale)
    }
    plmatrix(ldata[,vars,drop=FALSE],lres, panel=lpanel, pch=llab, range.=lylim,
             reference=FALSE,
             nrows=multnrows, ncols=multncols, main=main) # clrsmooth=colors[3]
    return()
  }
## --- loop ---
  for (lj in 1:nvars) {
    lv <- vars[lj]
    lcmpj <- terminmodel[lj] && reflines&&lnnls
    lci <- if (lcmpj) lcompy[, lv] else 0 ## !!!
    rr <- lres
    if (partial.resid)
      if (addcomp && lcmpj) rr <- rr+lcompdt[, lv]
## - plot range
##-     if (llimy) rr <- plcoord(rr, range=lylim, limext=ylimext)  else
##-       if (llimyrob) rr <- plcoord(rr, limfac=ylimfac, limext=ylimext)
##-     ylims <- if (partial.resid) range(rr, na.rm = TRUE) else
##-                 range(lci, na.rm = TRUE)
##-     if (rug)  ylims[1] <- ylims[1] - 0.07 * diff(ylims)
## ---
    if (is.fac[lv]) { # ---
## factors
      ff <- factor(ldata[, lv])
      ll <- levels(ff)
      lnl <- length(ll)
      if (mbox)
          plmboxes(rr~ff, data=ldata, xlab = xlabs[lv], ylab = ylabs[lv],
                         refline=0) else {
      xx <- as.numeric(ff)+runif(ln,-jitter,jitter)
      xlims <- c(0.5,lnl+0.5)
      i.plotlws(xx, rr, xlab = xlabs[lv], ylab = ylabs[lv],
        "", outer.margin, cex.title,
        colors, lty, lwd, lpch, col, lpty, llab, cex.lab, ltxt,
        lwsymbols, lwgt, lweights, liwgt, lsyinches,
        FALSE, smooth, smooth.par, smooth.iter, 1,
        lylim, ylimfac, ylimext,
        reflinex=NULL, lnsims=0, axes=2,
                condprobrange=condprobrange)
      axis(1,labels=ll,at=1:lnl)
      axis(2)
      box(lty=3)
  }
      ## -
      if (lcmpj) {
        lx <- seq(along = ll)
##-         ww <- if (addcomp) match(ll,as.character(ff)) else 1:lnl
        lcil <- lci[1:lnl]
        if (!is.factor(lcompx[,lv]))  # binary (non)factor: fitcomp did not know
          lcil <- lci[c(1,last(which(!is.na(lcompx[,lv]))))]
        if (llimy) {
          lcilp <- plcoord(lcil, lylim, ylimfac, ylimext)
          if (attr(lcilp,"nmod")) {
            liout <- lcilp!=lcil
            segments(lx[liout]-0.4, lcilp[liout], lx[liout]+0.4, lcilp[liout],
                 lwd = lwd.term/2, lty=lty.term, col=colors[5])
            lcil[liout] <- NA
          }
        }
        segments(lx-0.4, lcil, lx+0.4, lcil,
                 lwd = lwd.term, lty=lty.term, col=colors[5])
        if (se) {
          wid <- qnt * lcompse[1:lnl, lv]
          lines(c(rbind(lx-0.1,lx+0.1,lx+0.1,lx-0.1,lx-0.1,NA)),
                c(rbind(lcil-wid,lcil-wid,lcil+wid,lcil+wid,lcil-wid,NA)),
                lty = lty[6], col = colors[6], lwd=lwd[6])
        }
      }
    } else { # ---
## --- continuous explanatory variable
      lrefx <- NULL
      lrefyw <- NULL
      if (lcmpj) {
        lrefx <- lcompx[,lv]
        if (se) lrefyw <- qnt*lcompse[,lv]
    }
      i.plotlws(as.numeric(ldata[, lv]), rr, xlab = xlabs[lv], ylab = ylabs[lv],
        "", outer.margin, cex.title,
        colors, lty, lwd, lpch, col, lpty, llab, cex.lab, ltxt,
        lwsymbols, lwgt, lweights, liwgt, lsyinches,
        ldosm, smooth, smooth.par, smooth.iter, smooth.power=1,
        ylim=if(llimy) lylim else NULL, ylimfac, ylimext, yaxp,
        lrefx, lci, lrefyw, lnsims, lsimres,
        smooth.group=smooth.group, smooth.col=smooth.col,
        smooth.pale=smooth.pale, smooth.legend=smooth.legend,
        condprobrange=condprobrange)
    }
##    if (llimy|llimyrob) abline(h=attr(rr,"range"),lty=lty[2],col=colors[2])
##-     if (rug) {
##-         n <- length(xx)
##-         lines(rep(jitter(xx), rep(3, n)), rep(ylims[1] +
##-             c(0, 0.05, NA) * diff(ylims), n))
##- ##-         if (partial.resid)
##- ##-             lines(rep(xlims[1] + c(0, 0.05, NA) * diff(xlims),
##- ##-               n), rep(rr, rep(3, n))
##-     }
    i.main(if (lonetitle) lmain else lmain[lj], cex = cex.title,
           outer.margin = outer.margin)
    stamp(sure=FALSE)
  }
  invisible(nvars)
}
## =======================================================================
## =======================================================================
plmatrix <-
function(x, y=NULL, data=NULL, panel=l.panel,
         nrows=0, ncols=nrows, save=TRUE, robrange.=FALSE, range.=NULL,
         pch=NULL, col=1, reference=0, ltyref=3,
         log="", xaxs="r", yaxs="r", xaxmar=NULL, yaxmar=NULL,
         vnames=NULL, main="", cex=NA, cex.points=NA, cex.lab=0.7, cex.text=1.3,
         cex.title=1, bty="o", oma=NULL, mar=rep(0.2,4), keepmf=FALSE,
         axes=TRUE, ...)
{
## Purpose:    pairs  with different plotting characters, marks and/or colors
##             showing submatrices of the full scatterplot matrix
##             possibly on several pages
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: 23 Jul 93; minor bug-fix+comments:
  ## M.Maechler
  l.panel <- function(x,y,indx,indy,pch=1,col=1,cex=cex.points,...) {
    if (is.character(pch)) text(x,y,pch,col=col,cex=cex) else
    points(x,y,pch=pch,col=col,cex=cex,...)
  }
  oldpar <- par(c("mfrow","mar","cex","mgp"))
  lmfg <- par("mfg")
  on.exit(if (!keepmf) par(oldpar))
##---------------------- preparations --------------------------
## data
  if (is.formula(x))  {
    if (length(x)==2)
    x <- model.frame(x,data, na.action=NULL)  else {
      ld <- model.frame(x[c(1,3)],data, na.action=NULL)
      ld <- cbind(ld, model.frame(x[1:2],data, na.action=NULL))
      x <- ld
    }
  }
  if (is.data.frame(x)) {
    for (jj in 1:length(x)) x[[jj]] <- as.numeric(x[[jj]])
    x <- as.matrix(x)
  } else x <- cbind(x)
##  stop("!plmatrix! first argument must either be a formula or a data.frame or matrix")
  nv1 <- dim(x)[2]
  lv1 <- lv2 <- 0
  if (is.null(y)) {
    ldata <- x
    if (save) { nv1 <- nv1-1; lv2 <- 1 }
    nv2 <- nv1
  } else { # cbind y to data for easier preparations
    save <- FALSE
    if (is.formula(y))  {
      ld <- model.frame(x[c(1,3)],data, na.action=NULL)
    if (length(x)>2)
      ld <- cbind(ld, model.frame(x[1:2],data, na.action=NULL))
    x <- ld
  }
    if (is.formula(y)) {
      if (length(y)==2)
        y <- model.frame(y,data, na.action=NULL)  else {
          ld <- model.frame(y[c(1,3)],data, na.action=NULL)
          ld <- cbind(ld, model.frame(y[1:2],data, na.action=NULL))
          y <- ld
        }
    }
    if (is.data.frame(y)) {
      for (jj in 1:length(y)) y[[jj]] <- as.numeric(y[[jj]])
      y <- as.matrix(y)
    }
    ldata <- cbind(x, as.matrix(y))
    nv2 <- ncol(ldata)-nv1 ; lv2 <- nv1 }
  nvv <- ncol(ldata)
  tnr <- nrow(ldata)
## variable labels
  lvn <- dimnames(ldata)[[2]]
  if (is.null(lvn)) lvn <- paste("V",1:nvv)
  lvnm <- lvn
  if (!is.null(vnames)) {
    vnames <- rep(vnames,length=nvv)
    lvnm[!is.na(vnames)] <- vnames[!is.na(vnames)]
  }
  vnames <- lvnm
## plotting characters
  if (length(pch)==0) pch <- 1
  if (length(pch)>tnr) pch <- pch[1:tnr]
## range
  rg <- matrix(nrow=2,ncol=nvv,dimnames=list(c("min","max"),lvn))
  if(is.matrix(range.)) {
    if (is.null(colnames(range.))) {
      if (ncol(range.)==ncol(rg)) rg[,] <- range.  else
      warning("argument  range.  not suitable. ignored")
    } else {
      lj <- match(colnames(range.),lvn)
      if (any(is.na(lj))) {
        warning("variables", colnames(range.)[is.na(lj)],"not found")
        if (any(!is.na(lj))) rg[,lj[!is.na(lj)]] <- range.[,!is.na(lj)]
      }
    }
  }
  else
    if (length(range.)==2&&is.numeric(range.)) rg[,] <- matrix(range.,2,nvv)

  lna <- apply(is.na(rg),2, any)
  if (any(lna))
    rg[,lna] <- apply(ldata[,lna,drop=FALSE],2,
      if(robrange.) robrange else range, na.rm=TRUE, finite=TRUE)
  colnames(rg) <- lvn
## reference lines
  tjref <- (length(reference)>0)&&!(is.logical(reference)&&!reference)
  if (tjref) {
    if(length(reference)==1) lref <- rep(reference,length=nvv) else {
      lref <- rep(NA,nvv)
      lref[match(names(reference),lvn)] <- reference
    }
    names(lref) <- lvn
  }
## plot
  jmain <- !is.null(main)&&main!=""
  lpin <- par("pin")
  lnm <- if (lpin[1]>lpin[2]) {
    if (nv1==6 && nv2==6) c(6,6) else c(5,6) } else c(8,5)
  if (is.na(nrows)||nrows<1) nrows <- ceiling(nv1/((nv1-1)%/%lnm[1]+1))
  if (is.na(ncols)||ncols<1) ncols <- ceiling(nv2/((nv2-1)%/%lnm[2]+1))
  if (is.null(xaxmar)) xaxmar <- 1+(nv1*nv2>1)
  if (any(is.na(xaxmar))) xaxmar <- 1+(nv1*nv2>1)
  xaxmar <- ifelse(xaxmar>1,3,1)
  if (is.null(yaxmar)) yaxmar <- 2+(nv1*nv2>1)
  if (any(is.na(yaxmar))) yaxmar <- 2+(nv1*nv2>1)
  yaxmar <- ifelse(yaxmar>2,4,2)
  if (length(oma)!=4)
    oma <- c(2+(xaxmar==1), 2+(yaxmar==2),
             1.5+(xaxmar==3)+cex.title*2*jmain,
             2+(yaxmar==4))
#    oma <- 2 + c(0,0,!is.null(main)&&main!="",1)
  if (!keepmf) par(mfrow=c(nrows,ncols))
##-   if (!is.na(cex)) par(cex=cex)
##-   cex <- par("cex")
##-   cexl <- cex*cexlab
##-   cext <- cex*cextext
  par(oma=oma*cex.lab, mar=mar, mgp=cex.lab*c(1,0.5,0))
  if (keepmf) par(mfg=lmfg, new=FALSE)
  if (!is.na(cex)) cex.points <- cex
  if (is.na(cex.points)) cex.points <- max(0.2,min(1,1.5-0.2*log(tnr)))
##
  ## log
  if (length(grep("x",log))>0) ldata[ldata[,1:nv1]<=0,1:nv1] <- NA
  if (length(grep("y",log))>0) ldata[ldata[,lv2+1:nv2]<=0,lv2+1:nv2] <- NA
  npgr <- ceiling(nv2/nrows)
  npgc <- ceiling(nv1/ncols)
##----------------- plots ----------------------------
  for (ipgr in 1:npgr) {
    lr <- (ipgr-1)*nrows
  for (ipgc in 1:npgc) {
    lc <- (ipgc-1)*ncols
    if (save&&((lr+nrows)<=lc)) break
  for (jr in 1:nrows) { #-- plot row [j]
    jd2 <- lr+jr
    j2 <- lv2 + jd2
    if (jd2<=nv2)  v2 <- ldata[,j2]
    for (jc in 1:ncols) { #-- plot column  [j2-lv2] = 1:nv2
      jd1 <- lc+jc
      j1 <- lv1 + jd1
    if (jd2<=nv2 & jd1<=nv1) {
      v1 <- ldata[,j1]
      plot(v1,v2, type="n", xlab="", ylab="", axes=FALSE,
           xlim <- rg[,j1], ylim <- rg[,j2],
           xaxs=xaxs, yaxs=yaxs, log=log, cex=cex.points)
      usr <- par("usr")
      if (axes) {
        if ((jr==nrows||jd2==nv2)) {
          if (xaxmar==1) axis(1)
          mtext(vnames[j1], side=1, line=(0.5+1.2*(xaxmar==1))*cex.lab,
                cex=cex.lab, at=mean(usr[1:2]))
        }
        if (jc==1) {
          if (yaxmar==2) axis(2)
          mtext(vnames[j2], side=2, line=(0.5+1.2*(yaxmar==2))*cex.lab,
                cex=cex.lab, at=mean(usr[3:4]))
        }
        if (jr==1&&xaxmar==3) axis(3,xpd=TRUE)
        if (jc==ncols||jd1==nv1) if (yaxmar==4) axis(4,xpd=TRUE)
      }
      box(bty=bty)
      if (any(v1!=v2,na.rm=TRUE)) { # not diagonal
        panel(v1,v2,jd1,jd2, pch, col, ...)
        if (tjref) abline(h=lref[j1],v=lref[j2],lty=ltyref)
      }
      else { uu <- par("usr") # diagonal: print variable name
             text(mean(uu[1:2]),mean(uu[3:4]), vnames[j1], cex=cex.text) }
    }
      else frame()
    }
  }
  if (jmain) mtext(main,3,oma[3]*0.9-2*cex.title,outer=TRUE,cex=cex.title)
##-   stamp(sure=FALSE,line=par("mgp")[1]+0.5)
  stamp(sure=FALSE,line=oma[4]-1.8) # ??? why does it need so much space?
  }}
  "plmatrix: done"
}

## ====================================================================
plmbox <- function(x, at=0, probs=NULL, outliers=TRUE, na.pos=NULL,
  width=1, wfac=NULL, h0=NULL, adj=0.5, extquant=TRUE, 
  ilim=NULL, ilimext=0.05, widthfac=c(max=2, med=1.3, medmin=0.3, outl=NA),
  colors=c(box="lightblue2",med="blue",na="gray90"))
{
  ## Purpose:   multi-boxplot
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 14 Dec 2013, 10:09
##-   if (is.null(col))
##-     col <- "gray70"
  f.box <- function(wid, quant, col) {
    if (wid>lwmax) {
      polygon(at+lwmax*lpos, quant, col="black")
      if (is.na(col)||col==0) col="white"
      polygon(at+lpos*lwmax^2/wid, quant, col=col)
    }  else 
    if(wid>0) polygon(at+wid*lpos, quant, col=col)
}
  
  if (length(x)==0) {
    warning(":plmbox: no data")
    return()
  }
  stopifnot(length(width)==1,length(wfac)<=1)
  if (is.null(probs))
      probs <- if (sum(!is.na(x))<20) c(0.1,0.5,1)/2 else
               c(0.05,0.1,0.25,0.50,0.75,1)/2
  lprobs <- if (all(probs<=0.5))  c(probs,1-probs)  else c(probs)
  lprobs <- sort(unique(lprobs))
  colors <- as.list(colors)
  box.col <- colors[["box"]]
  if (length(box.col)==1)
    box.col <- ifelse(0.25<=last(lprobs,-1) & lprobs[-1]<=0.75, box.col, NA)
  ## values for degenerate case
  lfac <- if (is.null(wfac)) width*2*IQR(x, na.rm=TRUE) else wfac*length(x)
                                        # was mad/dnorm(0)
  lq <- NULL
  lmed <- median(x, na.rm=TRUE)
  lwmed <- width
  lhtot <- diff(range(x, na.rm=TRUE))
  if (lhtot>0) { ## non-degenerate
  if (is.null(h0)) h0 <- lhtot*0.02
  lrg <- range(x, na.rm=TRUE)
  lqy <- lq <- quinterpol(x, probs=lprobs, extend=extquant)
  if (length(ilim)) {
      lrg <- plcoord(lrg, range=ilim, limext=ilimext)
      lqy <- plcoord(lqy, range=ilim, limext=ilimext)
  }
  loutl <- x[x<min(lq)|x>max(lq)]
  ## ---
  lwid <- lfac*diff(lprobs)/pmax(diff(lq),h0)
  lwmax <- widthfac["max"]*lfac*0.5/IQR(x, na.rm=TRUE)
  lwmed <- max(widthfac["med"]*min(lwmax,max(nainf.exclude(lwid))),
               widthfac["medmin"],na.rm=TRUE)
  lpos <- c(-adj,-adj,1-adj,1-adj)
  lwoutl <- widthfac["outl"]
  if (is.na(lwoutl)) lwoutl <- 0.1*lwmax
  ## ---
  BR()
  for (li in 1:(length(lprobs)-1)) 
      f.box(lwid[li], lqy[li+c(0,1,1,0)], box.col[li])
  if (!is.null(na.pos)) {
      lmna <- mean(is.na(x))
      if (lmna) {
          ldna <- diff(na.pos)
          if (length(ldna)==0 || is.na(ldna) || ldna==0)
              stop("!plmbox! argument 'na.pos' not suitable")
          lwidna <- lfac*lmna/abs(ldna)
          f.box(lwidna, na.pos[c(1,2,2,1)], colors[["na"]])
      }
  }
  lines(c(at,at), # +linepos*0.01*diff(par("usr")[1:2])*(0.5-adj),
        lrg, lwd=2)
  if (outliers&&length(loutl)) {
      lat <- rep(at,length(loutl))
      segments(lat-lwoutl*adj, loutl, lat+lwoutl*(1-adj), loutl)
  }
}
  ## median
  lines(at+lwmed*c(-adj,1-adj), rep(lmed,2), col=colors[["med"]], lwd=3)
##
  invisible(structure(lfac/length(x), attributes=list(q=lq,width=lwid)))
}
## ====================================================================
plmboxes <- function(formula, data, width=1, at=NULL, 
    probs=NULL, outliers=TRUE, na=FALSE, 
    refline=NULL, add=FALSE, ilim=NULL, ilimext=0.05,
    xlim=NULL, ylim=NULL, axes=TRUE, xlab=NULL, ylab=NULL, 
    labelsvert=FALSE, mar=NULL,
    widthfac=c(max=2, med=1.3, medmin=0.3, outl=NA, sep=0.003),
    colors=c(box="lightblue2",med="blue",na="gray90",refline="magenta"),...)
{
  ## Purpose:    multibox plot
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 14 Dec 2013, 23:38
##  widthfac <- 1
  f.ylim <- function(ylm, ext)
    c((1+ext)*ylm[1]-ext*ylm[2], (1+ext)*ylm[2]-ext*ylm[1])

  formula <- as.formula(formula)
  if (length(formula)<3) stop("!plmboxes! formula must have left hand side")  
  ## widths
  lwfac <- c(max=2, med=1.3, medmin=0.3, outl=NA, sep=0.003)
  if (is.null(names(widthfac)))
      names(widthfac) <- names(lwfac)[1:length(widthfac)]
  if (any(names(widthfac)%nin%names(lwfac)))
    warning(":plmboxes: argument 'widthfac' has unsuitable names") else
    lwfac[names(widthfac)] <- widthfac
  ## colors
  lcol <- list(box="lightblue",med="blue",na="gray90",refline="magenta")
  if (is.null(names(colors))) names(colors) <- names(lcol)[1:length(colors)]
  colors <- as.list(colors)
  if (any(names(colors)%nin%names(lcol)))
    warning(":plmboxes: argument 'colors' has unsuitable names") else
    lcol[names(colors)] <- colors
  ## data
  if (length(dim(data))!=2||nrow(data)==0)
    stop("!plmboxes! Argument 'data' has dimension   ",
         paste(dim(data),colapse=" ")) 
  ldt <- model.frame(formula, data, na.action=na.pass)
  ly <- ldt[,1]
  ##
  if (length(formula[[3]])>1 && as.character(formula[[3]][[2]])=="1") 
      ldt <- data.frame(ldt[,1],0,ldt[,2])
  ## preliminary 
  lx <- ldt[,2] <- factor(ldt[,2]) # unused levels are dropped
  lxx <- ldt[,-1]
  llr <- ncol(ldt)>2
  llist <- split(ly,lxx)
  llev <- levels(lx)
  lng <- length(llev)
  lnn <- sapply(llist,length)
  lsd <- mean(sapply(llist,mad,na.rm=TRUE),na.rm=TRUE)
  width <- rep(width, length=lng)
  lfac <- width*lsd/(max(lnn)*(1+llr))
  ## labels
  if (is.null(xlab)||is.na(xlab)) {
      xlab <- as.character(formula[[3]])
      if (length(xlab)>1) xlab <- xlab[2]
      if (xlab=="1") xlab <- ""
  }
  if (is.null(ylab)) ylab <- as.character(formula[[2]])
  ## position
  if (is.null(at)) at <- 1:lng else
    if (length(at)!=lng) {
      warning(":plmboxes: 'x' has wrong length")
      at <- at[1:lng] ## may produce NAs
  }
  ## probabilities
  if (is.null(probs))
      probs <- if (sum(!is.na(ly))/(lng*(1+llr))<20) c(0.1,0.5,1)/2 else
               c(0.05,0.1,0.25,0.50,0.75,1)/2
  ## graphics
  if (is.null(na)||is.na(na)||(is.logical(na)&&!na)) na.pos <- NULL else 
      if (is.logical(na))
          na.pos <- c(min(ly, na.rm=TRUE)*(1-0.3)-0.3*max(ly, na.rm=TRUE))
  if (length(na.pos)==1)
      na.pos <- na.pos+ 0.03*diff(range(ly,na.rm=TRUE))*c(-1,1)
  if (!add) {
  if (is.null(mar)) mar <-
      c(ifelse(labelsvert, min(7,1+1.1*max(nchar(llev))), 4), 4,4,1)
  oldpar <- par(mar=mar)
  if(is.null(xlim)) xlim <- 
    range(at, na.rm=TRUE)+ max(width[c(1,length(width))])*c(-1,1)*0.5
  lrg <- range(ly,na.pos, na.rm=TRUE)
  if (!is.null(ilim)) if (length(ilim)!=2 || ilim[1]>=ilim[2]) {
      warning(":plmboxes: unsuitable argument 'ilim'")
      ilim <- NULL
  }
  if (is.null(ilim)) ilim <- lrg
  ljlim <- ilim[1]<lrg[1] | ilim[2]<lrg[2] ## inner range is actif
  if(is.null(ylim)) ylim <- f.ylim(ilim,ilimext)
  ## ---------------------------------
  plot(xlim, ylim, type="n", axes=FALSE, xlab="", ylab=ylab, mar=mar, ...)
  lusr <- par("usr")
  if (axes) {
    axis(1,at=at,labels=llev,las=1+2*labelsvert)
    lat <- pretty(f.ylim(lrg, ilimext)) #, n=7,n.min=5
    if(!is.null(na.pos)) {
        lat <- lat[lat>max(na.pos)]
        mtext("NA",2,1,at=mean(na.pos),las=1)
    }
    axis(2, at=lat)
    if (ljlim) { ## inner and outer box
        box(lty=3)
        lines(lusr[c(1,2,2,1,1)],ilim[c(1,1,2,2,1)])
    } else box()
  }
  mtext(xlab,1,par("mar")[1]-1)
} # if (!add)
  ## ---
  if (!is.null(refline))
      abline(h=refline, col=lcol[["refline"]], lty=3, lwd=1.5)
  ## ---
  lusrd <- diff(par("usr")[1:2])
  lsep <- lwfac["sep"]*llr*lusrd
  lwoutl <- lwfac["outl"]
  if (is.na(lwoutl)) {
      lwoutl <- 0.05*lusrd
      lwfac["outl"] <- lwoutl/lng
  }
  if (llr) lwfac[c("medmin","outl")] <- lwfac[c("medmin","outl")] /2
  ## ------------
  for (li in 1:lng) {
    if (is.na(at[li])) next
    if (length(lli <- llist[[li]])) 
    plmbox(lli,at[li]-lsep, probs=probs, outliers=outliers, wfac=lfac[li],
           adj=0.5*(1+llr), na.pos=na.pos, extquant=TRUE,
           ilim=if(ljlim) ilim, ilimext=ilimext, 
           widthfac=lwfac, colors=lcol)
    if (llr) ## second half of asymmetrix  mbox
      if (length(llir <- llist[[li+lng]]))
        plmbox(llir,at[li]+lsep,probs=probs, outliers=outliers, wfac=lfac[li],
               adj=0, na.pos=na.pos, extquant=TRUE,
               ilim=ilim, ilimext=ilimext, 
               widthfac=lwfac, colors=lcol)
  }
  par(oldpar)
  invisible(at)
}
## ===========================================================================
plres2x <-
  function(formula=NULL, reg=NULL, data=reg, restricted=NULL, size = 0,
  slwd = 1, scol = 2, xlab = NULL, ylab= NULL, xlim=NULL, ylim=NULL,
  main = NULL, cex.title= NULL, ...)
{
## Purpose:  plot residuals vs. two x`s
## Author:   ARu , Date:  11/Jun/91
## Aenderungen: MMae, 30/Jan/92, Dez.94
## --------------------------------------------------------------------------
## Arguments:
##   formula    z~x+y, where
##              x, y  coordinates of points given by two vector arguments.
##              z     gives orientation (by sign)
##                   and size (by absolute value) of symbol.
##   reg        regression results
##   data       data
##              you must specify either  reg  or  data
##   restricted absolute value which truncates the size.
##              The corresponding symbols are marked by stars.
##   size       the symbols are scaled so that "size" is the size of
##              the largest symbol in cm.
##   main       main title, defaults to the formula
##   ...        additional arguments for the S-function `plot`
## the function currently only plots  z  for the first two terms of the
## right hand side of  formula
## --------------------------------------------------------------------------
  lform <- as.formula(formula)
  if (length(reg)==0) {
    if (length(data)==0) stop("either  reg  or  data  must be specified")
    ldata <- data
    if (length(lform)<3)
      stop ("left hand side of formula is missing. Did you mean to use reg results?")
  }  else
  if (inherits(reg,"lm")) {
    ldata <- eval(reg$call$data)
##-     ldata <- eval(parse(text=as.character(reg$call[3])))
    lftext <- deparse(formula(reg))
    if (length(formula)==0) lform <- formula(reg)[c(1,3)]
    if (length(lform)<3) {
      lform <- update.formula(lform,residuals~.)
      lrs <- resid(reg)
      if (length(lrs)!=nrow(ldata)) {
          ldata <- ldata[names(lrs),]
          if (nrow(ldata)!=length(lrs)) stop("!plres2x! residuals and data incompatible")
      }
      ldata <- data.frame(ldata,residuals=resid(reg))
    }
  } else stop("!plres2x! unsuitable argument reg")
  lftext <- deparse(lform)
  if (length(main)==0) main <- lftext
  if (is.logical(main)) main <- if (main) lftext else ""
  main <- as.character(main)
  if (length(cex.title)==0) cex.title <- max(0.5, min(1.2,
      par("mfg")[4]*par("pin")[1]/(par("cin")[1]*nchar(main))))
  if (!is.data.frame(ldata)) {
    if(is.matrix(data)) ldata <- as.data.frame(data) else
      stop("data is not a data.frame") }
  ld <- model.frame(lform,ldata)
  ld <- nainf.exclude(ld)
  z <- ld[,1]
  x <- as.numeric(ld[,2])
  y <- as.numeric(ld[,3])
  if (length(xlim)==0) xlim <- range(x)
  if (length(ylim)==0) ylim <- range(y)
##--- restrict z values: ---
  if(length(restricted)==0)   restr <- FALSE else {
    restr <- abs(z) > restricted
    z <- pmin( pmax( z, -restricted), restricted) }
## size
  if (is.null(size)||is.na(size)||size<=0)  size <- 5/log10(length(x))
  lpin <- par("pin")
  fx <- (size * diff(xlim))/100
  fy <- fx/diff(xlim)*diff(ylim)/lpin[2]*lpin[1]
##--
  if (length(xlab)==0) xlab <- names(ld)[2]
  if (length(ylab)==0) ylab <- names(ld)[3]
  plot(x, y, xlim = xlim + c(-1,1)* fx, ylim = ylim + c(-1,1)* fy, pch = ".",
       xlab=xlab, ylab=ylab, main="", ...)
##---------------
##--- draw symbols: ---
  z <- z/max(abs(z), na.rm = TRUE)
  usr <- par("usr")
  sxz <- fx * abs(z)
  syz <- fy * z
  segments(x - sxz, y - syz,  x + sxz, y + syz, lwd = slwd, col=scol)
##--- mark restricted observations: ---
  if(any(restr)) {
    points((x - sxz)[restr], (y - syz)[restr], pch= 8, mkh = 1/40)
    points((x + sxz)[restr], (y + syz)[restr], pch= 8, mkh = 1/40)
  }
  if (length(main)>0) mtext(main, 3, 1, cex=cex.title*par("cex"))
  stamp(sure=FALSE)
"plres2x done"
}
## ==========================================================================
plfitpairs <- function(object, ssize=0.02, main=NULL) #, pch=NULL
{
  ## Purpose:   pairs plot of fitted values for multinomial regression
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  5 Aug 2004, 10:54
  if (is.null(main)) main <- paste("fitted prob.",object$formula)
  lpr <- object$fitted.values
  lny <- ncol(lpr)
  ly <- object$y
  if(length(ly)==0) stop("!plfitpairs! no response values found")
  ly <- as.numeric(factor(object$y))
##-   if (is.factor(ly)) ly <-  as.numeric(factor())
  if (max(ly)!=lny)
    stop("!plfitpairs! ncol of fitted values != number of levels in y")
##  if (length(pch)<lny) pch <- 1:lny
  lmx <- max(lpr)
  l.panel <- function(x,y,indx,indy,ly,clr, ssize) {
    lix <- indx==ly
    liy <- indy==ly
    x[!(lix|liy)] <- NA
    segments(x-ssize*lix,y-ssize*liy,x+ssize*lix,y+ssize*liy,col=clr)
    abline(1,-1,lty=3)
  }
  plmatrix(lpr, panel=l.panel, pch=ly, range.=c(0,lmx), main=main, ssize=ssize)
  "plfitpairs done"
}
## ===================================================================
plTA.polr <- function(object, colbars=grey(0.7), colref=grey(0.7),
                         ploty=FALSE)
{
  ## Purpose:   plot "conditional median" residuals against fit for
  ##            cumulative logit model
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   object   result of polr
  ##   colbars  color to be used for plotting residuals
  ##   colref   color to be used for plotting the reference line
  ##   ploty    if TRUE, the latent response will be plotted instead of the
  ##            residuals
  ## Remark:    experimental function !!!
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  2 Oct 2007, 12:18
  if (length(object$call$weights))
      warning(":plTA.polr: Weigths are not visualized in the plot.")
  lyy <- residuals.polr(object)
  lfit <- lyy[,"fit"]
  if (ploty) lyy <- lyy[,c("median","lowq","uppq")]+lfit
  plot(range(lfit),range(lyy[,1:3]),type="n", # robrange(lyy[,1:3],fac=4),
       xlab="fit",ylab=if (ploty) "latent variable" else "residual")
  if (ploty) abline(0,1,col=colref) else
    abline(h=c(0,-1,1)*qlogis(0.975),col=colref)
  for (lk in 1:length(object$zeta))
    if (ploty) abline(h=object$zeta[lk],col=colref,lty=5) else
       abline(object$zeta[lk],-1,col=colref,lty=5)
  segments(lfit,lyy[,"lowq"],lfit,lyy[,"uppq"],col=colbars)
  points(lfit,lyy[,"median"],pch="-",cex=1.5)
  ls <- loess(lyy[,"median"]~lfit,span=0.7)
  lx <- seq(min(lfit),max(lfit),length=51)
  lsy <- predict(ls,newdata=data.frame(lfit=lx),family="symmetric")
  lines(lx,lsy,col="red")
}
## ===========================================================================
## additional useful functions
## ===========================================================================
dropdata <- function(data, rowid=NULL, incol="row.names", colid=NULL)
{
  ## Purpose:   drop observations from a data frame
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel
  li <- lj <- NULL
  lattr <- attributes(data)
  lattr <- lattr[is.na(match(names(lattr),
                             c("dim","dimnames","row.names","names")))]
  ln <- NROW(data)
  if (!is.null(rowid)) {
    lrn <- RNAMES(data)
    if (is.null(lrn)) lrn <- as.character(1:NROW(data))
    if (incol=="row.names")
      li <- match(as.character(rowid),lrn,nomatch=0)
    else {
      incol <- if (is.numeric(incol)) (1:ncol(data))[incol] else
      match(incol, colnames(data))
      if (is.na(incol)) stop("misspecified argument `incol`")
      li <- match(rowid,data[,incol],nomatch=0)
    }
    if (any(li==0)) warning(":dropdata: observations",
              paste(rowid[li==0],collapse=", "),"not found")
    li <- li[li>0]
    if (!is.null(li)) {
      data <- cbind(data)[-li,]
      names(li) <- lrn[li]
    }
  }
  ## drop variables
  if (!is.null(colid)) {
    lj <- match(as.character(colid),names(data),nomatch=0)
    if (any(lj==0)) warning(":dropdata: variables  ",
              paste(colid[lj==0],collapse=", "),"  not found")
    lj <- lj[lj>0]
    if (!is.null(lj)) data <- data[,-lj,drop=FALSE]
  }
  if (length(li)==0&length(lj)==0) {
      warning(":dropdata: no data to be dropped")
      return(data)
    }
  if (length(li)) {
    if (length(li)==NROW(data)) warning(":dropobs: no observations left")
    if (length(lattr$na.action))  {
      lin <- which(naresid(lattr$na.action, 1:ln%in%li))
      names(lin) <- lrn[li]
      li <- c(lattr$na.action, lin)
    }
    class(li) <- "exclude"
    lattr$na.action <- li
  }
  attributes(data) <- c(attributes(data),lattr)
  data
}
## ======================================================================
subset <- function(x, ...) {
  lattr <- attributes(x)
  lattr <- lattr[!names(lattr)%in%c("dim","dimnames","row.names","names")]
  lsubs <- base::subset(x, ...)
  attributes(lsubs) <- c(attributes(lsubs),lattr)
  lsubs
}
## ======================================================================
showd <- function(data, first=3, nrow.=4, ncol.=NULL)
{
## print some rows (and columns) of a matrix or data.frame
  ldoc <- getUserOption("doc")
  if (length(ldoc)>0 && ldoc && length(tit(data))>0) {
    cat("tit: ",tit(data),"\n")
  }
  lldim <- length(dim(data))
  if (lldim>2) stop("!showd not yet programmed for arrays")
  if (lldim>0) cat("dim: ",dim(data),"\n") else
    if (is.factor(data)) data <- as.character(data)
  ldata <- cbind(data)
  l.nr <- nrow(ldata)
  l.nc <- ncol(ldata)
  ## select columns
  l.ic <- if (length(ncol.)==0) 1:l.nc  else {
    if (length(ncol.)==1) {
      if (l.nc>ncol.)
        c(seq(1,by=l.nc%/%ncol.,length=ncol.-1),l.nc) else 1:l.nc
    } else  {
      lic <- ncol.[ncol.>0&ncol<=l.nc]
      if (length(lic)>0) lic else 1:l.nc
    }
  }
  ## select rows
  if (l.nr<=nrow.+first)  l.dc <- format(ldata[,l.ic])  else {
    l.ir <- c(1:first,round(seq(first,l.nr,length=nrow.+1))[-1])
    l.ir <- unique(c(last(l.ir,-1),l.nr))
    l.dc <- data.frame(u.merge(format(ldata[l.ir,l.ic]),"",after=first),
                       stringsAsFactors=FALSE)
    names(l.dc) <- colnames(ldata)[l.ic]
    lrn <- row.names(ldata)
    if (is.null(lrn)) lrn <- paste("r",1:l.nr,sep=".")
    row.names(l.dc) <- c(lrn[1:first],"...", lrn[l.ir[-(1:first)]])
  }
  ## was vector or array with only 1 column
  if (l.nc==1) {
    if (lldim>0) cat("     transposed column\n")
    row.names(l.dc) <-
      format(rbind(row.names(l.dc),l.dc[,1]),justify="right")[1,]
    l.dc <- t(l.dc)
  }
  print(l.dc,quote=FALSE)
  if (length(ldoc)&&ldoc&&length(doc(data)))
    cat("\ndoc:  ",paste(doc(data),collapse="\n  "),"\n")
  invisible(l.dc)
}
## -------------------------------------------------------------------------
mframe <-
function(mfrow=NULL, mfcol=NULL, mft=NULL, row=TRUE, oma=c(0,0,2,1),
                 mar=getUserOption("mar"), mgp=getUserOption("mgp"), ...)
{
## Purpose:    par(mfrow...)
## Author: Werner Stahel, 1994 / 2001
  if (is.null(mft)) {
    if (is.null(mfrow)) mfrow <- 1
    if (is.null(mfcol)) mfcol <- 1
  } else {
    t.din <- par("din")
    if (is.null(mfrow)) mfrow <- max(1,ceiling(sqrt(mft*t.din[2]/t.din[1])))
    mfcol <- ceiling(mft/mfrow)
  }
  mfrow <- max(1,mfrow)
  mfcol <- max(1,mfcol)
  t.oma <- if (mfrow*mfcol>1) oma else rep(0,4)
  if (length(mar)==0) mar <- c(3,3,1,1)+0.5 else
  mar <- rep(mar, length=4)
  if (length(mgp)!=3) mgp <- c(2,0.8,0) else
  invisible(if(row)
            par(mfrow=c(mfrow,mfcol), oma=oma, mar=mar, mgp=mgp, ...) else
            par(mfcol=c(mfrow,mfcol), oma=oma, mar=mar, mgp=mgp, ...) )
}
## ==========================================================================
robrange <-
  function(data, trim=0.2, fac=3, na.rm=TRUE)
{
  lna <- any(!is.finite(data))
  if (lna) {
    if(!na.rm) stop("!robrange! 'data' contains NAs")
    data <- data[is.finite(data)]
  }
  ln <- length(data)
  if (is.character(data)|length(data)==0) stop("!robrange! invalid data")
  trim <- c(trim, 0.2)[1]
  if (!is.finite(trim)) trim <- 0.2
  lmn <- mean(data,trim=trim)
  lds <- sort(abs(data-lmn))
  lnt <- ceiling((1-trim)*ln)
  if (lnt<3 | lnt==ln) {
    warning(":robrange: not enough valid data. returning ordinary range")
    lsd <- Inf } else {
    lsd <- fac*sum(lds[1:lnt]/(lnt-1))
    if (lsd==0) {
      warning(":robrange: robust range has width 0. returning ordinary range")
      lsd <- Inf }
  }
  c(max(lmn-lsd,min(data)), min(lmn+lsd,max(data)))
}
## ==========================================================================
stamp <- function(sure=TRUE, outer.margin = NULL,
                  project=getUserOption("project"), step=getUserOption("step"),
                  stamp=getUserOption("stamp"), ...)
{
## Purpose:   plot date and project information
## -------------------------------------------------------------------------
## Arguments:
##   sure     if F, the function only plots its thing if  getOption("stamp")>0
##   outer    if T, the date is written in the outer margin
##   project  project title
##   step     title of step of data analysis
##   ...      arguments to  mtext , e.g., line=3
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: 13 Aug 96, 09:00
  if (length(stamp)==0) {
    message("stamp() setting userOptions(stamp=1)")
    userOptions(stamp=1)
    stamp <- 1
  }
  if (length(outer.margin)==0) outer.margin <- par("oma")[4]>0
  t.txt <- date()
  t.txt <- paste(substring(t.txt,5,10),",",substring(t.txt,22,23),"/",
                 substring(t.txt,13,16),sep="")
    if (length(project)>0) t.txt <- paste(t.txt,project,sep=" | ")
    if (length(step)>0) t.txt <- paste(t.txt,step,sep=" | ")
  if( sure | stamp==2 | ( stamp==1 & (
##     last figure on page
     { t.tfg <- par("mfg") ; all(t.tfg[1:2]==t.tfg[3:4]) }
       || (is.logical(outer.margin)&&outer.margin) ))  )
       mtext(t.txt, 4, cex = 0.6, adj = 0, outer = outer.margin, ...)
  invisible(t.txt)
}
## =======================================================================
quinterpol <- function(x, probs = c(0.25,0.5,0.75), extend=TRUE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 15 Nov 2014, 16:04
  lx <- x[!is.na(x)]
  ltb <- table(lx)
  ln <- length(lx)
  lnn <- length(ltb)
  ln1 <- lnn+1
  lxx <- as.numeric(names(ltb))
  lxm <- (lxx[-1]+lxx[-lnn])/2
  lx0 <- if(extend) 2*lxx[1]-lxm[1] else lxx[1]
  lx1 <- if(extend) 2*lxx[lnn]-lxm[lnn-1] else lxx[lnn]
  lxe <- c(rbind(c(lx0,lxm),lxx),lx1)
  lp <- c(0,cumsum(ltb)/ln)
  lpp <- (lp[-1]+lp[-ln1])/2
  lpe <- c(rbind(lp,c(lpp,1)))  ## last element (1) is ineffective
  ld <- outer(probs,lpe,"-")
  li <- apply(ld>0,1,sum)
  lii <- 1:length(probs)
  ldd <- cbind(ld[cbind(lii,li)],ld[cbind(lii,li+1)])
  lh <- ldd[,1]/(ldd[,1]-ldd[,2])
  lxe[li]*(1-lh) + lxe[li+1]*lh
}
## ===========================================================================
getUserOption <- function (x, default = NULL)
{
    if (is.null(default))
        return(userOptions(x))
    if (x %in% names(userOptions()))
        userOptions(x)
    else default
}
userOptions <- function (x=NULL, default=NULL, list=NULL, ...)
    {
##-     lpos <- find("UserOptions")
##-     luopt <- get("UserOptions", pos=lpos)
    luopt <- if (exists("UserOptions", where=1)) get("UserOptions", pos=1) else
         UserOptions            
    if ((!is.null(x)&&is.character(x))) ## asking for options
      return(if(length(x)==1) luopt[[x]] else luopt[x])
    if (!is.null(default)) {
        default <- as.character(default)
      if (default=="TRUE"|default=="all") return(userOptions(list=UserDefault))
      if (default=="unset")
        userOptions(list=UserDefault[names(UserDefault)%nin%names(luopt)])
      if (!is.character(default))
        stop("!userOptions! Unsuitable argument  default .")
      return(userOptions(list=UserDefault[
                           default[default%in%names(UserDefault)]]))
    }
    lop <- c(list,list(...))
    ## show all options
    if (length(lop)==0) return(luopt)
    ## set options
    lold <- luopt[names(lop)]
    for (li in names(lop))
      luopt[li] <- list(lop[[li]])
    assign("UserOptions", luopt,pos=1)
        ## assignInMyNamespace does not work
    invisible(lold)
  }
## -----------------------------------------------------
UserDefault <- UserOptions <- 
  list(stamp=1, project="", step="", doc=TRUE, show.dummy.coef=TRUE,
       colors.ra = c("black","gray4","blue4","cyan","darkgreen","green",
         "burlywood4","burlywood3","burlywood4"),
       mar=c(3,3,3,1), mgp=c(2,0.8,0), digits=4,
       regr.contrasts=c(unordered="contr.sum", ordered="contr.poly"),
       regr.testcoefcol=c("coef", "stcoef", "df", "R2.x", "signif", "p.value"),
       debug=0
       )
##- if (!exists("UserOptions")) UserOptions <- UserDefault  else
##-     userOptions(default="unset")
## ===========================================================================
last <-
function(data,n = NULL, ncol=NULL, drop=is.matrix(data))
{
  ldim <- dim(data)
  if (is.null(ldim)) {
    if (is.null(n)) n <- 1
    ldt <- length(data)
    if (is.null(n)) n <- ldt
    return(data[sign(n)*((ldt-abs(n)+1):ldt)])
  }
  if (length(ldim)!=2)
    stop ("!last! not programmed for arrays of dimension >2")
  if (is.null(n)&is.null(ncol)) n <- 1
  if (is.null(n)) n <- ldim[1]
  if (is.null(ncol)) ncol <- ldim[2]
  data[sign(n)*((ldim[1]-abs(n)+1):ldim[1]),
       sign(ncol)*((ldim[2]-abs(ncol)+1):ldim[2]), drop=drop]
}
## ==============================================================
nainf.exclude <- function (object, ...)
  ## na.omit, modified to omit also Inf and NaN values
{
  if (is.atomic(object)) {
    i <- is.finite(object)
    if (length(dim(i))) ## matrix 
      return( object[apply(i,1,all),,drop=FALSE] )
    else return( object[i] )
  }
  ## list
  n <- length(object)
    omit <- FALSE
    vars <- seq_len(n)
    for (j in vars) {
        x <- object[[j]]
        if (!is.atomic(x))
            next
##-         x <- is.na(x)
        x <- if (is.numeric(x)) !is.finite(x) else is.na(x)
        d <- dim(x)
        if (is.null(d) || length(d) != 2)
            omit <- omit | x
        else for (ii in 1:d[2]) omit <- omit | x[, ii]
    }
    xx <- object[!omit, , drop = FALSE]
    if (any(omit > 0L)) {
        temp <- seq(omit)[omit]
        names(temp) <- attr(object, "row.names")[omit]
        attr(temp, "class") <- "exclude"
        attr(xx, "na.action") <- temp
    }
    xx
  }
## ===================================================
nna <- function(x,inf=TRUE) if (inf) x[is.finite(x)] else x[!is.na(x)]
## ===================================================
sumna <- function(object,inf=TRUE)
{
  ## Purpose:   count NAs along columns
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   object   data.frame, matrix or vector
  ##   inf      treat Inf as NA
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 10 Oct 2007, 08:18
  ff <- if(inf) {
    function(x)
    if(is.numeric(x)) sum(!is.finite(x)) else sum(is.na(x)) }
      else function(x) sum(is.na(x))
  if (is.matrix(object)) apply(object,2,ff)  else {
    if (is.list(object)) sapply(object,ff)
    else if(is.atomic(object)) ff(object)
  }
}
## ==========================================================================
logst <- function(data, calib=data, threshold=NULL, mult=1)
{
  ## Purpose:   logs of data, zeros and small values treated well
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  3 Nov 2001, 08:22
  data <- cbind(data)
  calib <- cbind(calib)
  lncol <- ncol(calib)
  ljthr <- length(threshold)>0
  if (ljthr) {
    if (!length(threshold)%in%c(1, lncol))
      stop("!logst! length of argument `threshold` is inadequate")
    lthr <- rep(threshold, length=lncol)
    ljdt <- !is.na(lthr)
  } else {
    ljdt <- rep(TRUE, lncol)
    lthr <- rep(NA, lncol)
    for (lj in 1:lncol) {
      lcal <- calib[,lj]
      ldp <- lcal[lcal>0&!is.na(lcal)]
      if(length(ldp)==0) ljdt[lj] <- FALSE else {
        lq <- quantile(ldp,probs=c(0.25,0.75),na.rm=TRUE)
        if(lq[1]==lq[2]) lq[1] <- lq[2]/2
        lthr[lj] <- lc <- lq[1]^(1+mult)/lq[2]^mult
      }
    }
  }
  ## transform data
  for (lj in 1:lncol) {
    ldt <- data[,lj]
    lc <- lthr[lj]
    li <- which(ldt<lc)
    if (length(li))
      ldt[li] <- lc * 10^((ldt[li]-lc)/(lc*log(10)))
    data[,lj] <- log10(ldt)
  }
  if (length(colnames(data)))
    lnmpd <- names(ljdt) <- names(lthr) <- colnames(data)  else
    lnmpd <- as.character(1:lncol)
  if (ncol(data)==1) data <- data[,1]
  attr(data,"threshold") <- lthr
  if (any(!ljdt)) {
    warning(":logst: no positive data for variables ",lnmpd[!ljdt],
            ". These are not transformed")
    attr(data,"untransformed") <- c(ljdt)
  }
  data
}
## ===========================================================================
asinperc <- function(x) asin(sqrt(x/100))/asin(1)
## ===========================================================================
plcoord <-
function(x, range=NULL, limfac=3.0, limext=0.1)
{
  ## Purpose:    values for plot with limited "inner" plot range
  lrg <- if (length(range)==0||any(is.na(range)))
    robrange(x, fac=limfac)  else  range(range)
  if (diff(lrg)==0) lrg <- c(-range,range)
  rr <- pmax(pmin(x,lrg[2]),lrg[1])
  lxd <- x-rr
  lnmod <- sum(lxd!=0,na.rm=TRUE)
  if (lnmod>0) rr <- rr+(lxd)/(1+abs(lxd)/(diff(lrg)*limext))
  attr(rr,"range") <- lrg
  attr(rr,"nmod") <- lnmod
  class(rr) <- class(x)
  rr
}
## ===========================================================================
legendr <- function(x=0.05,y=0.95,legend, ...) {
  lusr <- par("usr")
  lx <- lusr[1] + x*diff(lusr[1:2])
  ly <- lusr[3] + y*diff(lusr[3:4])
  legend(lx,ly,legend, ...)
}

## ===========================================================================
doc <- function(x) attr(x,"doc")
## ---
"doc<-" <- function(x, value)
{
  ##-- Create doc attribute or  PREpend  new doc to existing one.
  value <- as.character(value)
  attr(x, "doc") <- if (length(value)==0) NULL else
  if(value[1]=="^") value[-1] else c(value, attr(x, "doc"))
  x
}
## ---
tit <- function(x) attr(x,"tit")
## ---
"tit<-" <- function(x, value) ## ! argument must be `value`. demanded by attr
{
  attr(x, "tit") <- value
  x
}
## ---
is.formula <- function(object)
  length(class(object))>0 && class(object)=="formula"
## =================================================================
## auxiliary functions
## ============================================================
nafalse <- function(x) if (is.null(x)) FALSE else ifelse(is.na(x), FALSE, x)

u.true <- function(x) length(x)>0 && x[1]
u.debug <- function() u.true(getUserOption("debug"))
u.merge <- function(dd1, dd2 = NA, which=NULL, after=NULL,
                    length=NULL, names=NULL)
{
## Purpose:   merge two vectors or expand a vector by NA s
## -------------------------------------------------------------------------
## Arguments:
##   dd1      first vector or matrix or data.frame (?),
##   dd2      second vector, ...
##   which    is T for indices for which first vector is used
##   after    elements of  dd2  will be inserted after "after" in dd1
##   length   length of the result (will be expanded if necessary)
##   names    names of the result (if length is adequate)
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: 11 Mar 93, 13:50, and later
  llen <- length
  n1 <- length(dd1)
  nc1 <- ncol(dd1)
  nc2 <- ncol(dd2)
  if (length(nc1)>0) {
    n1 <- nrow(dd1)
    if (!( length(dd2)==1 || is.null(nc2) || nc2==nc1 ))
      stop("unsuitable second argument")
    }
## --- generate  which  vector for all cases
  if (length(which)==0) {
## - after  specified
      if (length(after)==0) stop("specify either  which  or  after")
      if (is.logical(after))  after <- which(after)
      wh <- rep(TRUE,n1+length(after))
      wh[after+1:length(after)] <- FALSE }
  else {
## - which  specified
    if(is.logical(which)) wh <- which
    else {
      if (length(llen)==0)  llen <- n1+length(which)
        wh <- rep(TRUE, llen)
        wh[which] <- FALSE }
  }
## --- merge
  nn <- length(wh)
  n2 <- nn-n1
  if (!(is.null(names)|length(names)==nn))
    warning("argument  names  not used (unsuitable length)")
  if (length(nc1)>0) {
    if (!(length(dd2)==1 || NROW(dd2)==n2))
      stop("unsuitable number of rows")
    rr <- matrix(NA,nn,nc1)
    rr[wh,] <- as.matrix(dd1)
    rr[!wh,] <- if (is.data.frame(dd2)) as.matrix(dd2) else dd2
##-     if (length(names)>0) row.names(rr) <- names else {
##-       if (length(lrn1 <- row.names(dd1))>0)
  }
  else {
    rr <- rep(NA,nn)
    rr[wh] <- dd1
    rr[!wh] <- dd2
    if (length(names)>0) names(rr) <- names
  }
  rr
}
## ============================================================
i.main <- function(main, line=1-outer.margin, cex=NULL, adj=NULL,
                   outer.margin=FALSE, col="black",
                   doc=getOption("doc"))
{
  ## Purpose:   title
## ----------------------------------------------
  ladj <- if (length(adj)==0) 0.5 else adj
  if (length(cex)==0) {
    cex <- max(0.5, min(1.2,
      par("mfg")[4]*par("pin")[1]/(par("cin")[1]*nchar(main))))
    if (cex==0.5) ladj <- if (length(adj)==0) 0 else adj
  }
  if (is.na(outer.margin)) outer.margin <- par("oma")[3]>0
  if (outer.margin && 1!=prod(par("mfg")[1:2])) return()
  if (length(main)!=0)
    mtext(main, 3, line, cex = cex*par("cex"), adj=ladj, outer = outer.margin,
          col=col)
  if ((!is.null(doc))&&doc&&length(tit(main)))
    mtext(tit(main), 3, line+1, outer = outer.margin, col=col)
}
## ==========================================================================
##- is.R <- function ()
##- exists("version") && !is.null(vl <- version$language) && vl == "R"
RNAMES <- function(x) if (!is.null(dim(x))) row.names(x) else names(x)
"%nin%" <- function(x,y) !x%in%y
getmeth <- function(fn,mt)  getS3method(as.character(substitute(fn)),
                                        as.character(substitute(mt)))
warn <- function()   table(names(warnings()))
BR <- browser
DB <- function(on=TRUE) options(error=if(on) recover else NULL)
options(show.dummy.coef=TRUE)
IR <- function(condition) { if (condition) {
  cat("INTERRUPT: ",as.character(substitute(condition)))
  traceback()
  browser()
}
}
## ===========================================================================
if (length(getUserOption("colors"))==0)
  userOptions(colors = c("black","firebrick3","deepskyblue3","springgreen3",
         "darkgoldenrod3","olivedrab3","purple3","orange3","palegreen3"))
if (length(getUserOption("colors.ra"))==0)
  userOptions(colors.ra =
          c("black","gray","blue","cyan","darkgreen","green",
            "burlywood4","burlywood3","burlywood4"))
c.weekdays <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
            "Saturday", "Sunday")
c.months <- c("January", "February", "March", "April", "May", "June",
        "July", "August", "September", "October", "November",
        "December")
c.mon <- substring(c.months,1,3)

drevgumbel <- function (x, location = 0, scale = 1)
{ # from VGAM
    if (!nafalse(scale>0))
        stop("\"scale\" must be positive")
    E <- exp((x - location)/scale)
    E * exp(-E)/scale
}
prevgumbel <- function (q, location = 0, scale = 1)
{
    if (!nafalse(scale>0))
        stop("\"scale\" must be positive")
    -expm1(-exp((q - location)/scale)) # expm1(u) = exp(u)-1, accurately also for |u| << 1
}
qrevgumbel <- function (p, location = 0, scale = 1)
{
    if (!nafalse(scale>0))
        stop("\"scale\" must be positive")
    location + scale * log(-log(p))
}

qrevgumbelexp <- function (p) exp(qrevgumbel(p))

rrevgumbel <- function (n, location = 0, scale = 1)
{
    if (!nafalse(n>=1))
        stop("bad input for argument \"n\"")
    if (!nafalse(scale>0))
        stop("\"scale\" must be positive")
    location + scale * log(-log(runif(n)))
}
## ===========================================================================
## repaired versions of drop functions
## ===========================================================================
##- if(getRversion() <= "2.7.1") {

add1.default <- function(object, scope, scale = 0, test=c("none", "Chisq"),
			 k = 2, trace = FALSE, ...)
{
    if(missing(scope) || is.null(scope)) stop("no terms in scope")
    if(!is.character(scope))
	scope <- add.scope(object, update.formula(object, scope))
    if(!length(scope))
	stop("no terms in scope for adding to object")
##     newform <- update.formula(object,
##                               paste(". ~ . +", paste(scope, collapse="+")))
##     data <- model.frame(update(object, newform)) # remove NAs
##     object <- update(object, data = data)
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2,
                  dimnames = list(c("<none>", scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
    n0 <- length(object$residuals)
    env <- environment(formula(object))
    for(i in seq(ns)) {
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying +", tt, "\n", sep="")
	    utils::flush.console()
	}
	nfit <- update(object, as.formula(paste("~ . +", tt)),
                       evaluate = FALSE)
	nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)
	ans[i+1, ] <- extractAIC(nfit, scale, k = k, ...)
        if(length(nfit$residuals) != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[,1] - ans[1,1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[,2])
    test <- match.arg(test)
    if(test == "Chisq") {
	dev <- ans[,2] - k*ans[, 1]
	dev <- dev[1] - dev; dev[1] <- NA
	nas <- !is.na(dev)
	P <- dev
	P[nas] <- pchisq(dev[nas], dfs[nas], lower.tail=FALSE)
	aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    head <- c("Single term additions", "\nModel:",
	      deparse(as.vector(formula(object))),
	      if(scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
## ==================================================================
drop1.default <- function(object, scope, scale = 0, test=c("none", "Chisq"),
			  k = 2, trace = FALSE, ...)
{
    tl <- attr(object$terms, "term.labels")
    if(missing(scope)) scope <- drop.scope(object)
    else {
	if(!is.character(scope))
	    scope <- attr(terms(update.formula(object, scope)), "term.labels")
	if(!all(match(scope, tl, 0L) > 0L))
	    stop("scope is not a subset of term labels")
    }
##    data <- model.frame(object) # remove NAs
##    object <- update(object, data = data)
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2,
                  dimnames =  list(c("<none>", scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
    n0 <- length(object$residuals)
    env <- environment(formula(object))
    for(i in seq(ns)) {
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying -", tt, "\n", sep="")
	    utils::flush.console()
        }
        nfit <- update(object, as.formula(paste("~ . -", tt)),
                       evaluate = FALSE)
	nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)
	ans[i+1, ] <- extractAIC(nfit, scale, k = k, ...)
        if(length(nfit$residuals) != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[1,1] - ans[,1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[,2])
    test <- match.arg(test)
    if(test == "Chisq") {
        dev <- ans[, 2] - k*ans[, 1]
        dev <- dev - dev[1] ; dev[1] <- NA
        nas <- !is.na(dev)
        P <- dev
        P[nas] <- pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
        aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    head <- c("Single term deletions", "\nModel:",
	      deparse(as.vector(formula(object))),
	      if(scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
# ==================================================================
i.polrfit <-
function (formula, data, weights, start, ..., subset, na.action,
          contrasts = NULL, Hess = FALSE, model = FALSE, x = TRUE,
          method = c("logistic", "probit", "cloglog", "cauchit"))
  ## copy of polr from MASS.
  ## argument  x  added: keep model.matrix
  ## 1 line added by WSt to keep eta in the result
  ## ::: change argument name  method:::
{
    logit <- function(p) log(p/(1 - p))  ## !!!???!!! use qlogis() --- why?
    fmin <- function(beta) {
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))),
            100)
        eta <- offset
        if (pc > 0)
            eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        if (all(pr > 0))
            -sum(wt * log(pr))
        else Inf
    }
    gmin <- function(beta) {
        jacobian <- function(theta) {
            k <- length(theta)
            etheta <- exp(theta)
            mat <- matrix(0, k, k)
            mat[, 1] <- rep(1, k)
            for (i in 2:k) mat[i:k, i] <- etheta[i]
            mat
        }
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))),
            100)
        eta <- offset
        if (pc > 0)
            eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        p1 <- dfun(gamm[y + 1] - eta)
        p2 <- dfun(gamm[y] - eta)
        g1 <- if (pc > 0)
            t(x) %*% (wt * (p1 - p2)/pr)
        else numeric(0)
        xx <- .polrY1 * p1 - .polrY2 * p2
        g2 <- -t(xx) %*% (wt/pr)
        g2 <- t(g2) %*% jacobian(theta)
        if (all(pr > 0))
            c(g1, g2)
        else rep(NA, pc + q)
    }
    m <- match.call(expand.dots = FALSE)
    method <- match.arg(method)
    pgumbel <- function(q) exp(pweibull(log(q))) # ???
    dgumbel <- function(q) stop("BUG: dgumbel not programmed")
    pfun <- switch(method, logistic = plogis, probit = pnorm,
        cloglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logistic = dlogis, probit = dnorm,
        cloglog = dgumbel, cauchit = dcauchy)
    if (is.matrix(eval.parent(m$data)))
        m$data <- as.data.frame(data)
    m$start <- m$Hess <- m$method <- m$model <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    x.ret <- x  ## ! Wst need to copy x because x is used for the matrix
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    n <- nrow(x)
    pc <- ncol(x)
    cons <- attr(x, "contrasts")
    if (xint > 0) {
        x <- x[, -xint, drop = FALSE]
        pc <- pc - 1
    }
    else warning("an intercept is needed and assumed")
    wt <- model.weights(m)
    if (!length(wt))
        wt <- rep(1, n)
    offset <- model.offset(m)
    if (length(offset) <= 1)
        offset <- rep(0, n)
    y <- model.response(m)
    if (!is.factor(y))
        stop("response must be a factor")
    lev <- levels(y)
    if (length(lev) <= 2)
        stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- length(lev) - 1
    Y <- matrix(0, n, q)
    .polrY1 <- col(Y) == y
    .polrY2 <- col(Y) == y - 1
    if (missing(start)) {
      q1 <- length(lev)%/%2
      y1 <- (y > q1)
      X <- cbind(Intercept = rep(1, n), x)
      fit <- switch(method,
                    logistic = glm.fit(X, y1, wt, family = binomial(),
                      offset = offset),
                    probit = glm.fit(X, y1, wt, family = binomial("probit"),
                      offset = offset),
                    cloglog = glm.fit(X, y1, wt, family = binomial("probit"),
                      offset = offset),
                    cauchit = glm.fit(X, y1, wt, family = binomial("cauchit"),
                      offset = offset)
                    )
      if (!fit$converged) { ## new attempt
        fit <- lm.fit(X, as.numeric(y), wt)
        fit$coefficients <- fit$coef/sqrt(mean(fit$resid^2))
        warning("attempt to find suitable starting values may have failed")
        }
      coefs <- fit$coefficients
      if (any(is.na(coefs))) {
        warning("design appears to be rank-deficient, so dropping some coefs")
        keep <- names(coefs)[!is.na(coefs)]
        coefs <- coefs[keep]
        x <- x[, keep[-1], drop = FALSE]
        pc <- ncol(x)
      }
      spacing <- logit((1:q)/(q + 1))
      if (method != "logistic")
        spacing <- spacing/1.7
      gammas <- -coefs[1] + spacing - spacing[q1]
      thetas <- c(gammas[1], log(diff(gammas)))
      s0 <- c(coefs[-1], thetas)
    }
    else if (length(start) != pc + q)
        stop("\"start\" is not of the correct length")
      else {
        s0 <- if (pc > 0)
            c(start[seq_len(pc + 1)], diff(start[-seq_len(pc)]))
        else c(start[1], diff(start))
    }
    res <- optim(s0, fmin, gmin, method = "BFGS", hessian = Hess,
        ...)
    beta <- res$par[seq_len(pc)]
    theta <- res$par[pc + 1:q]
    zeta <- cumsum(c(theta[1], exp(theta[-1])))
    deviance <- 2 * res$value
    niter <- c(f.evals = res$counts[1], g.evals = res$counts[2])
    names(zeta) <- paste(lev[-length(lev)], lev[-1], sep = "|")
    if (pc > 0) {
        names(beta) <- colnames(x)
        eta <- drop(x %*% beta)
    }
    else {
        eta <- rep(0, n)
    }
    cumpr <- matrix(pfun(matrix(zeta, n, q, byrow = TRUE) - eta),
        , q)
    fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
    dimnames(fitted) <- list(row.names(m), lev)
    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance,
                linear.predictor = eta,  ## added by WSt
        fitted.values = fitted, lev = lev, terms = Terms, df.residual = sum(wt) -
            pc - q, edf = pc + q, n = sum(wt), nobs = sum(wt),
        call = match.call(), method = method, convergence = res$convergence,
        niter = niter)
    if (Hess) {
        dn <- c(names(beta), names(zeta))
        H <- res$hessian
        dimnames(H) <- list(dn, dn)
        fit$Hessian <- H
    }
    if (model)  fit$model <- m
    if (x.ret)  fit$x <- x  ## !WSt  return model.matrix
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    class(fit) <- "polr"
    fit
}
## -------------------------------------------------------
i.multinomfit <- 
function (formula, data, weights, subset, na.action, contrasts = NULL, 
    Hess = FALSE, summ = 0, censored = FALSE, model = FALSE, x = TRUE,
    ...)
    ## copy of multinom from MASS. Argument x added by WSt
{
##  require(nnet)
    class.ind <- function(cl) {
        n <- length(cl)
        x <- matrix(0, n, length(levels(cl)))
        x[(1L:n) + n * (as.integer(cl) - 1L)] <- 1
        dimnames(x) <- list(names(cl), levels(cl))
        x
    }
    summ2 <- function(X, Y) {
        X <- as.matrix(X)
        Y <- as.matrix(Y)
        n <- nrow(X)
        p <- ncol(X)
        q <- ncol(Y)
        Z <- t(cbind(X, Y))
        storage.mode(Z) <- "double"
        z <- .C(nnet:::VR_summ2, as.integer(n), as.integer(p), as.integer(q), 
            Z = Z, na = integer(1L))
        Za <- t(z$Z[, 1L:z$na, drop = FALSE])
        list(X = Za[, 1L:p, drop = FALSE], Y = Za[, p + 1L:q])
    }
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$summ <- m$Hess <- m$contrasts <- m$censored <- m$model <- m$... <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    XX <- X <- model.matrix(Terms, m, contrasts) ## !WSt  store X
    cons <- attr(X, "contrasts")
    Xr <- qr(X)$rank
    Y <- model.response(m)
    if (!is.matrix(Y)) 
        Y <- as.factor(Y)
    w <- model.weights(m)
    if (length(w) == 0L) 
        if (is.matrix(Y)) 
            w <- rep(1, dim(Y)[1L])
        else w <- rep(1, length(Y))
    lev <- levels(Y)
    if (is.factor(Y)) {
        counts <- table(Y)
        if (any(counts == 0L)) {
            empty <- lev[counts == 0L]
            warning(sprintf(ngettext(length(empty), "group %s is empty", 
                "groups %s are empty"), paste(sQuote(empty), 
                collapse = " ")), domain = NA)
            Y <- factor(Y, levels = lev[counts > 0L])
            lev <- lev[counts > 0L]
        }
        if (length(lev) < 2L) 
            stop("need two or more classes to fit a multinom model")
        if (length(lev) == 2L) 
            Y <- as.integer(Y) - 1
        else Y <- class.ind(Y)
    }
    if (summ == 1) {
        Z <- cbind(X, Y)
        z1 <- cumprod(apply(Z, 2L, max) + 1)
        Z1 <- apply(Z, 1L, function(x) sum(z1 * x))
        oZ <- order(Z1)
        Z2 <- !duplicated(Z1[oZ])
        oX <- (seq_along(Z1)[oZ])[Z2]
        X <- X[oX, , drop = FALSE]
        Y <- if (is.matrix(Y)) 
            Y[oX, , drop = FALSE]
        else Y[oX]
        w <- diff(c(0, cumsum(w))[c(Z2, TRUE)])
        print(dim(X))
    }
    if (summ == 2) {
        Z <- summ2(cbind(X, Y), w)
        X <- Z$X[, 1L:ncol(X)]
        Y <- Z$X[, ncol(X) + 1L:ncol(Y), drop = FALSE]
        w <- Z$Y
        print(dim(X))
    }
    if (summ == 3) {
        Z <- summ2(X, Y * w)
        X <- Z$X
        Y <- Z$Y[, 1L:ncol(Y), drop = FALSE]
        w <- rep(1, nrow(X))
        print(dim(X))
    }
    offset <- model.offset(m)
    r <- ncol(X)
    if (is.matrix(Y)) {
        p <- ncol(Y)
        sY <- Y %*% rep(1, p)
        if (any(sY == 0)) 
            stop("some case has no observations")
        if (!censored) {
            Y <- Y/matrix(sY, nrow(Y), p)
            w <- w * sY
        }
        if (length(offset) > 1L) {
            if (ncol(offset) != p) 
                stop("ncol(offset) is wrong")
            mask <- c(rep(FALSE, r + 1L + p), rep(c(FALSE, rep(TRUE, 
                r), rep(FALSE, p)), p - 1L))
            X <- cbind(X, offset)
            Wts <- as.vector(rbind(matrix(0, r + 1L, p), diag(p)))
            fit <- nnet.default(X, Y, w, Wts = Wts, mask = mask, 
                size = 0, skip = TRUE, softmax = TRUE, censored = censored, 
                rang = 0, ...)
        }
        else {
            mask <- c(rep(FALSE, r + 1L), rep(c(FALSE, rep(TRUE, 
                r)), p - 1L))
            fit <- nnet.default(X, Y, w, mask = mask, size = 0, 
                skip = TRUE, softmax = TRUE, censored = censored, 
                rang = 0, ...)
        }
    }
    else {
        if (length(offset) <= 1L) {
            mask <- c(FALSE, rep(TRUE, r))
            fit <- nnet.default(X, Y, w, mask = mask, size = 0, 
                skip = TRUE, entropy = TRUE, rang = 0, ...)
        }
        else {
            mask <- c(FALSE, rep(TRUE, r), FALSE)
            Wts <- c(rep(0, r + 1L), 1)
            X <- cbind(X, offset)
            fit <- nnet.default(X, Y, w, Wts = Wts, mask = mask, 
                size = 0, skip = TRUE, entropy = TRUE, rang = 0, 
                ...)
        }
    }
    fit$formula <- attr(Terms, "formula")
    fit$terms <- Terms
    fit$call <- call
    fit$weights <- w
    fit$lev <- lev
    fit$deviance <- 2 * fit$value
    fit$rank <- Xr
    edf <- ifelse(length(lev) == 2L, 1, length(lev) - 1) * Xr
    if (is.matrix(Y)) {
        edf <- (ncol(Y) - 1) * Xr
        if (length(dn <- colnames(Y)) > 0) 
            fit$lab <- dn
        else fit$lab <- 1L:ncol(Y)
    }
    fit$coefnames <- colnames(X)
    fit$vcoefnames <- fit$coefnames[1L:r]
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    fit$edf <- edf
    fit$AIC <- fit$deviance + 2 * edf
    if (model) 
        fit$model <- m
    class(fit) <- c("multinom", "nnet")
    if (Hess) 
        fit$Hessian <- nnet:::multinomHess(fit, X)
    if (x) fit$x <- XX ## !Wst return design matrix
    fit
}
## ================================================================
colorpale <- function(col=NA, pale=0.3, ...)
  {
    if (is.na(col)) {
      col <- palette()[2]
      warning(":colorpale: Argument 'col' is NA. I assume  ", col)
    }
    crgb <- t(col2rgb(col)/255)
    rgb(1-pale*(1-crgb), ...)
  }
