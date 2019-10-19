##  regr.R  Functions that are useful for regression, W. Stahel,
## ==========================================================================
regr <-
  function (formula, data=NULL, family=NULL, 
            robust = FALSE, method=NULL, 
            nonlinear = FALSE, start=NULL,
            subset=NULL, weights=NULL, offset=NULL, ...)
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
  ## --- a. get arguments for  regr.control
  lcall <- match.call()
  lac <- as.list(lcall)[-1]
  laRegr <- c("formula", "data", "family", "robust", "method",
              "nonlinear", "start", "subset", "weights", "offset")
  laControl <- c("contrasts","factorNA", "na.action", "calcdisp", "suffmean",
                 "model", "x", "termtable", "vif", "testlevel", "leveragelimit",
                 "dist", "control")
  latransfer <- c("na.action", "calcdisp", "suffmean",
                 "model", "x", "termtable", "vif")
  lextra <- setdiff(names(lac), c(laRegr, laControl))
  if (length(lextra))
    warning(":regr: argument(s)  ",paste(lextra, collapse=", "),
            "  not used")
  lac <- lac[names(lac)%nin%lextra]
  lcl <- c(list(quote(regr.control)), 
            lac[match(laControl,names(lac), nomatch=0)])
  mode(lcl) <- "call"
  largs <-eval(lcl)
  ## ------------------------------------------------------------------  
  ## b. ===   preparation:
  ## convert character formula to formula
  lform <- lformula <- as.formula(formula)
  lcall$formula <- lformula
  ## nonlinear: drop constants
  nonlinear <- i.def(nonlinear, FALSE, TRUE, FALSE)
  if (as.character(nonlinear)!="FALSE") 
    lform <- setdiff(all.vars(lform), names(start))
  ## nonlinear <- i.nonlincheck(nonlinear, lformula, ldata)
  ## ------------------------------------------------------------------  
  ## d. === data
  lextrav <- as.list(lcall)[names(lcall)%in%c("weights", "offset", "subset")]
  lcgetv <- c(list(quote(plgraphics::getvariables)),
              list(formula = lform, data = data, transformed = FALSE),
              lextrav)
  mode(lcgetv) <- "call"
  lallvars <- eval(lcgetv)
  ## ---------------------
  lallvars <- nainf.exclude(lallvars)
  lnaaction <- attr(lallvars, "na.action")
  ## --------------------------------------
  ## f. === compose call of fitting function
  lcl <- lcall
  ## !!! extras: check, replace names
  if (length(lextrav)) {
    lextras <- 
      c(weights=".weights.", offset=".offset.", subset=".subset.")[names(lextrav)]
    lcl[names(lextras)] <- lallvars[lextras] ## was intersect(names(lallvars),c(".weights.",".offset.",".subset."))]
  }
  ## missing response
  if (length(lform)>2) {
    lvy <- all.vars(lform[1:2])
    ly0 <- lallvars[,lvy]
    linna <- apply(cbind(ly0), 1, is.finite)
  } else linna <- TRUE
  if (length(linna)>1 && sum(linna)<2)
    stop("!regr! Less than 2 non-missing response values") 
  ## --- convert character to factor, drop unused levels, generate .NA. level
  lfacna <- i.def(largs$factorNA, TRUE)
  lfnalabel <- if(is.character(lfacna)) lfacna else ".NA."
  for (lvn in 1:ncol(lallvars)) {
    lv <- lallvars[[lvn]]
    if (is.logical(lv))
      lallvars[[lvn]] <- as.numeric(lv) ## logical -> numeric!!!
    if (is.character(lv)|is.factor(lv)) {
      if (lfacna) lv <- factorNA(lv, lfnalabel)[[1]]
      lallvars[[lvn]] <- factor(lv, levels=levels(factor(lv[linna])))
    }
  }
  ## g. --- check for variables with a single value
##-   lv1 <- which( apply(lmodelframe, 2, function(x) all(x==x[1]) ) )
##-   if (length(lv1)) {  ## adjust formula
##-     lfac <- attr(terms(lmodelframe),"fac")
##-     lt1 <- names(which(apply(lfac[lv1,,drop=F],2,any)))
##-     warning("!regr! formula contains single valued variables: ",
##-             paste(row.names(lfac)[lv1], collapse=". "),
##-             "\n     I drop the following terms from the formula:\n    ",
##-             paste(lt1, collapse=", "))
##-     lfupd <- as.formula( paste( ".~ .- ",paste(lt1, collapse=" - ") ) )
##-     lcl$formula <- update(lcl$formula, lfupd)
##-   }
  ## -------------------------------------------
  ## h. === response type
  if (length(lformula)==2) { # nonlinear called with formula of type ~...
    ly <- rep(0,NROW(lallvars))
    lytype <- "numeric"
  } else {
##-     attr(formula, "response")
    lyf <- model.frame(lformula[1:2], lallvars, na.action=na.pass)
    ## I tried to generate model.frame for x and y together. This failed
    ## because model.frame  needs adequate method (when y is matrix)
    ltrm <- attr(lyf, "terms")
    lytype <- substring(attr(ltrm, "dataClasses"),1,5)
##  lysimple <- lytype!="nmatr" ## not a matrix
    lyy <- lyf[[1]]
    lysimple <- length(dim(lyy))==0
##    ly <- na.omit(lyy)
    if (lysimple&&length(unique(dropNA(lyy)))==2 &&
        all(as.numeric(lyy)%in%0:1)) ## FALSE for numeric !={0,1}
      lytype <- "binary" 
    if (inherits(lyy,"Surv"))  {
        lytype <- "survival"
    }
## strange variables
##-   l1v <- sapply(ldta, function(x) all(x==c(x[!is.na(x)],0)[1],na.rm=TRUE) )
##-                                 ## covers case of several or all NAs
##-   if (any(l1v)) {
##-     warning(paste(":regr: variable(s)", paste(lvars[l1v],collapse=", "),
##-                   "has (have) no distinct values")) #  -> dropped.
##-   }
  }
  ## ----------------------------------------------
  ## k. === family and fitting function
  lfam <- i.def(as.character(substitute(family)), largs$dist)[1]
  lcl$dist <- NULL
  if (lytype=="survival")
    lfam <- c( lfam, attr(lyy,"distribution"))[1]
  if (length(lfam)==0 || any(is.na(lfam)))
    lfam <- switch(substring(lytype,1,5),
                   numer="normal", nmatr="normal", binar="binomial",
                   binco="binomial", order="cumlogit",
                   facto="multinomial", survi="ph", "unknown")
  if (substring(lfam,1,7)=="multinom") lfam <- "multinomial"
  if (lfam=="multinom") lfam <- "binomial"
  ##
  lfitfun <-
      switch( lfam,
             gaussian="lm", normal="lm", binomial="glm", poisson="glm",
             Gamma="glm",
             cumlogit="polr", multinomial="multinomial",
             weibull="survreg", lognormal="survreg", loggaussian="survreg",
             loglogis="survreg", loglogistic="survreg", extreme="survreg",
             ph="survreg", prop.hazard="survreg",
             "unknown")
  if (lfitfun=="unknown") stop("!regr! Fitting function not identified")
  ## additional checks
  if (lytype=="survival") {
    if (!inherits(lyy,"Surv"))
      stop("!regr! bug: convert response to Surv object")
    ## !!! hier machen! lallv[,1] ersetzen durch Surv davon
    lfitfun <- "survreg"
  }
  else  if (lfitfun=="glm")
    lcl$control <- list(calcdisp=largs$calcdisp, suffmean=largs$suffmean,
                        lcl$control)
  ## 
  lfitname <- paste("i",lfitfun,sep=".")
  if (!exists(lfitname)||!is.function(get(lfitname)))
    stop (paste("!regr! Fitting function",lfitname, "not found"))
  ## -----------------------------------------------------
  ## m. === prepare call
  lcl$fname <- lfam
  ##  lcl$na.action <- substitute(largs$na.action)
  ## lcl <- c(list(quote(regr)), as.list(lcl[-1]), largs[latransfer])
  lcl[latransfer] <- largs[latransfer]
  lcl[[1]] <- ## hack --> eval(.) works also when call is source()d ...
      switch(lfitname,
	     "i.lm" = quote(regr::i.lm),
	     "i.glm" = quote(regr::i.glm),
	     "i.multinomial" = quote(regr::i.multinomial),
	     "i.polr" = quote(regr::i.polr), 
##	     "i.smooth" = quote(regr::i.smooth), ## ??
	     "i.survreg" = quote(regr::i.survreg),
	     ## default:
	     as.name(lfitname))
##!!!  if (lfitname=="i.glm") lcl$family <- lfam
  if (lfitname=="i.survreg") {
    lcl$yy <- lyy
    lcl$model <- TRUE  ## model needed, see below
  }
  ## --- contrasts
  lcontr <- largs$contrasts
  if(is.atomic(lcontr)&&length(lcontr)) {
    if(!is.character(lcontr))
      warning("!regr! invalid contrasts argument")
    else {
      loldopt <- options(contrasts=c(lcontr,getOption("contrasts")[2])[1:2])
      on.exit(options(loldopt))
      lcl$contrasts <- NULL
    }
    lcw <- lcontr==c("contr.wsum","contr.wpoly")
    if (ncol(lallvars)>1) {
      if (any(lcw))
        lyna <- if (length(dim(lyy)))
                  c(0,NA)[1+apply(is.na(as.matrix(lyy)),1,any)] else lyy 
      for (lj in 2:ncol(lallvars)) { ## no contrasts for y {
        if(lcw[1]&&class(lallvars[,lj])[1]=="factor")
          attr(lallvars[,lj],"contrasts") <-
            contr.wsum(lallvars[,lj], y=lyna)
        if(lcw[2]&&class(lallvars[,lj])[1]=="ordered")
          attr(lallvars[,lj],"contrasts") <-
            contr.wpoly(lallvars[,lj], scores=NULL, y=lyna)
      }
    }
##    if (lcontr[1]=="contr.wsum") lallvars <- contr.wsum(lallvars, y=lyy)
  }
  lcl$data <- lallvars ## must be evaluated!
  lcall$na.action <- lcl$na.action <- largs$na.action
  mode(lcl) <- "call"
  ## === --------------------------------------------
  ##-   lreg <- eval(lcl, envir=environment(formula))
  lreg <- eval(lcl)
  ## === --------------------------------------------
  if (is.null(lreg$distrname)) lreg$distrname <- lfam
  if (length(lreg$AIC)==0) {
    laic <- try(extractAIC(lreg), silent=TRUE)
    if (class(laic)!="try-error") lreg$AIC <- laic
  }
  lreg$response <- lyy
##-     if (length(lnaaction)) {
##-       if (is.matrix(lyy)) lyy[-lnaaction,] else lyy[-lnaaction]} else lyy
  lreg$allvars <- lallvars
  ## if (length(lnaaction)) lallvars[-lnaaction,] else lallvars
  ## needed more than $model since $model contains transformed variables
  ## --- recover some arguments to effective function call
  lfc <- lreg$call
  la <- intersect(c("data", "weights", "offset", "subset","na.action"), names(lfc))
  ## these arguments should be restored because otherwise,
  ## add1 does not work if they have changed.
  lfc[la] <- lcall[la]
  lreg$funcall <- lfc
  lcall$formula <- formula(lreg) # hope this never damages anything
  lreg$call <- lcall
  tit(lreg) <- if (length(largs$tit)==0) attr(data,"tit") else largs$tit
  doc(lreg) <- attr(data,"doc")
  if (largs$model&&length(lreg$model)==0) {
      if (nonlinear) warning(":regr: no $model available for nonlinear regr.")
      else lreg$model <- lm(lformula, data, method="model.frame")
  }
  lterms <- if (nonlinear) NULL else terms(lreg)
  if ((!nonlinear) && is.null(attr(lterms, "predvars")))  ## needed for survreg
    attr(lreg$terms,"predvars") <- attr(attr(lreg$model,"terms"),"predvars")
  ## -----------------------------------------------------------------
  ## r. === leverages, standardized res
  ## get residuals if missing (as for  polr  objects
  if (!inherits(lreg, "multinom")) {
    lj <- match("resid", names(lreg), nomatch=0)
    if (lj>0) names(lreg)[lj] <- "residuals"
    if ("residuals" %nin% names(lreg)) {
      lreg$residuals <- residuals(lreg)
##-       if (length(lnaaction) && class(lnaaction)=="exclude")
##-         lres <- if (is.matrix(lres)) lres[-lnaaction,] else lres[-lnaaction]
##-       lreg$residuals <- lres
    }
    lsigma <- c(lreg$sigma, lreg$scale)[1]
    if (length(lsigma)==0) lsigma <- sqrt(c(lreg$dispersion,1)[1])
    if (!inherits(lreg, "nls")) {
##-       lrg <- lreg
##-       if(length(lnaaction)) lrg$na.action <- structure(lnaaction, class="omit")
      lstr <- stdresiduals(lreg, sigma=lsigma, leveragelimit = largs$leveragelimit)
      lat <- attributes(lstr)
      attributes(lstr)[c("leverage","stdresratio","stddev")] <- NULL
      lreg$stdresiduals <- lstr
      lreg$stdresratio <- lat$stdresratio
      lreg$leverage <- lat$leverage
    }
  }
  ## --- misc
  if (nonlinear) lreg$r.squared <- 1-lreg$sigma^2/var(lyy,na.rm=TRUE)
  if (!largs$x) lreg$x <- NULL
  if (length(lnaaction)&length(largs$na.action))
    class(lnaaction) <-
      sub("na.","", sub("nainf.","", as.character(largs$na.action)))
  lreg$na.action <- lnaaction
  class(lreg) <- if (class(lreg)[1]=="orig")  ##  nls shall not be regr
    class(lreg)[-1] else c("regr",class(lreg))
  ## ------------------------------------------------------------------
  ## result of regr
  lreg
}
## -----------------------------------------------------------------------
regr.control <-
  function (contrasts=i.getoption("regr.contrasts"), factorNA = NULL, 
           na.action=as.name("nainf.exclude"), calcdisp=NULL, suffmean=3,
           dist=NULL,
           model = FALSE, x = TRUE, termtable=NULL, vif=NULL,
           testlevel = NULL, leveragelimit=NULL, tit=NULL,
           control = NULL
           )
{
  factorNA <- i.getopt(factorNA)
  testlevel <- i.getopt(testlevel)
  termtable <- i.getopt(termtable)
  vif <- i.getopt(vif)
  list(contrasts=contrasts, factorNA=factorNA,
       na.action=as.name("nainf.exclude"), calcdisp=calcdisp,
       suffmean=suffmean,
       dist=dist, model=TRUE, x=x, termtable=termtable, vif=vif,
       testlevel=testlevel, leveragelimit=leveragelimit, tit=tit,
       control=control
       )
  ## flicken !!! model=T needed in i.lm_ for getting ly
}
## ===================================================================
## =========================================================================
i.lm <-
  function (formula, data, family, fname="gaussian", nonlinear=FALSE,
            robust=FALSE, method=NULL, control=NULL, 
            vif=TRUE, termtable=TRUE, testlevel=0.05, call = NULL, ...)
{
  ## Purpose:  internal: fit lm
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
  ## b. --- method
  ## NEW version feb 20
  lcall <- match.call()
  lmeth <- c(lcall$method,"")[1]
  lfn <- if (nonlinear) {
    lcall$contrasts <- NULL ## lcall[-match("contrasts",names(lcall),nomatch=0)]
    "nls" } else {
    if (lmeth=="lmrob") robust <- TRUE
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
        lcall$setting <- "KS2014"
      }
      lcall$x <- TRUE
  }
    if (lfn=="rlm") {
##      require(MASS)    ##  !?!
      lcall$method <- c(method,"MM")
      lcall$x.ret <- TRUE
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
	     lmrob = quote(robustbase::lmrob),
	     rlm = quote(MASS::rlm),
             rq = quote(quantreg::rq),
             lm = quote(stats::lm),
             nls = quote(stats::nls),
	     ## default:
	     as.name(fn))
  }
  ##  lcall[[1]] <- mkFn(lfn)
  if(lfn!="lmrob") lcall$control <- NULL
  lcall <- lcall[setdiff(names(lcall),
                         c("fname","family","vif","nonlinear","robust",
                           "calcdisp","suffmean","termtable"))] #,"control"
  lcl <- c(list(mkFn(lfn)),as.list(lcall)[-1])
  mode(lcl) <- "call"
  ## --------------------------
  lreg <- eval(lcl, envir=environment())
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
##- ##-       if (length(lhat)!=NROW(lreg$stdres))
##- ##-         if (length(lreg[["w"]])==NROW(lreg$stdres))
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
  if (is.null(lsig)) lsig <- sd(lreg$resid) # !!! used for rq
  lreg$sigma <- lsig
##-   ## standardized residuals
##-   if (is.finite(lsig)&&lsig>0) {
##-     lreg$stdres <- lreg$residuals/lsig
##-     if (length(lreg$weights)) lreg$stdres <- lreg$stdres*sqrt(lreg$weights)
  if (class(lreg)=="lmrob") lreg1$cov.unscaled <- lreg$cov/lsig^2 ## !!!
  ## from summary
  lcomp <- c("r.squared","fstatistic","colregelation","aliased",
             "df","cov.unscaled")
  lreg[lcomp] <- lreg1[lcomp]
  if (lfn=="lm") lreg$AIC <- extractAIC(lreg)[2]
  ## degrees of freedom
  if (is.null(lreg$df)) # needed for rq
    lreg$df <- c(length(coef(lreg))-attr(terms(lreg),"intercept"),
                 length(lreg$residuals)-length(coef(lreg)))
  lreg$df.residual <- ldfr <- df.residual(lreg)
  ## coef table
  lcftab <- lreg1$coefficients
  if (NCOL(lcftab)>1) {
    lcf <- lcftab[,1]
    attr(lcf, "se") <- lcftab[,2]
    lreg$coefficients <- lcf
    lreg$coeftable <- ciSignif(lcftab, df=lreg$df.residual, testlevel=testlevel)
##-     if (nonlinear) {
##-       lreg$termtable <- 
##-     }
#   lreg$r.squared <- 1-(lsig/lsdy)^2
  }
  ##
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
      lreg[names(ltt)] <- ltt
    }  else  class(lreg) <- c("orig",class(lreg))
  }
  ## result of i.lm
  lreg
}
## -----------------------------------------------------------------------
i.mlmsum <- function (object, termtable=TRUE)
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
  lres <- object$residuals
##-   if (all(lsig>0)) {
##-     object$stdres <- sweep(lres,2,lsig,"/")
##-     if (length(object$weights))
##-       object$stdres <- object$stdres*sqrt(object$weights)
##-   }
  object$resmd <- mahalanobis(lres,0,var(lres))
  ldfr <- object$df.residual
  object$r.squared <- lr2 <- lts["r.squared",]
  object$adj.r.squared <- 1-(1-lr2)*(nrow(object$residuals)-1)/ldfr
  lcomp <- c("aliased","df","cov.unscaled")
  object[lcomp] <- lreg1[[1]][lcomp]
  object$drop1 <- if (lmodel) drop1.mlm(object)
##  class(lreg) <- c("mregr","mlm","lm")
  object
} # i.mlmsum
## -----------------------------------------------------------------------
i.glm <-
  function (formula, data, family, fname,
            control=NULL, vif=TRUE, termtable=TRUE, call=NULL, ...)
{
  ## Purpose:  internal: fit glm
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
  lfamily <- get(fname)
##-   environment(formula) <- environment()
  lcall <- match.call() 
  lcall$x <- TRUE
  lcall <- lcall[setdiff(names(lcall),
                         c("fname","vif","nonlinear","robust",
                           "calcdisp","suffmean","termtable","control"))]
  lcall <- c(list(quote(stats::glm)),as.list(lcall[-1]))
  mode(lcall) <- "call"
  ## ---------------
  lreg <- eval(lcall, envir=environment())
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
  if (ldisp>1) {
      lstr <- residuals(lreg, type="pearson")/sqrt(ldisp)
      lnaa <- lreg$na.action
      if (class(lnaa)=="exclude") lstr <- lstr[-lnaa]
      lreg$stdres <- lstr
  }
  ## bug? leverage not taken into account
  lcomp <- c("deviance","AIC","df.residual","null.deviance", # "family",
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
##-     lcmpn <- c("termtable","termeffects","leverage")
    ##-     lreg[lcmpn[lcmpn%in%names(ltt)]] <- ltt
    lreg[names(ltt)] <- ltt
  }
  ## result of i.glm
  lreg
}
## -----------------------------------------------------------------------
i.multinomial <-
  function (formula, data, family, fname,
            model=TRUE, vif=TRUE, termtable=TRUE, call=NULL, ...)
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
  lcall[[1]] <- quote(regr::i.multinomfit)
  lcall$fname <- lcall$family <- lcall$control <- lcall$vif <- NULL
  lcall$trace <- FALSE
  lreg <- eval(lcall, envir=environment())
  ## ---------------
  if (length(lreg$na.action)) {
    lnaact <- attr(lreg$na.action,"class")
    attr(lreg$na.action,"class") <- "omit"
  } else  lnaact <- NULL ##  summary does not work with  exclude
  lreg$call$formula <- formula
  lreg1 <- summary(lreg)
  lreg$dispersion <- lreg$sigma <- 1
  lres <- lreg1$residuals
  lreg$residuals <- lres
  lcf <- lreg1$coefficients
  lreg$coefficients <- lcf
  lreg$AIC <- lreg1$AIC
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
  } else { ##xxx
    ldr1 <- ldr1[-1,]}
  ## signif :
  ldr1 <- cbind( ldr1[,1:3], sqrt(ldr1[,3]/qchisq(0.95,ldr1[,1])), ldr1[,4] )
  names(ldr1) <- c("df", "AIC", "Chisq", "signif", "p.value") 
  lreg$termtable <- lreg$drop1 <- ldr1
  if (length(lnaact)) attr(lreg$na.action,"class") <- lnaact
## result of i.multinomial
  lreg
}
## -----------------------------------------------------------------------
i.polr <-
  function (formula, data, family, fname, weights = NULL, 
            model=TRUE, vif=TRUE, termtable=TRUE, call=NULL, ...)
{
  ## Purpose:  internal: fit ordered y
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
##  require(MASS)   ## !?!
  lcall <- match.call()
  ladrop <- c("fname","family","vif","nonlinear","robust",
              "calcdisp","suffmean","termtable")
  lcl <- as.list(lcall[setdiff(names(lcall),ladrop)])
  lcl[1] <- list(quote(regr::i.polrfit))
##  lcl <- c(list(quote(regr::i.polrfit), lcl)
  mode(lcl) <- "call"
  lcl$Hess <- TRUE
  lenv <- environment()
  lcl$envir <- lenv
## ---
  lreg <- eval(lcl, envir=lenv)
##  lreg$call$formula <- formula
  lreg$w <- data$.weights.
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
  lreg$stdres <- NULL
  ## cov of estimates!
  lreg$covariance <- lcov <- vcov(lreg)
  lse <- sqrt(diag(lcov))
  lreg$correlation <- lcov/outer(lse, lse)
  ## --- deviances
  lreg$fitfun <- "polr"
  if (termtable) {
    ltt <- i.termtable(lreg, lreg1$coef, data, lcov, ltesttype="Chisq",
                       lsdy=1, vif=vif, leverage=TRUE)
    lreg[names(ltt)] <- ltt
##-     lcmpn <- c("termtable","termeffects","leverage")
##-     lreg[lcmpn[lcmpn%in%names(ltt)]] <- ltt
  }
  lreg$dispersion <- 1
  ##  lreg$residuals <- residuals(lreg)
  llp <- fitted.polr(lreg, type="link")
  if (length(lnaaction <- lreg$na.action)) llp <- llp[-lnaaction]
  lreg$linear.predictors <- llp
  ## result of i.polr
  lreg
}
## -----------------------------------------------------------------------
i.survreg <-
  function (formula, data, family, yy, fname="ph", method, control,
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
  lyy <- lcall$yy
  if (is.null(lcall$dist)) lcall$dist <- lcall$family
  lcall <- lcall[names(lcall)%nin%c("yy","fname","vif","family",
                                    "calcdisp","suffmean","termtable")]
##  lcall$yy <- lcall$fname <- lcall$family <- lcall$vif <- NULL
  lcall$y <- lcall$x <- TRUE
  ## ---
  lreg <- eval(lcall, envir=environment())
  lreg$n.obs <- nobs(lreg)
  ## ---
  class(lreg$y) <- "Surv"
  attr(lreg$y, "type") <- attr(lyy, "type")
##  lreg$call$formula <- formula
  lreg1 <- if (u.debug()) summary(lreg) else
           try(summary(lreg), silent=TRUE)
  if (class(lreg1)[1]=="try-error") {
    warning(paste(":regr/i.survreg: summary did not work. ",
                  "I return the survreg object"))
##    lreg$call$data <- call$data
    class(lreg) <- c("orig",class(lreg))
    return(lreg)
  }   ## ---
  lreg$resid.orig <- lreg$residuals
  lreg$stdres <- NULL
  lcf <- lreg1$coefficients
  ## --- deviances
  ## lreg$scale
  if (lfitfun=="survreg") {
    attr(lreg$scale,"fixed") <- length(lcall$scale)>0
    ldf <- sum(lreg$df) - lreg$idf
    ldfr <- length(lreg$linear.predictors)-sum(lreg$df)
  }
  if (lfitfun=="coxph") {
    lreg1$table <- lreg1$coefficients
    ldf <- length(lreg$coefficients)
    ldfr <- length(lreg$residual)-ldf-1
  }
  lreg$df.residual <- ldfr
  lreg$AIC <- extractAIC(lreg)[2]
  lreg$deviance <- -2*lreg$loglik
  lchi <- 2*diff(lreg1$loglik)
  ltbd <- cbind(deviance=c(lchi,-2*lreg1$loglik[2]),
                df=c(ldf, ldf+ldfr),
                p.value=c(pchisq(lchi,ldf,lower.tail=FALSE),NA))
  dimnames(ltbd)[[1]] <- c("Model","Null")
  lreg$devtable <- ltbd
  lreg$covariance <- lcov <- lreg$var
  lsd <- sqrt(diag(lcov))
  lreg$correlation <- lcov/outer(lsd,lsd)
  lreg$fitfun <- lfitfun
  lnaaction <- lreg$na.action
  lfit <- lreg$linear.predictors
  lres <- if (inherits(lreg, "survreg")) {
            residuals.regrsurvreg(lreg)  ## includes NAs
            } else residuals.regrcoxph(lreg)
  if (length(lnaaction) && class(lnaaction)=="exclude")
    lres <- lres[-lnaaction]
  ly <- lreg$y
##-   lreg$n.censored <-
##-     if (attr(ly,"type")%in%c("right","left"))
##-       table(ly[,2])[2] else  sum(ly[,2]!=1)    #interval
  ltype <- attr(ly,"type")
  ##-  lreg$n.censored <- sum(lres[,"prob"]>0, na.rm=TRUE)
  ltb <- table(ly[,2])
  llimit <- attr(ly,"limit")
  lreg$n.censored <- NA
  if (ltype=="left") {
    lreg$n.censored <- structure(ltb[2], names="left")
    if (length(llimit))
      lreg$n.fitout <- structure(sum(lfit<llimit, na.rm=TRUE), names="left")
  } else {
    if (ltype=="right") {
      lreg$n.censored <- structure(ltb[1], names="right")
    if (length(llimit))
    lreg$n.fitout <- structure(sum(lfit>llimit, na.rm=TRUE), names="right")
    } 
  }
  if (termtable) {
    ltt <- i.termtable(lreg, lreg1$table, data, lcov, ltesttype="Chisq",
                       lsdy=1, vif=vif)
    lreg[names(ltt)] <- ltt
    ## log(scale): signif<-NA. no! log(scale)==0 means
    ##    exp.distr for weibull/gumbel
##-     lcmpn <- c("termtable","termeffects","leverage")
##-     lreg[lcmpn[lcmpn%in%names(ltt)]] <- ltt
  }
##-   ## lreg$df <- c(model=ldf, residual=ldfr, original=lreg$df) ## !!! lreg has df = ldf+object$idf !
##-     ## do not modify before calling i.termtable
  lreg$distrname <- if (lfitfun=="coxph") "prop.hazard" else lreg$dist
  lreg$residuals <- lres
  ## result of i.survreg
  lreg
}

## -----------------------------------------------------------------------
i.polrfit <-
function (formula, data, weights, start, ..., subset, na.action,
          contrasts = NULL, Hess = FALSE, model = FALSE, x = TRUE,
          method = c("logistic", "probit", "cloglog", "cauchit"),
          envir = parent.frame())
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
    m$start <- m$Hess <- m$method <- m$model <- m$... <- m$envir <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, envir=environment())
    Terms <- attr(m, "terms")
    x.ret <- x  ## ! Wst need to copy x because x is used for the matrix
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    n <- nrow(x)
    pc <- ncol(x)
    asgn <- attr(x, "assign") ## ws
    cons <- attr(x, "contrasts")
    if (xint > 0) {
        x <- x[, -xint, drop = FALSE]
        pc <- pc - 1
        asgn <- asgn[-xint] ## ws
        attr(Terms, "intercept") <- 0 ## ws
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
      if (anyNA(coefs)) {
        warning("design appears to be rank-deficient, so dropping some coefs")
        keep <- names(coefs)[!is.na(coefs)]
        coefs <- coefs[keep]
        x <- x[, keep[-1], drop = FALSE]
        pc <- ncol(x)
        asgn <- asgn[keep[-1]] ## ws
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
    ##
    res <- optim(s0, fmin, gmin, method = "BFGS", hessian = Hess,
                 ...)
    ##
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
    ##
    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance,
                linear.predictor = eta,  ## added by WSt
                fitted.values = fitted, lev = lev, terms = Terms,
                df.residual = sum(wt) - pc - q, edf = pc + q, n = sum(wt),
                nobs = sum(wt),
                call = match.call(), method = method,
                convergence = res$convergence,
                niter = niter)
    if (Hess) {
        dn <- c(names(beta), names(zeta))
        H <- res$hessian
        dimnames(H) <- list(dn, dn)
        fit$Hessian <- H
    }
    if (model)  fit$model <- m
    if (x.ret)
      fit$x <- structure(x, assign=asgn) ## !WSt  return model.matrix
##        structure(cbind("(Intercept)"=1, x), assign=c("(Intercept)"=0, asgn) ) 
    ## , contrasts=cons  probably not needed
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
i.termtable <-
  function (lreg, lcoeftab, ldata, lcov, ltesttype="F",
            lsdy, vif=TRUE, leverage=vif)
{
  ## Purpose:  generate term table for various models
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 15:37
##-   if (length(lnaaction <- lreg$naaction) && class(lnaaction)=="exclude")
##-     lreg$na.action <- structure(lnaaction, class="omit")
##-  lreg$na.action <- NULL
  lterms <- terms(lreg)
  if(length(attr(lterms,"term.labels"))==0)
    return(list(termtable = data.frame(
      coef=c(lreg$coef,NA)[1], se=NA, ciLow=NA, ciUp=NA, 
      df=1, testst=NA, signif=NA, p.value=NA, p.symbol="", stcoef=NA, R2.x=NA,
      stringsAsFactors=FALSE)
                ))
## degrees of freedom
  ldfr <- df.residual(lreg)
  if (ldfr<1) {
    warning(":regr/i.termtable: no degrees of freedom left.")
    return(list(termtable = data.frame(
      coef=c(lreg$coef,NA)[1], se=NA, ciLow=NA, ciUp=NA, 
      df=1, testst=NA, signif=NA, p.value=NA, p.symbol="", stcoef=NA, R2.x=NA,
      stringsAsFactors=FALSE)
                ))
  }    
  pvCutpoints <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  pvSymbols <- c("***", "**", "*", ".", " ", "")
  pvLegend <- paste(rbind(pvCutpoints,pvSymbols), collapse="  ")
  ## drop1
  ldr1 <-
      if (class(lreg)[1]%in%c("lm","lmrob")) {
          if (u.debug()) 
              drop1Wald(lreg, test=ltesttype, scope=lterms) else
              try(drop1Wald(lreg, test=ltesttype, scope=lterms),
              silent=TRUE) } else {
          if (u.debug()) 
              drop1(lreg, test=ltesttype, scope=lterms) else
              try(drop1(lreg, test=ltesttype, scope=lterms),
                  silent=TRUE)
                           }
  if (class(ldr1)[1]=="try-error") {
    warning(paste(":regr: drop1 did not work. I return the table produced by ",
                  lreg$fitfun))
##-     lsum <- summary(lreg)
##-     lcft <- lsum$coef
##-     if (length(lcft)==0) lcft <- lsum$parameters ## nls
##-     return(list(test=lcft)) # !!! noch reparieren
    return(list(termtable=lcoeftab))
  }
  ldr1 <- ldr1[-1,]
  ldr1$RSS <- NULL # same ncol for lm and glm
  if (inherits(lreg,"rlm"))  ldr1[,4] <- ldr1[,2]/ldr1[,1] ## !!!
  if (inherits(lreg,"mlm")||inherits(lreg,"manova"))
    return(list(termtable=ldr1))  ## !!! needs much more
  ltstq <- if (ltesttype=="F") qf(0.95,c(1,ldr1[,1]),ldfr) else {
    if (ltesttype=="Chisq") qchisq(0.95,c(1,ldr1[,1])) else NA }
  ltstq1 <- sqrt(ltstq[1]) ## 1 degree of freedom
  ltstq <- ltstq[-1]
## coefficients
  lcoef <- lreg$coefficients
## model.matrix
  lmmt <- lreg[["x"]]
  if (length(lmmt)==0)
      lmmt <- model.matrix(lreg)
  lasg <- attr(lmmt,"assign")[!is.na(lcoef)]
##  if (class(lreg)[1]%in%c("polr")) lasg <- lasg[-1] ## ,"coxph"
  ## terms without factor involvement
  lfactors <- attr(lterms,"factors")
  lvcont <- !attr(lterms,"dataClasses")[row.names(lfactors)] %in%
    c("numeric","logical") ## [...] excludes .weights. and possibly others
  ## terms only containing continuous variables
  lcont <- which( lvcont %*% lfactors ==0 ) 
  ## licasg <- which(lasg%in%lcont)
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
  lpv <- ldr1[,lpvcol]
  ltb <- data.frame(coef=NA, se=NA, ciLow=NA, ciUp=NA, 
                    df=ldr1[,1], testst=ldr1[,lpvcol-1], signif=NA,
                    p.value=lpv, p.symbol="", stcoef=NA, R2.x=lr2,
                    stringsAsFactors=FALSE)
  row.names(ltb) <- row.names(ldr1)
  ## intercept
  ljint <- "(Intercept)"==names(lcoef)[1]
  if (ljint) {
##    ltstint <- # if(class(lreg)[1]%in%c("lm","nls","rlm"))
      lcoeftab[1,3]^2 # else lcoeftab[1,3]
    ltb <- rbind(
      "(Intercept)"=
        data.frame(coef=NA, se=NA, ciLow=NA, ciUp=NA, 
                   df=1, testst=NA, signif=NA,
                   p.value=NA, p.symbol="", stcoef=NA, R2.x=NA,
                   stringsAsFactors=FALSE),
      ltb)
    ltstq <- c(ltstq1, ltstq)
    lcont <- c(0, lcont)
  }
  lcont1 <- lcont+ljint  # row number in dr1
## p.symbol and signif
  ltb$signif <- sqrt(pmax(0,ltb$testst)/ltstq)
## coefficients and statistics for terms with 1 df
  if (length(lcont)) { ## lcont refers to assign
    ltlb <- dimnames(ltb)[[1]]
    lclb <- ltlb[lcont1] ## lcont1 is the row in the coef table of lreg1
    ljc <- match(lcont,lasg) # index of coefs for cont variables
    lcf <- lcoef[ljc]
    ## fill in
    ltb$coef[lcont1] <- lcf
    ltb$se[lcont1] <- lse <- lcoeftab[ljc,2]
    lci <- lcf+outer(ltstq1*lse, c(-1,1))
    ## confint(lreg,row.names(ltb)[lcont1]) does not always work...
    ltb[lcont1,c("ciLow","ciUp")] <- lci
    ltb[lcont1,"signif"] <- sign(lcf)*ltb[lcont1,"signif"]
    ## standardized coefficients
    lstcf <- lcf[lcont>0] *  # exclude intercept term
      sqrt(apply(lmmt[,names(lcf[lcont>0]),drop=FALSE],2,var)) / lsdy
    ltb$stcoef[lcont1[lcont>0]] <- lstcf
  }
  if (row.names(lcoeftab)[nrow(lcoeftab)]=="Log(scale)") { # survreg
    ltsc <- lcoeftab[nrow(lcoeftab),]
    lcont1 <- c(lcont1, nrow(lcoeftab))
    if (!u.true(lreg$dist=="weibull")) ltsc[2:4] <- NA
    ltb <- rbind(ltb,"log(scale)"=
                   c(ltsc[1],ltsc[2],ltsc[1]+c(-1,1)*qnorm(0.975)*ltsc[2],
                     1, ltsc[3], ltsc[3]/qnorm(0.975), ltsc[4], NA, NA, NA))
  }
## p-symbol
  lipv <- as.numeric(cut(ltb$p.value, pvCutpoints))
  ltb[,"p.symbol"] <- pvSymbols[lipv]
  attr(ltb, "legend") <- pvLegend
  class(ltb) <- c("termtable", "data.frame")
  ## --- termeffects (dummy coef)
  lallcf <- termeffects(lreg) 
  if (inherits(lreg,"polr")) lreg$coefficients <- c("(Intercept)" = NA, lcoef)
  rr <- list(termtable = ltb, termeffects = lallcf)
  if (leverage) rr <- c(rr, leverage=list(hat(lmmt)))
  rr
}
## --------------------------------------------------------------------------
ciSignif <- function (estimate, se=NULL, df=Inf, testlevel=0.05) {
  if (inherits(estimate, regrModelClasses))
    estimate <- summary(estimate)$coefficients
  if (is.null(se))
    if (NCOL(estimate)>1) {
      se <- estimate[,2]
      estimate <- estimate[,1]
    }
  if (is.null(se)) se <- attr(estimate, "se")
  if (is.null(se))
    stop("!ciSignif! no standard errors found")
  ltq <- qt(1-testlevel/2, df)
  lci <- estimate+outer(ltq*se, c(ciLow=-1,ciUp=1))
  ltst <- estimate/se
  lsgf <- ltst/ltq
  lpv <- 2*pt(-abs(ltst), df)
  lipv <- as.numeric(cut(lpv, c(0, 0.001, 0.01, 0.05, 0.1, 1)))
  lsst <- c("***", "**", "*", ".", " ")[lipv]
  data.frame(estimate=estimate, se=se, lci, testst=ltst,
             signif=lsgf, p.value=lpv, p.symbol=lsst)
}
  
## ==========================================================================
contr.wsumpoly <- 
  function (n, scores = NULL, y = NULL, w = NULL,
           contrasts = TRUE, sparse = FALSE, poly = NA) 
{  ## provide weighted sum contrasts
  if (is.data.frame(n)) {
    for (lj in 1:ncol(n))
      if (is.factor(n[,lj])) 
        attr(n[,lj],"contrasts") <-
          contr.wsumpoly(n[,lj], scores=scores, y=y,
                         contrasts=contrasts, sparse=sparse)
    return(n)
  }
  ## not a data.frame, but...
  if (is.character(n)) n <- factor(n)
  if (is.factor(n)) { ## ... a factor
    if (length(y)) {
      if (length(y)!=length(n)) {
        warning(":contrasts.wsum: unequal lengths of arguments. ",
                "I ignore argument 'y'")
        ## y only used to eliminate NAs in target variable
      } else  n <- n[apply(is.finite(cbind(y)), 1, all)]
    }
    w <- c(table(factor(n))) ## exclude unused levels
    if (is.na(poly)) poly <- is.ordered(n)
  }
  nn <- length(w)
  if (is.na(poly)) poly <- FALSE
  if (nn<1) {
    if (!(is.atomic(n)&&is.numeric(n)&&length(n)==1))
      stop ("!contr.wsumpoly! Provide either 'n' or 'w'")
    nn <- n
##    w <- rep(1,nn)
  }
  contr <-
    if (poly) {
      scores <- if(length(scores)) scores else 1:nn
      contr.poly(nn, scores = scores, contrasts=contrasts, sparse=sparse)
    } else
      contr.sum(nn, contrasts=contrasts, sparse=sparse)
  if (is.null(w) || anyNA(w) || any(w<=0)) {
    warning(":contr.wsum: weights 'w' not suitable.",
            " Returning unweighted sum contrast")
    return(contr)
  }
  if (poly) {
    contr <- make.poly( nn, scores=scores, w=w)
    if (contrasts) {
      dn <- colnames(contr)
      dn[2:min(4, nn)] <- c(".L", ".Q", ".C")[1:min(3, nn - 1)]
      colnames(contr) <- dn
      contr <- contr[, -1, drop = FALSE]
    }
    else {
      contr[, 1] <- 1
      contr
    }
  } else  contr[nn,] <- - w[-nn]/w[nn]
##-   if (sparse) 
##-     contr <- .asSparse(contr)
  structure(contr, w=w)
}
## --------------------------------------------------------------------
contr.wsum <-
  function (n, scores = NULL, y=NULL, w = NULL, contrasts = TRUE,
            sparse = FALSE)
{
  if (is.ordered(n)) n <- factor(n)
  contr.wsumpoly (n, y=y, w=w, contrasts=contrasts, sparse=sparse, poly=FALSE)
}
## --------------------------------------------------------------------
contr.wpoly <-
  function (n, scores = NULL, y = NULL, w = NULL, contrasts = TRUE,
            sparse = FALSE)
{
  if (is.factor(n)) n <- ordered(n)
  contr.wsumpoly (n, scores = scores, y = y, w = w,
                  contrasts = contrasts, sparse = sparse, poly=TRUE)
}
## -----------------------------------------------------------------
make.poly <- function (n, scores, w) {
        y <- scores - sum(scores*w)/sum(w)
        X <- sqrt(w)*outer(y, seq_len(n) - 1, "^")
        QR <- qr(X)
        z <- QR$qr
        z <- z * (row(z) == col(z))
        Z <- qr.qy(QR, z)  ## raw <- 
##-         Z <- sweep(raw, 2L, apply(raw, 2L, function(x) sqrt(sum(x^2))), 
##-             "/", check.margin = FALSE)
      ##      do not standardize. WSt
      Z <- Z / sqrt(w)
        colnames(Z) <- paste0("^", 1L:n - 1L)
        Z
}
## ==========================================================================
summary.regr <- function (object, ...)  object
plot.regr <- plgraphics::plregr
## envirnment(plot.regr) <- environment(plgraphics::plregr)
##  ===================================================================
print.regr <-
  function (x, call=TRUE, correlation = FALSE,
    termeffects = i.getoption("show.termeffects"),
    termcolumns = i.getoption("termcolumns"),
    termeffcolumns = i.getoption("termeffcolumns"),
    coefcolumns = i.getoption("coefcolumns"),
    digits = i.getoption("digits"), 
    symbolic.cor = p > 4, signif.stars = getOption("show.signif.stars"),
    na.print = i.getoption("na.print"),
    residuals=FALSE, niterations=FALSE, ...)
{
  ##
  if (is.null(na.print)) na.print <- "."
  ## doc
  ldoc <- i.getoption("doc")
  if (length(ldoc)==0) ldoc <- 1
  if (ldoc>=1) if (length(tit(x)))
    cat("\n ",tit(x),"\n")
  if (ldoc>=2) if (length(doc(x)))
    cat(" ",paste(doc(x),"\n "))
  ## mlm
  if (inherits(x,"mlm"))
    return(invisible(print.mregr(x, na.print=na.print, ...)))
  ## preparation
  lItermeff <- i.def(termeffects, TRUE)
  ## call, fitting fn, residuals
  if (call) {
    if(!is.null(x$call)) {
      cat("\nCall:\n")
      cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),"\n", sep = "")
    }
    cat("Fitting function: ",x$fitfun,
        if (length(lfam <- x$family))
          paste("  Family:",
                if(inherits(lfam, "family"))
                  paste(lfam$family,
                        if (length(llink <- lfam$link))
                          paste("  Link:",llink)
                        ) else lfam
                ), "\n")
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
      print(rq, digits = digits, na.print=na.print, ...)
    } else {
      if (rdf > 0) print(resid, digits = digits, na.print=na.print, ...)
      else  cat("ALL", df[1],
                "residuals are 0: no residual degrees of freedom!\n")
    }
  }
  ##
  nsingular <- df[3] - df[1]
  if ((!is.na(nsingular))&&nsingular>0)
    cat("\nCoefficients: (", nsingular,
        " not defined because of singularities)\n", sep = "")
  if (!is.null(x$sigma))
    if((!is.finite(x$sigma))||x$sigma<=0)
      cat("\n!!! Error variance is 0 !!!")
  ## termtable
  lttab <- x$termtable
  if (length(lttab)>0) {
    if (inherits(lttab, "termtable")) {
      lIttab <- TRUE
      if (signif.stars) termcolumns <- union(termcolumns, "p.symbol")
      if(!is.null(termcolumns)) {
        if (all(termcolumns=="")) lIttab <- FALSE else {
          ljp <- match(termcolumns,colnames(lttab), nomatch=0)
          if (sum(ljp)!=0) 
            ##        warning(":print.regr: no valid columns of  termtable  selected") else
            lttab <- lttab[,ljp,drop=FALSE]
        }
      }
      if (lIttab) {
        cat("\nTerms:\n")
        ## round R2.x, signif, p.value
        ljrp <- colnames(lttab)[dropNA(pmatch(c("R2","signif","p.v"),
                                             colnames(lttab)))]
        if (length(ljrp))
          lttab[,ljrp] <- round(as.matrix(lttab[,ljrp]),max(3,digits))
        if ("signif"%in%ljrp) lttab$signif <- round(lttab$signif,last(digits)-1)
        lttabf <- format(lttab, na.encode=FALSE)
        lttabp <- data.frame(lapply(lttabf, function(x) sub("NA",na.print,x)),
                             row.names=row.names(lttab))
        print(lttabp, quote=FALSE, na.print=na.print)
        if (signif.stars>=1) 
          cat("---\nSignif. codes:  ", attr(x$termtable, "legend"),"\n", sep = "")
      } ## end if(lIttab)
    } else print(lttab, digits=digits, na.print=na.print) ## !inherits(.,"termtable")
  ## --- error block
  } else {
    if (length(x$coeftable)) {
      cat("\nCoefficients:\n")
      lcftb <- x$coeftable
      if("coef"%in%coefcolumns & "estimate"%in%names(lcftb))
        names(lcftb)[match("estimate", names(lcftb))] <- "coef"
      print(lcftb[,coefcolumns], na.print=na.print, digits=digits)
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
        if (length(lr2a <- x$adj.r.squared))
          paste("   Adjusted R^2: ", formatC(lr2a, digits = digits)),
        if (length(lAIC <- x$AIC)&&!is.na(lAIC))
          paste("   AIC: ", formatC(lAIC, digits = log10(abs(lAIC))+3),"\n" ) )
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
    cat ("\n")
    if (length(x$devtable)) print(x$devtable, na.print=na.print, digits=digits)
    cat("\nDistribution: ",x$distrname)
    if (length(x$dispersion))
      cat(".  Dispersion parameter: ",
          if ((!is.null(attr(x$dispersion,"fixed")))&&
              attr(x$dispersion,"fixed"))
             "fixed at ", format(x$dispersion))
    else {
      if (length(x$scale))
        cat(".  Shape parameter ('scale'): ",
            if ((!is.null(attr(x$scale,"fixed")))&&
                attr(x$scale,"fixed"))
              "fixed at ", format(x$scale))
      else cat("\n")
    }
  }
  ## --- additional coefficients
  if (x$distrname=="multinomial") {
    cat("\nCoefficients:\n")
    print(t(x$coefficients), na.print=na.print)
  } else {
    if (length(lttab)&lItermeff) {
      if (lItermeff==1) {
##-         lidf <- match("df",colnames(x$termtable))
##-         if (is.na(lidf)) {
##-           if (getOption("verbose"))
##-             warning(":print.regr: df of coef not available")
##-         } else { ## dummy coefficients
        mterms <- row.names(x$termtable)[is.na(x$termtable[,1])]
        mt <- if (length(mterms)>0 & length(x$termeffects)>0) {
                imt <- mterms%in%names(x$termeffects)
                x$termeffects[mterms[imt]]
              } else NULL
      } else mt <- x$termeffects
      if (length(mt)>0) {
        cat("\nEffects of factor levels:\n")
        print.termeffects(mt, digits=digits, na.print=na.print,
                          columns=i.getoption("termeffcolumns")) }
    } ## else  cat("\n")
  }
  if (length(x$n.censored)) {
    lnc <- 100*x$n.censored/x$n.obs
    cat(paste("\ncensored         ",
              paste(paste(names(lnc), "=", round(lnc,1), "%"), collapse=" ;  ")))
    if (length(x$n.fitout)) {
      lnf <- 100*x$n.fitout/x$n.obs
      cat(paste("\nfit outside limit",
        paste(paste(names(lnf), "=", round(lnf,1), "%"), collapse=" ;  ")))
    }
  }
  ## cat("\nAIC: ", format(x$AIC, digits = digits + 1), "\n", sep = "  ")
  if (niterations&&length(x$iter)>0)
    cat("Number of iterations:", x$iter, "\n")
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
        print(correl[-1, -p, drop = FALSE], digits = digits,
              na.print = na.print)
      }
    }
  }
  cat("\n")
  invisible(x)
}
## ==========================================================================
## currently only called from print.regr():
print.mregr <- function (x, na.print=i.getoption("na.print"), ...)
{
  ## Purpose:   collect results for mregr object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: Feb 2008
  f.prv <- function(x) paste(paste(names(x),x,sep=" = "),collapse=", ")
  if (is.null(na.print)) na.print <- "."
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),"\n", sep = "")
  cat("\nCoefficients:\n")
  print(x$coefficients, na.print=na.print)
  cat("\nP-values:\n")
  print(round(x$pvalues,4), na.print=na.print)
  cat("\nStatistics for individual response models:\n")
  print(x$stats)
  cat("\nResidual degrees of freedom: ",x$df,"\n")
  ldr <- x$drop1
  if (!is.null(ldr)) {
    cat("\nMultivariate tests for all responses\n  Degrees of freedom: ",
        f.prv(attr(ldr,"df")),"\n")
    print(ldr[,], na.print=na.print)
  }
  invisible(x)
}
## =====================================================================
residuals.regr <- function (object, type=NULL, standardized=FALSE, ...)
{
##-   if (!is.na(pmatch(type,"condquant"))) {
##-     ## this seems to apply only if residuals.regr is called explicitly
  lcall <- match.call()
  lcall$type <- if (is.null(type)||is.na(type)) NULL else type
  lff <- object$fitfun
  if (lff=="glm" && !is.null(type) && substr(type,1,4)=="cond") lff <- "polr"
  lcall[[1]] <-
    switch(as.character(lff),
           "polr" = quote(residuals.regrpolr),
           "survreg" = quote(residuals.regrsurvreg),
           "coxph" = quote(residuals.regrcoxph),
           quote(residuals))
  class(object) <- setdiff(class(object), "regr")
  lcall$object <- object
  rr <- structure( eval(lcall, envir=parent.frame()), type=type)
  if (standardized) {
    lstr <- stdresiduals(object, residuals=rr)
    lat <- attributes(lstr)
    attributes(lstr)[c("leverage","stdresratio","stddev")] <- NULL
    attr(rr,"stdresiduals") <- lstr
    attr(rr,"leverage") <- lat$leverage
    attr(rr,"stdresreatio") <- lat$stdresratio
  }
  rr
}
## ===================================================================
fitted.regr <-
  function (object, type=NULL, ...)
{
  if (is.null(type)&&pmatch("fitted",names(object),nomatch=0))
    return( naresid(object$na.action, object$fitted) )
  lres <- object$residuals
  if (inherits(lres, "condquant"))
    structure(lres[,"fit"], names=row.names(lres))
  else  {
    class(object) <- c(paste("regr",class(object)[1], sep=""),
                       setdiff(class(object), "regr"))
    predict(object, type=type, ...)
  }
}
## ==========================================================================
predict.regr <-
  function (object, newdata = NULL, scale = object$sigma,
            df=object$df.residual, type = NULL, ...)
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
  if (is.null(ldt)) return(
    predict(object, type=type, scale=object$sigma,
            dispersion=object$dispersion^2, ... ) )
    ## analyze variables
##  if (!is.null(ldt)) {
    ## terms with logst -> need original thresholds
  ltl <- attr(terms(object),"term.labels")
  ltll <- grep("logst\\(",ltl)
  lvlogst <- NULL
  if (length(ltll)) {
    lvlogst <- unique( gsub(".*logst\\((.*)\\).*","\\1", ltl[ltll]) )
    lmodel <- model.frame(object)
    lform <- as.character(as.expression(formula(object)))
  }
  for (lvn in names(ldt)) {
    lv <- ldt[[lvn]]
    ## factors
    if (is.factor(lv))  ldt[[lvn]] <- 
      if (match(lvn,names(object$binlevels),nomatch=0)>0) ## binary
        match(lv,object$binlevels[[lvn]])-1
      else  factor(lv)
    ## logst
    if (lvn %in% lvlogst) {
      lt <- ltl[grep(lvn, ltl)[1]]
      lvv <- lmodel[[lt]]
      lth <- attr(lvv, "threshold")
      if (is.null(lth))
        stop("!predict.regr! variable in term ",lt,
             "  not found in model.frame or threshold not available. \n",
             "  Prediction with 'logst' would fail.",
             "  Store transformed variable in data.frame")
      lth <- round(lth,5) ## get it from there,
        ## since it may have been set by the call to regr
      lform <- gsub(paste("logst *\\(",lvn,"\\)",sep=" *"),
                   paste("logst(",lvn,", threshold=",lth,")",sep=""), lform)
    }
  }
  if (length(lvlogst))
    object$terms <- terms(as.formula(lform), data=object$data)
  predict(object, newdata=ldt, type=type, scale=object$sigma, 
          df=df, dispersion=object$dispersion^2, ... )
}
## ==========================================================================
confint.regr <- function (fitted, ...)
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
##- extractAIC.regr <- function (fit, scale = 0, k = 2, ...)
##- {  ##- #  AIC, divided by n to allow for comparing models with different n
##- ##-   lres <- fit$residuals
##- ##-   if (is.null(lres)) {
##- ##-     lfit <- fit$fitted
##- ##-     fit$fitted.values <- dropNA(lfit)
##- ##-   } else  fit$residuals <- dropNA(lres)
##-   class(fit) <- setdiff(class(fit),"regr")
##-   extractAIC(fit, scale = scale, k = k, ...)
##- }
## -------------------------------------------------------------------------
vcov.regr <- function (object, ...) {
  cov <- object$covariance 
  if (is.null(cov)) {
    class(object) <- setdiff(class(object),"regr")
    vcov(object)
  }
  cov
}
## --------------------------------------------------------------------
df.residual.regr <- function (object, ...) {
  df <- object$df.residual
  if (is.null(df)) df <- object$df
  if (length(df)>=2) df <- df[2]
  if (is.null(df)) {
    sry <- summary(object)
    df <- sry$df
    if (is.null(df))
      df <- NROW(object$residual)-NROW(object$coefficients)
  }
  if (is.null(df)) df <- Inf
  df
}
## ------------------------------------------------------------------------
vif.regr <- function (mod, cov, mmat)
{
  ## Purpose:   vif.lm  of library  car
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
    else if (mod$fitfun%nin%c("polr","coxph","survreg"))
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
## ==========================================================================
nobs.survreg <- function (object, use.fallback = TRUE) {
  lnobs <- length(object$linear.predictors)
  if (lnobs==0) lnobs <- NROW(residuals(object))
  lnobs
}
nobs.coxph <- function (object, use.fallback = TRUE) {
  object$n
}
## ==========================================================================
createNAvars <-
  function (data, vars=NULL, na.prop=0.1, na.label=".NA.",
           na.values=NULL, name.suffix=c(".X",".NA"), append=TRUE, ...)
{
  if (is.matrix(data)) data <- as.data.frame(data)
  if (!(is.data.frame(data)))
    stop("!createNAvars! unsuitable first argument")
  ldt <- if (length(vars)) {
           if (length(lwr <- setdiff(vars, colnames(data))))
             stop("!createNAvars! Variables  ",paste(lwr, collapse=", "),
                  "  not in 'data'.")
           data[,vars]} else data
  if ((!is.numeric(na.prop))||na.prop>=1)
    stop("!createNAvars! unsuitable argument 'na.prop'")
  lvna <- sumNA(data) > max(na.prop,0) * nrow(data)
  ldclass <- sapply(data, is.numeric) + 2*sapply(data, is.factor)
  if (any(lwr <- ldclass==0))
    stop("!createNAvars! Funny variables  ", paste(lwr, collapse=", "))
  if (any(lnum <- lvna & ldclass==1)) 
    rrx <- xNA(ldt[,lnum, drop=FALSE], na.values=na.values,
              name.suffix=name.suffix)
  if (any(lfac <- lvna & ldclass==2)) {
    rrf <- factorNA(ldt[,lfac, drop=FALSE],na.label=na.label,
                ...)
  }
  if (!any(c(lnum,lfac))) {
    warning(":createNAvars: no variables had enough NAs to be modified")
    return( if(append) data else NULL )
  }
  if (append) { if (any(lfac)) data[names(rrf)] <- rrf
    if (any(lnum)) data.frame(data, rrx) else data
  }
  else data.frame(rrf,rrx)
}
## --------------------------------------------------------------------  
factorNA <- function (data, na.label=".NA.", na.prop=0, ...)
{
  if (missing(data)||length(data)==0)
    stop("!xNA! Argument 'data' missing, with no default, or NULL")
  data <- data.frame(data)
##  if (!is.data.frame(data)) stop("!xNA! Argument 'data' ...")
  lnalabel <- as.character(na.label)
  lvn <- names(data)[sapply(data,is.factor)] 
  ldt <- data[,lvn, drop=FALSE]
  lnanum <- na.prop*nrow(data)
  for (lv in seq_along(lvn)) {
    lfac <- factor(data[,lvn[lv]], ...)
    lna <- is.na(lfac)
    if (sum(lna)>lnanum) {
      levels(lfac) <- c(levels(lfac), lnalabel)
      lfac[lna] <- lnalabel
      ldt[,lv] <- lfac
    }
  }
  structure(ldt, NA.label = na.label)
}
## -------------------------------------------------------------------
xNA <-
  function (data, na.values=NULL, na.prop=0.1, name.suffix=c(".X",".NA"))
{
  if (missing(data)||length(data)==0)
    stop("!xNA! Argument 'data' missing, with no default, or NULL")
  data <- data.frame(data)
##  if (!is.data.frame(data)) stop("!xNA! Argument 'data' ...")
  lvn <- names(data)
  lnsuff <- if (length(name.suffix)==1) c(name.suffix,".NA") else
    c(name.suffix,".X",".NA")[1:2]
  lvnx <- paste(lvn, lnsuff[1], sep="")
  lvnna <- paste(lvn, lnsuff[2], sep="")
  na.values <-
    if (is.null(na.values)) apply(data,2,median,na.rm=TRUE)  else
      rep(na.values, length=length(lvn))
  names(na.values) <- lvn
  ldt <- data[,1,drop=FALSE]
  names(ldt) <- "."
  for (lv in seq_along(lvn)) {
    lna <- is.na(data[,lv])
    if (any(lna)) {
      ldt[,lvnx[lv]] <-
        structure(ifelse(lna, na.values[lv], data[,lv]), na.value=na.values[lv])
      ldt[,lvnna[lv]] <- lna
    } else ldt[,lvn[lv]] <- data[,lv]
  }
  ##structure(ldt[,-1, drop=FALSE], xNA.values = na.values)
  ldt[, -1, drop=FALSE]
}
## ===========================================================================
termeffects <-
  function (object, se = 2, df = df.residual(object), ...)
    ## --------------------------------------------------------------
{
  if (is.atomic(object)||is.null(terms(object)))
      stop("!termeffects! inadequate first argument")
 ##  xl <- object$xlevels
  Terms <- delete.response(terms(object))
  tl <- attr(Terms, "term.labels")
  dcl <- attr(Terms,"dataClasses")[-1]
  if (all(dcl=="numeric")) 
    return(as.list(coef(object)))
  ## result already available?
  allc <- object$termeffects
  if ((!is.null(allc))&&length(allc)==length(tl)&&
      (is.matrix(allc[[length(allc)]])|!se)) return(allc) ## !!! check!
  int <- attr(Terms, "intercept")
  facs <- attr(Terms, "factors")
  mf <- object$model  ##! d.c used all.vars
  if (is.null(mf)) mf <- model.frame(object)
  xtnm <- dimnames(facs)[[1]]  ## names  ##! replaces vars
  xtlv <- lapply(mf[,xtnm, drop=FALSE],function(x) levels(x)) ## levels
  lcontr <- object$contrasts
  imat <- which(substr(dcl,1,7)=="nmatrix") ## resulting from bs()
  if (length(imat)) {
    xtlv[imat] <-
        lapply(as.list(dcl[imat]),
               function(x) as.character(1:as.numeric(substr(x,9,12))))
    ##    lcontr <- c(lcontr, structure(rep(contr.id,length(tl)), names=tl)[imat])
    lctr <- list()
    for (li in seq_along(imat))
      lctr <- c(lctr, list(diag(length(xtlv[[li]]))))
    names(lctr) <- names(dcl)[imat]
    lcontr <- c(lcontr, lctr)
  }
  xtnl <- pmax(sapply(xtlv,length),1) ## number of levels
  termnl <- apply(facs, 2L, function(x) prod(xtnl[x > 0])) ##! lterms
  nl <- sum(termnl)
## --- df.dummy: data frame of simple terms
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
      rnn[pos + 1L:termnl[j]] <-
        apply(as.matrix(tmp), 1L, function(x) paste(x, collapse = ":"))
    }
    pos <- pos + termnl[j]
  }
  ## attributes
  attr(df.dummy,"terms") <- attr(mf,"terms")
  lci <- sapply(df.dummy,is.factor)
  lcontr <- lcontr[names(lci)[lci]] ## factors with 1 level have disappeared (?)
  if (lIpolr <- inherits(object, "polr")) {
    attr(Terms, "intercept") <- 1
    mm <- model.matrix(Terms, df.dummy, contrasts.arg=lcontr, xlev=xtlv)
    asgn <- attr(mm, "assign")[-1]
    mm <- mm[,-1]
  } else {
    mm <- model.matrix(Terms, df.dummy, contrasts.arg=lcontr, xlev=xtlv)
    asgn <- attr(mm, "assign")
  } 
  if (anyNA(mm)) {
    warning("some terms will have NAs due to the limits of the method")
    mm[is.na(mm)] <- NA
  }
  ## calculate dummy coefs
  coef <- object$coefficients ##!!! cf <-
##-   if (!use.na) 
  ##-     coef[is.na(coef)] <- 0
  lnna <- is.finite(coef)
  names(asgn) <- colnames(mm)
  if (any(!lnna)){
    coef <- coef[lnna]
    mm <- mm[,lnna]
    asgn <- asgn[lnna]
  }
  if (se) {
    cov <- vcov(object)
    if (is.null(cov)) {
      warning(":termeffects: no covariance matrix of coefficients found.",
              " Returning coefficients only")
      se <- FALSE
    } else
      if (inherits(object, "polr"))
        cov <- cov[1:length(coef),1:length(coef)]
  }
  licf <- pmatch(colnames(mm), names(coef))
##  asgn <- asgn[names(coef)] ## !!!
  res <- setNames(vector("list", length(tl)), tl)
  ljfail <- NULL
  for (j in seq_along(tl)) {
    mmr <- rn == tl[j]  ## rows corresponding to the term
    mmc <- asgn==j ## & !is.na(coef)
    lcf <- coef[licf[mmc]]
##    mmc <- names(asgn)[asgn == j & !is.na(coef)]  ## columns (logical fails for polr, vcov() too large) !!! was  which
    mmpart <- mm[mmr, mmc, drop=FALSE]
    rrj <- setNames(drop(mmpart %*% lcf), rnn[mmr]) ## coef[mmc]
    if (se) {
      if (any(is.na(rrj))) {
        warning(":termeffects: missing coef for term '", tl[j],
                "'. no standard errors etc")
        ljfail <- c(ljfail, tl[j])
      } else {
        sej <- sqrt(diag(mmpart %*% cov[mmc,mmc] %*% t(mmpart)))
        rrj <- ciSignif(rrj, sej, df) ## already calculated in  coeftable
      }
    }
    res[[j]] <- rrj
  }
  if (length(ljfail))
    warning(":termeffects: error calculating se for terms ",
            paste(ljfail, collapse=", "))
  if (int > 0) {
    res <- c(list(`(Intercept)` = coef[int]), res)
  }
  if (lIpolr) 
    res <- c(res,list("(Intercepts)"=ciSignif(object$intercepts[,1],
                                         object$intercepts[,2], df) ))
  ##  class(res) <- "termeffects" ## don't do that:
  ##                                 want to be able to print the whole table
  res
}
## --------------------------------------------------------------------
print.termeffects <- function (x, columns=NULL, transpose=FALSE, ...)
{
  if (is.null(columns)) columns <- "all"
  columns[columns=="estimate"] <- "coef"
  csymb <- "coefsymb"%in%columns
  if ("all"%in%columns)  columns <-
      if(csymb)
        c("coefsymb", "se", "ciLow", "ciUp", "testst",
          "signif", "p.value") else
        c("coef", "se", "ciLow", "ciUp", "testst",
          "signif", "p.value")
  for (li in seq_along(x)) {
    xi <- x[[li]]
    if (is.null(dim(xi))) next
    if (csymb)
      xi$coefsymb <-
        if ("p.symbol"%in%names(xi)) {
          lps <- as.character(xi[,"p.symbol"])
          lps[is.na(lps)] <- "   "  ##  misleading?
          paste(format(xi[,1],...), lps)
        } else  xi[,1]
    xif <- format(xi[,intersect(columns,names(xi)), drop=FALSE],...)
    xif <- if (ncol(xif)==1 || (nrow(xif)>1 & transpose)) t(xif) else xif
    if (nrow(xif)==1) row.names(xif) <- " "  ## drop row name
    if (ncol(xif)==1) colnames(xif) <- " "  ## drop col name
    if (prod(dim(xif))==1) xif <- as.character(xif[1,1])
    x[li] <- list(xif)
  }
  print(unclass(x), quote=FALSE, ...)
}
## ====================================================================
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
  } else {
##-     if (!is.character(scope))
##-             scope <- attr(terms(update.formula(object, scope)),
##-                 "term.labels")
    if (is.character(scope))
      as.formula(paste("~",paste(scope, collapse="+"))) else scope
  }
  if (length(scope)==0) {  ## || !is.formula(scope) ## drop.scope is character
    warning(":drop1/add1.regr: no valid scope")
    ldr1 <- data.frame(Df = NA, "Sum of Sq" = NA, RSS =NA, AIC = NA,
                      row.names = "<none>")
    return(ldr1)
  }
  ##  lform <- update(formula(object), scope !!!)
  ## !!! model.frame for finding the valid rows
  class(object) <- setdiff(class(object), "regr")
  fcall <- object$funcall
  if (!is.null(fcall)) object$call <- fcall
##-   dfm <-  object$df
  ##-   object$df <- dfm[setdiff(names(dfm),"residual")] ## survreg !
##-   if (inherits(object, c("survreg","coxph")))
##-     object$df <- object$df["original"]
##
  dr1 <- if (add) { ## ------------ add
    if (class(object)[1]=="lmrob")
        stop("!add1.regr! 'add1' not (yet) available for 'lmrob' objects")
    ldata <- eval(object$call$data, envir=environment(formula(object)) )
    lres <- object$residuals
    li <-
      if (length(lres)) row.names(ldata)%in%rownames(cbind(lres)) else rep(TRUE, NROW(ldata))
    if (length(ldata[li,])==0) stop("!drop1.regr! no data found ")
    lvars <-unique(c(all.vars(formula(object)),
                     if (is.formula(scope)) all.vars(scope) else scope))
    lvars <- lvars[lvars%in%names(ldata)]
    linotNA <- li & !apply(is.na(ldata[,lvars]),1,any)
    lnobs <- sum(linotNA)
    lnrd <- sum(li)
    lfc <- object$funcall ## the call to the effective R function
    if (is.null(lfc)) lfc <- object$call
    if (lnobs!= lnrd) {
      notice("add1.regr: refitting object to ",lnobs," / ",lnrd,
              " observations due to missing values")
      if(!is.null(lsubs <-
                    eval(lfc$subset, envir=environment(formula(object))))) {
        lnsubs <- rep(FALSE,length(linotNA))
        lnsubs[lsubs] <- TRUE
        linotNA <- linotNA &!lnsubs
      }
      lfc$subset <- linotNA
      object <- eval(lfc, envir=environment(formula(object)))
##-       object$call[[1]] <-
##-         if (is.null(lfc)) as.name(class(object)[1]) else
##-            lfc[[1]]
##-       object <- update(object, subset=linotNA)
#      environment(object$call$formula) <- environment()
##-       class(object) <- setdiff(class(object), "regr")
    }
    if (!all(linotNA)) { ## needed if NA's have been generated by transformations
      lfc$subset <- linotNA
      object <- eval(lfc, envir=environment(formula(object)))
    }
    add1(object, scope=scope, scale=scale, test=test, k=k, ...)
  } else {  ## ---------------------  drop
    if (class(object)[1]%in%c("lmrob")) ## to be expanded
       drop1Wald(object, test="F", ...) else {
    ldata <- object$allvars # eval(object$call$data)
    if (is.null(ldata)) stop("!drop1.regr! no data found ")
    ## all predictors must get the same missing observations
    lina <- apply(is.na(ldata),1,any)  
    if (any(lina)) ldata[lina,] <- NA
    object$call$data <- ldata
    drop1(object, scope=scope, scale=scale, test=test, k=k, ...)
    }
  }
##-   rnm <- row.names(dr1)
##-   row.names(dr1) <- paste(ifelse(substring(rnm,1,1)=="<","",
##-                                  if (add) "+ " else "- "),rnm,sep="")
  attr(dr1,"drop") <- !add
##-   if(add) attr(dr1,"ndropped") <- lndiff
  if (sorted) {
    lsrt <- dropNA(match(c("AIC","p.value"),colnames(dr1)))
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
           k = 2, ...) 
{
    x <- model.matrix(object)
    offset <- model.offset(model.frame(object))
    n <- nrow(x)
    asgn <- attr(x, "assign")
    lterms <- terms(object)
    tl <- attr(lterms, "term.labels")
    attr(lterms, "order") <-  rep(1,length(tl))
    if (is.null(scope)) 
      scope <- tl # drop.scope(lterms)
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
    lsig <- c(object$sigma, object$scale)[1]
    chisq <- lsig^2 * rdf
    ## sum(weighted.residuals(object)^2, na.rm = TRUE)
    ## deviance.lm(object)
    dfs <- numeric(ns)
    RSS <- numeric(ns)
    cov <- object$cov.unscaled
    if (is.null(cov)) cov <- object$covariance/lsig^2
    if (length(cov)==0) stop("!drop1Wald! no covariance matrix found")
    cf <- object$coefficients
    ##- jj <- match(names(cf),colnames(cov), nomatch=0)
    ##-     if (!(any(jj==0)&&all(is.na(cf[jj==0]))))
    jj <- match(colnames(cov),names(cf), nomatch=0)
    if (any(jj==0))
      warning(":drop1Wald: coefficient(s) and cov. matrix may not correspond")
    coef <- cf[jj]
    asgn <- asgn[jj]
    if (any(names(coef[!is.na(coef)])%nin%names(coef)))
      stop("!drop1Wald! coefficient(s) not appearing in covariance matrix")
##-     y <- object$residuals + predict(object)
    for (i in 1:ns) {
        ii <- which(asgn==ndrop[i]) ## seq_along(asgn)[asgn == ndrop[i]]
        RSS[i] <- if (length(ii)==1) coef[ii]^2/cov[ii,ii] else
          coef[ii]%*%solve(cov[ii,ii])%*%coef[ii]  ## !!! REPLACE THIS
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
    dfs <- c(c(object$rank,object$df)[1], dfs)
    RSS <- chisq + c(0, RSS)
    if (scale > 0) 
        AIC <- RSS/scale - n + k * dfs
    else AIC <- n * log(RSS/n) + k * dfs
##-     dfs <- dfs[1] - dfs
##-     dfs[1] <- NA
    aod <- data.frame(Df = dfs, "Sum of Sq" = c(NA, RSS[-1] - 
        RSS[1]), RSS = RSS, AIC = AIC, row.names = scope, check.names = FALSE)
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
##- step <- function (object, ...)
##-   UseMethod("step")
##- 
##- step.default <- stats::step
## step.default <- get("step", pos="package:stats")
#### !!! sollte das anders heissen? step.default <- stats::step  ???

step.regr <- function (object, scope=NULL, expand=FALSE, scale = 0,
  direction = c("both", "backward", "forward"), trace = FALSE, keep = NULL,
  steps = 1000, k = 2, ...)
{
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
  } ## end step.results
  ## ---
  Terms <- terms(object)
  object$call$formula <- object$formula <- Terms
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  ## scope
  lform <- formula(object)
  if (is.null(scope))  ## !! was missing(scope)
    scope <- if (expand) terms2order(object) else lform
  if (!is.list(scope)) {
    if (is.character(scope))
      scope <- as.formula(paste("~.+",paste(scope, collapse="+")))
    scope <- list(lower=NULL, upper=scope)
  }
  if (is.null(names(scope))&length(scope)==2)
    names(scope) <- c("lower","upper")
  fdrop <- if (length(fdrop <- scope$lower)) {
             attr(terms(update.formula(lform, fdrop)), "factors")
           } else numeric()
  if (length(fadd <- scope$upper)) {
    lform <- update.formula(lform, fadd)
    fadd <- attr(terms(lform), "factors")
  } 
  ## data
  lcalldata <- object$call$data
  ldata <- object$allvars 
  lvars <- all.vars(lform)
  if (any(lvars%nin%names(ldata)))
    ldata <- eval(object$call$data, envir=environment(formula(object)) )
  if (any(li <- lvars%nin%names(ldata)))
    stop("!step.regr! variable  ", paste(lvars[li], collapse=", "), "  not available")
  stepdata <- na.omit(ldata[,lvars])
  object$call$data <- stepdata
  if (length(object$funcall)) object$funcall$data <- stepdata
  models <- vector("list", steps)
  if (!is.null(keep))
    keep.list <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  if (inherits(fit, "survreg")) fit$residuals <- NULL
  bAIC <- extractAIC(fit, scale, k = k, ...)
  edf <- bAIC[1L]
  bAIC <- bAIC[2L]
  if (is.na(bAIC))
    stop("AIC is not defined for this model, so `step` cannot proceed")
  nm <- 1
  if (trace) {
    cat("Start:  AIC=", format(round(bAIC, 2)), "\n",
        cut.string(deparse(formula(fit))), "\n\n", sep = "")
    utils::flush.console()
  }
  models[[nm]] <-
    list(deviance = mydeviance(fit), df.resid = n - edf, change = "",
         AIC = bAIC)
  if (!is.null(keep))
    keep.list[[nm]] <- keep(fit, bAIC)
  usingCp <- FALSE
  ## ------------------------
  lrgopt <- regroptions(notices=FALSE)
  on.exit(regroptions(notices=attr(lrgopt, "old")))
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    ## backward
    if (backward && length(scope$drop)) {
      aod <- drop1(fit, scope$drop, scale = scale, trace = trace,
                   k = k, sorted=FALSE, ...)
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
        aod <- if (is.null(aod)) aodf  else {
          names(aodf) <- names(aod)
          rbind(aod, aodf[-1, , drop = FALSE])
        }
      }
    }
    ## backward or forward
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
    if (o[1L] == 1)  break
    change <- rownames(aod)[o[1L]]
    ## update
    usingCp <- match("Cp", names(aod), 0L) > 0L
    ##    if (is.null(change))  break  else {
    fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit <- eval.parent(fit)
    if (inherits(fit, "survreg")) fit$residuals <- NULL
##-    nnew <- nobs(fit, use.fallback = TRUE)
##-     if (all(is.finite(c(n, nnew))) && nnew != n) {
##-       warning(":step.regr: number of rows in use has changed: \n  ",
##-               nnew," observations instead of ", n)
##-       n <- nnew
##-     }
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
    ##        if (bAIC >= AIC + 1e-07)
    ##  else   break
    nm <- nm + 1
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - edf,
                         change = change, AIC = bAIC)
    if (!is.null(keep))
      keep.list[[nm]] <- keep(fit, bAIC)
  }
  if (!is.null(keep))
    fit$keep <- re.arrange(keep.list[seq(nm)])
  lv <- all.vars(formula(fit))
  lnobs <- nrow(na.omit(ldata[,lvars]))
  lno <- NROW(na.omit(ldata[,lv]))
  if (lnobs<lno) {
    notice(" step.regr: Step worked only on ", lnobs,
           " observations, whereas result has ", lno,
           " -- due to missing values in dropped terms")
    fit <- update(fit, data=ldata)
  }
  fit$call$data <- lcalldata
  step.results(models = models[seq(nm)], fit, object, usingCp)
}
## ==========================================================================
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
    res <- object$resid
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
      bss <- ladd * (crossprod(lrg$resid)-rss)
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
  function (object, scope=NULL,
           test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"), ...)
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
    r <- cbind(object$resid) ## resid(object)
    if (nrow(r)!=nrow(m)) {
      notice(gettextf("add1! using the %d/%d rows from a combined fit",
                nrow(m), nrow(r)), domain = NA)
      lna <- !row.names(r)%in%row.names(m)
    }
    else lan <- NULL
  }
## ========================================================================
terms2order <- function (object, squared = TRUE, interactions = TRUE)
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
    notice("terms2order: nothing to add")
    return(formula(object))
  }
  ltadd <- paste("~.+", ltadd)
##-     attr(terms(update.formula(object, ltadd)), "term.labels")
  update.formula(object, ltadd)
}
## ==========================================================================
## ==========================================================================
compareTerms <-
  function (..., list=NULL, seq=NULL)
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
modelTable <-function (models, seq=NULL) 
{
  ## Purpose:   collect several models into a table
  ## -------------------------------------------------------------------------
  ## Arguments:
  ##   models   character vector of names of model objects
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  9 Aug 2000, 11:50
  if (inherits(models, regrModelClasses))
    stop("!modelTable! First argument is a single model.",
         " I need a list of at least two models")
  t.islist <- is.list(models)
  t.nm <- length(models)
  t.mnm <- names(models) ## if (t.islist) names(models) else models
  if (length(t.mnm)==0)
    if (t.islist) names(models) <- t.mnm <- paste("model",1:t.nm,sep=".")
    else t.mnm <- models
  if (t.islist) {
    if (any(t.wr <- sapply(models,
                           function(x) !inherits(x, regrModelClasses))))
      stop("!modelTable! element  ",paste(t.mnm[t.wr],collapse=", "),
           "  has inadequate class")
  } else  if(!is.character(models))
      stop("!modelTable! Argument 'models' must be a list or a character vector")
  t.ls <- t.trm <- t.cf <- setNames(vector("list",t.nm), t.mnm)
  t.nobs <- t.df <- t.sig <- t.fitfun <- setNames(rep(NA,t.nm), t.mnm)
  ## -----------------------------------------------------------
  for (li in 1:t.nm) {
    lr <- if (t.islist) models[[li]] else get(models[li],envir=parent.frame())
##-     lfitfun <- NULL
##-     for (lc in c("lm","lmrob","glm","multinom","polr","survreg"))
##-       if (inherits(lr,lc)) lfitfun <- lc
##-     if (is.null(lfitfun))
##-       stop(paste("!modelTable! Model ",li," is not an adequate model"))
    ##-     t.fitfun[li] <- lfitfun
    if (!inherits(lr, "regr"))
      stop("!modelTable! ... only programmed for 'regr' objects")
    t.fitfun[li] <- t.ff <- lr$fitfun
    t.nobs[li] <- lnr <-
      NROW(if(t.ff=="survreg") lr$linear.predictors else lr$fitted.values)
    t.df[li] <- ldf <- lnr-df.residual(lr)
    lt <- terms(lr)
    ltnm <- c( if(attr(lt,"intercept")) "(Intercept)", attr(lt, "term.labels"))
    t.cf[[li]] <-
      lr$termtable[match(ltnm,row.names(lr$termtable),nomatch=0),
                  c("coef","p.value")]
    t.trm[[li]] <- ltnm
##    t.trmc <- c(t.trmc, attr(terms(lr),"dataClasses"))
    lsig <- if (t.ff=="survreg") lr$scale else summary(lr)$sigma
    t.sig[li] <- c(lsig,NA)[1]
  }
  if (length(unique(t.nobs))>1)
    notice("modelTable: Models have different numbers of observations")
## --- collect
  t.tr <- unique(unlist(t.trm))
  ## --- coefs and p values
  t.nt <- length(t.tr)
  t.pr <- t.coef <- matrix(NA,t.nt,t.nm, dimnames=list(t.tr,t.mnm))
  ## ---
  for (li in t.mnm) {
    t.t <- t.trm[[li]]
    if (length(t.t)) {
      lcf <- t.cf[[li]]
      t.coef[t.t,li] <- lcf[,1]
      t.pr[t.t,li] <- lcf[,2]
    }
  }
  ## reorder
  if (length(seq)>0) {
    li <- match(seq, t.tr)
    t.t <- c(t.tr[dropNA(li)],t.tr[!t.tr%in%seq])
    t.coef <- t.coef[t.t,]
    t.pr <- t.pr[t.t,]
    t.sd <- t.sd[t.t]
  }
#  attr(t.coef,"standardized") <- t.trn
  if (all(is.na(t.sig))) t.sig <- NULL
  t.r <- list(coef=t.coef, p=t.pr, sigma=t.sig, nobs=t.nobs,
              df=t.df, fitfun=t.fitfun) # , sd.terms=t.sd
  class(t.r) <- "modelTable"
  t.r
}
## ==========================================================================
"[.modelTable" <- function (object,rows=NULL,cols=NULL, reduce=TRUE) {
  if (is.null(rows)) rows <- 1:nrow(object$coef)
  if (is.null(cols)) cols <- 1:ncol(object$coef)
  lp <- object$p[rows,cols,drop=FALSE]
  li <- if(reduce) !apply(is.na(lp),1,all)  else 1:nrow(lp)
  if (length(li)==0) stop("![.modelTable! no terms left")
  lsd <- if (length(object$sd.terms)) object$sd.terms[rows][li]
  lsig <- if (length(object$sigma)) object$sigma[li]
  rr <- list(coef=object$coef[rows,cols,drop=FALSE][li,,drop=FALSE],
             p=lp[li,,drop=FALSE], sigma=lsig,
             nobs=object$nobs[cols], df=object$df[cols],
             fitfun <- object$fitfun[cols] ) # ,sd.terms=lsd
#  attr(rr$coef,"standardized") <- attr(object$coef,"standardized")
  class(rr) <- class(object)
  rr
}
## ==========================================================================
format.modelTable <-
  function (x, digits=i.getoption("digits"),
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
  t.ntrm <- nrow(x$coef)
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
  t.cfo <- format(t.cf, digits=digits, ...)
  t.nna <- is.na(x$coef)&!is.na(x$p)
  t.cfo[1:t.ntrm,][t.nna] <- paste(c(rep(" ",1+digits/2),"+++"),collapse="")
  t.cfo[".df.",] <- paste("  ",format(t.cf[".df.",]))
  if (lnobs) t.cfo[".nobs.",] <- paste("",format(t.cf[".nobs.",]))
  lff <- x$fitfun
  if (length(lff))
    if (length(unique(lff))>1) {
      t.cfo <- rbind(t.cfo, fitfun=substr(lff,1,digits+3) )
      t.st <- rbind(t.st, .fitfun.="" )
    }
  t.cfo[grep("NA",t.cfo)] <- paste(c(rep(" ",2+digits/2),"-"),collapse="")
  t.ii <- (!is.na(t.cf))&t.cf==0
  t.cfo[t.ii] <- t.cf0 <- gsub("0"," ",t.cfo[t.ii])
  t.o <- paste(sep, format(t.cfo), sep, t.st)
  t.sdo <- if (lsd)
             paste(sep, sub("NA","  ", format(t.sd, digits=digits))) else NULL
  t.out <- cbind(t.sdo, matrix(t.o, nrow=nrow(t.cfo)))
  t.nm <- row.names(t.cfo)
  if (lsd) t.nm <- paste(t.nm,ifelse(!is.na(t.sd),"@",""))
  dimnames(t.out) <- list(t.nm,c(if (lsd) "sd",dimnames(x$p)[[2]]))
  structure( t.out, class=c("modelTable", "modelTF"), nterms=t.ntrm)
}
## ==========================================================================
print.modelTable <- function (x, tex = FALSE, transpose=FALSE, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   x        a modelTable object
  ##   tex      if TRUE, the output will be suitable for pasting into
  ##            (la)tex source
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 23 Dec 2008, 10:09
  t.out <- unclass( if (t.fo <- inherits(x, "modelTF")) x else
                   format(x, ...) )
  if (tex) {
    sep <- "&"
    if (t.fo) {
      t.a <- t.out[1,1]
      t.jsep <- substr(t.a,nchar(t.a),nchar(t.a))==sep
      if (!t.jsep) {
        warning(":print.modelTable: The header of the tex table will "
             ,"most probably be wrong")
        for (lj in 1:(ncol(t.out)-1)) t.out[,lj] <- paste(t.out[,lj], sep)
      }
    } else   t.out <- unclass( format(x, sep=sep, ...) )
    ##    if (transpose) t.out <- t(t.out) ## would not work yet, need to select rows
    mc1 <- "\\mc{2}{|c}{"
    mc2 <- "}"
    end <- "\\\\"
    headend <- "\\\\ \\hline"
    tabstart <- "  {\\def\\mc{\\multicolumn}\n\\begin{tabular}" # \\providecommand{
    tabend <- "  \\end{tabular}\n}"
    lnt <- attr(t.out,"nterms")
    t.end <- c(rep(end,nrow(t.out)-1),"")
    t.end[lnt] <- headend
    cat(tabstart, "{l", rep("|rl",ncol(t.out)),"}\n  ",
        paste(sep,mc1,colnames(t.out),mc2,sep=""),headend,"\n")
    t.rn <- format(row.names(t.out))
    if (t.fo) t.rn <- paste(t.rn,sep)
    for (li in 1:nrow(t.out)) cat(t.rn[li],t.out[li,],t.end[li],"\n")
    cat(tabend,"\n")
  } else {
    if (transpose) t.out <- t(t.out)
    colnames(t.out) <- paste("",colnames(t.out))
    print(structure(t.out, nterms=NULL), quote=FALSE)
  } ## !!!
  invisible(x)
}
## ==========================================================================
## additional useful functions
## ===========================================================================
## ======================================================================
quinterpol <- function (x, probs = c(0.25,0.5,0.75), extend=TRUE)
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
## =======================================================================
quantilew <- function (x, probs=c(0.25,0.5,0.75), weights=1, na.rm=FALSE)
{
  ## Purpose:   quantile with weights, crude version
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: KSK Projekt, Date: 14 Dec 1999, 12:02
  probs <- probs[!is.na(probs)]
  if (length(weights)==1) return(quantile(x, probs))
  if (length(weights)!=length(x))
    stop("!quantilew! lengths of 'x' and 'weights' must be equal")
  if (any(t.ina <- is.na(x))) {
    if (!na.rm) stop("!quantilew! NAs not allowed while 'na.rm' is FALSE")
    x <- x[!t.ina]
    weights <- weights[!t.ina]
  }
  t.i <- order(x)
  t.d <- x[t.i]
  t.wg <- cumsum(weights[t.i])/sum(weights)
  t.i1 <- apply(outer(t.wg,probs,"<"),2,sum)+1
  t.i2 <- pmin(apply(outer(t.wg,probs,"<="),2,sum)+1,length(t.d))
  (t.d[t.i1]+t.d[t.i2])/2
}
## ======================================================================
factor.na <- function (x, ordered=FALSE, naname="NA") {
  if (ordered) x <- ordered(x)
  if (is.ordered(x)) {
    levels(x) <- c(levels(x), naname)  
    x[is.na(x)] <- naname
    return(x)
  }
  x <- as.character(x)
  x[is.na(x)] <- naname
  factor(x)
}
## ---------------------------------------
##- quantNA <- function (vn, data, na.value=NULL) {
##-   if (missing(data)) stop("!quantNA! Argument 'data' missing, woth no default")
##-   dname <- as.character(substitute(data))
##-   if (is.character(vn)) {
##-     if (any(lvna <- vn%nin%names(data))) 
##-       stop("!quantNA! Variable(s) ",lv[lvna]," not in 'data")
##-     } else   vn <- names(data)[vn]
##-   na.value <- if (is.null(na.value))
##-                 apply(data[,vn, drop=FALSE],2,median,na.rm=TRUE)  else
##-     rep(na.value, length=length(lv))
##-   names(na.value) <- vn
##-   for (lv in vn) {
##-     lna <- is.na(data[,lv])
##-     if (any(lna)) {
##-       data[,paste(lv, ".NA", sep="")] <- lna
##-       data[lna, lv] <- na.value[lv]
##-     }
##-   }
##-   assign(dname, data, pos=1)
##- }
## =================================================================
## auxiliary functions
## ============================================================
factor2character <- function (x) {
  for (lj in 1:ncol(x))
    if (is.factor(x[,lj])) x[,lj] <- as.character(x[,lj])
  x
}
## ============================================================
regrModelClasses <- c("regr","lm","glm","survreg","coxph","rq","polr")
## ===========================================================================
## ===========================================================================
regrAllEqns <-
  function (formula, data, weights = NULL, nbest = 50, nvmax = 20,
           force.in = NULL, force.out = NULL, codes=NULL, really.big=FALSE,
           ...) 
{
  ## Purpose:   all subsets
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 14 Oct 2017, 09:26
  really.big <- really.big | nbest<=50
  lcall <- match.call()
  mm <- lcall[c("","formula","data","weights")]
  mm[[1]] <- as.name("model.frame")
  lcall$formula <- eval(lcall$formula)
  mf <- eval(mm, sys.frame(sys.parent()))
  x <- model.matrix(terms(formula, data = data), mf)
  lasgn <- attr(x,"assign")
  lxnm <- colnames(x)
  names(lasgn) <- lxnm
  ljint <- "(Intercept)"==lxnm[1]
  ltermnm <- c(if(ljint) "(Intercept)", attr(terms(mf),"term.labels"))
  y <- model.extract(mf, "response")
  wt <- model.extract(mf, "weights")
  if (is.null(wt))  wt <- rep(1, length(y))
  ## force.in and force.out: allow for formulas
  if (is.formula(force.in)) {
    if (any(lvi <- !(lv <- all.vars(force.in))%in%names(data)))
      stop("!regrAllEqns! Variable  ", paste(lv[lvi], collapse=", "),
           "  of argument 'force.in' not in 'data'")
    force.in <- setdiff(colnames(model.matrix(force.in, data)),"(Intercept)")
  }
  if(is.character(force.in)) {
    lwrong <- setdiff(force.in, lxnm)
    force.in <- which(lxnm%in%force.in)
  } else lwrong <- force.in%nin%(1:ncol(x))
  if (length(lwrong)) stop(paste("!regrAllEqns! Term ", lwrong,
                                 " of 'force.in' not in model"))
  if (is.formula(force.out)) {
    if (any(lvi <- !(lv <- all.vars(force.out))%in%names(data)))
      warning("!regrAllEqns! Variable  ", paste(lv[lvi], collapse=", "),
           "  of argument 'force.out' not in 'data'")
    force.out <- setdiff(colnames(model.matrix(force.out, data)),"(Intercept)")
  }
  if(is.character(force.out)) {
    lwrong <- setdiff(force.out, lxnm)
    force.out <- which(lxnm%in%force.out)
  } else lwrong <- force.out%nin%1:ncol(x)
  if (length(lwrong)) stop(paste("!regrAllEqns! Term ", lwrong,
                                 " of 'force.out' not in model"))
  if (ljint) x <- x[,-1]
  ## ---------------
  ls <- leaps::regsubsets(x, y, weights=wt, nbest=nbest, nvmax=nvmax,
                   force.in=force.in, force.out=force.out, int=ljint,
                   really.big=really.big, ...)
  lss <- summary(ls)
  lwhich <- lss$which
  ## factors: identify models that are unsuitable
  lbl <- unique(lasgn[duplicated(lasgn)])
  liok <- rep(TRUE, nrow(lwhich))
  ljok <- structure(rep(TRUE, ncol(lwhich)),names=colnames(lwhich))
  for (lk in lbl) {
    lj <- which(lasgn==lk)
    lwh <- lwhich[,lxnm[lj]]
    liok <- liok & ( apply(lwh,1,sum)%in%c(0,ncol(lwh)) )
    ljok[lxnm[lj[-1]]] <- FALSE
  }
  lwhs <- lwhich[liok,ljok,drop=FALSE]
  colnames(lwhs) <- lnm <- ltermnm[ljint+lasgn[colnames(lwhich)]][ljok]
  ## codes
  if (is.null(codes) || (length(codes)==1 && is.na(codes)))
    codes <- c(LETTERS,letters)
  if (length(names(codes))) {
    if (length(lwr <- setdiff(lnm,names(codes)))) {
      warning(":regrAllEqns: terms ", paste(lwr, collapse=", "),
              " not in 'names(codes)'. names are not used.")
      codes <- unname(codes)
    }
    codes <- codes[lnm]
  }
  if (is.null(names(codes))) 
    codes <- 
      structure(c(if(ljint) "1", rep(codes, length=ncol(lwhs)-ljint)),
                names=lnm)
  llb <- apply(lwhs,1, function(x) paste(codes[x], collapse="") )
  dimnames(lwhs) <- list(llb, lnm)
  lout <- lwhs
  lout[,] <- c(" ","*")[lwhs+1]
  ## criteria
  ldf <- apply(lwhich, 1, sum)
  lcr <- data.frame(df=ldf, lss[c("rsq","rss","adjr2","cp","bic")])
  lcrs <- lcr[liok,]
  dimnames(lcrs)[[1]] <- llb
  lall <- list(criteria=lcr, modsuit=liok, df=apply(lss$which, 1, sum),
               lss[c("which", "rsq", "rss", "adjr2", "cp", "bic", "outmat")])
  structure(
    list(which=lwhs, criteria=lcrs, codes=codes, force.in=force.in,
         force.out=force.out, call=lcall, outmat=lout, allsubsets=lall,
         obj=lss$obj),
    class="regrAllEqns")
}
## ------------------------------------------------------------
regrAllEqnsXtr <- function (object, nbest=1, criterion="cp")
{
  ## Author: Werner Stahel, Date: 18 Oct 2017, 17:37
  lwh <- object$which
  lwh <- lwh[order(object$criteria[,criterion])[1:nbest],,drop=F]
  rr <- apply(lwh,1, function(x)
    update(formula(object),
           as.formula(paste("~", paste(setdiff(colnames(lwh)[x], "(Intercept)"),
                                collapse="+"))) ) )
  if (nbest==1) structure(rr[[1]], modelcode=names(rr)) else rr
}
## ---------------------------------------------------------------------
print.regrAllEqns <-
  function (x, nbest=20, criterion="cp", printcriteria=FALSE, printcodes=TRUE,
           ...)
{
  ## Author: Werner Stahel, Date: 14 Oct 2017, 14:46
  li <- order(x$criteria[,criterion])[1:min(nbest,nrow(x$criteria))]
  lout <- x$outmat[li,]
  colnames(lout) <- x$codes[colnames(lout)]
  lout <- cbind(code=row.names(lout), df=as.character(x$criteria[li,"df"]), lout)
  row.names(lout) <- 1:nrow(lout)
  print(lout, quote=FALSE, ...)
  if (printcriteria) {
    lcr <- x$criteria[li,] 
    lcr <- cbind(code=row.names(lcr), lcr)
    row.names(lcr) <- 1:nrow(lcr)
    print(lcr, ...)
  }
  if (printcodes) print(cbind(code=x$codes), quote=FALSE, ...)
}
## ------------------------------------------------------------------
plot.regrAllEqns <-
  function (x, criterion="cp", critrange=10, minnumber=10, nbest=10,
           codes=x$codes, ncharhorizontal=6, col="blue",
           legend=TRUE, mar=6, main="", cex=0.7*par("cex"),
           cex.lab = par("cex.lab"), ...)
{
  ## Author: Werner Stahel, Date: 14 Oct 2017, 14:43
  lmod <- x$criteria[,c("df",criterion)]
  lcr <- lmod[,2]
  limod <- 1:nrow(lmod)
  if (is.null(critrange)) critrange <- Inf
  if ((!is.na(critrange)) && critrange>0)
    limod <- lcr<min(lcr)+critrange
  if (sum(limod)<minnumber) limod <- order(lcr)[1:min(minnumber,nrow(lmod))]
  lmod <- lmod[limod,]
  ldf <- lmod[,1]
  lcr <- lmod[,2]
  lwh <- x$which[limod,]
  llab <- codes[colnames(lwh)]
  lmar <- par("mar")
  if (length(mar)) {
    if (length(mar)!=4) {
      lmar <- c(lmar[1:2],mar[1],lmar[4])
    }
    oldpar <- par(mar=lmar)
    on.exit(par(oldpar))
  }
  plot(ldf,lmod[,2], type="n", xlab="df", ylab=criterion, main="", ...)
  if (length(main)) mtext(main, 3, lmar[3]-1.2, cex=1.2)
  for (ls in unique(ldf)) {
    lii <- which(ldf==ls)
    llcr <- lcr[lii]
    if (length(lii)) {
      lk <- lii[order(llcr)[1:min(length(lii),nbest)]]
      lvf <- apply(lwh[lk,,drop=F],2,all) ## in all models
      lvft <- paste(c(llab[lvf],"+"),collapse="")
      mtext(lvft, 3, 1, at=ls, cex=cex.lab, las=(nchar(lvft)>ncharhorizontal)+1,
            col=col)
      lwhs <- lwh[lk,!lvf, drop=F]
      llbv <- llab[!lvf]
      llb <- apply(lwhs,1,function(x) paste(llbv[x],collapse=""))
      llcr <- lcr[lk]
      if (length(lk)==1) text(ls, llcr, "+", cex=2*cex, col=col) else
        text(rep(ls,length(lk)),lcr[lk],llb, cex=cex, col=col)
    }
  }
  ##  text(lmod[,1],lmod[,2], row.names(lmod))
  if (length(legend)) {
    llab <- x$codes
    if (is.logical(legend)&&legend) {
      lcmin <- sapply(split(lcr,ldf), min)
      legend <- if (lcmin[1]>last(lcmin)) "bottomleft" else "bottomright"
    }
    if (is.character(legend)) {
      if (legend%nin%c("topleft","topright","bottomleft","bottomright"))
        legend <- "bottomright"
      legend(legend, paste(llab,names(llab)),pch=rep("",length=length(llab)))
    }
    if (is.numeric(legend)) {
      if (length(legend)==2)
        legend(legend[1],legend[2],
               paste(llab,names(llab)),pch=rep("",length=length(llab)))
      else warning(":plot.regrAllEqns! Argument 'legend' not suitable")
    }
  }
  invisible(lmod)
}
## ==================================================================
i.findformfac <- function (formula) {
  ## find variable involved in explicit factor terms in formula
  lfo <- format(formula)
  lmf <- c(gregexpr("(factor *\\([^)]*\\))", lfo),
           gregexpr("(ordered *\\([^)]*\\))", lfo) )
  lf <- function(x) 
    if(x[1]!=-1) substring(lfo, x, x+attr(x,"match.length"))
  all.vars(as.formula(
        paste("~",paste(unlist(lapply(lmf, lf)), collapse="+"))))
}
## ===================================================================
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
            notice("predict.mlm: Predictions on current data refer to _future_ responses\n")
        if (missing(newdata) && missing(weights)) {
            w <- weights(object) ## .default
            if (!is.null(w)) {
                weights <- w
                notice("predict.mlm: Assuming prediction variance inversely proportional to weights used for fitting\n")
            }
        }
        if (!missing(newdata) && missing(weights) && !is.null(object$weights) &&
##-             missing(pred.var)
            is.null(pred.var)
            )
            notice("predict.mlm: Assuming constant prediction variance even though model fit is weighted\n")
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
## repaired versions of drop functions
## ===========================================================================
add1.default <-
  function (object, scope, scale = 0, test=c("none", "Chisq"),
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
drop1.default <-
  function (object, scope, scale = 0, test=c("none", "Chisq"),
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
    n0 <- nobs(object)
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
      if (nobs(nfit) != n0)
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
fitted.polr <- function (object, type="link", na.action=object, ...) {
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
    lfit <- object$fitted
  if (type=="class")
    lfit <- factor(max.col(lfit), levels = seq_along(object$lev),
            labels = object$lev)
  naresid(object$na.action, lfit)
}
## ==========================================================================
## ==========================================================================
##- getRegroption <- function (x, default = NULL)
##- {
##-   if ((!is.atomic(x))||!is.character(x)) {
##-     warning(":i.getopt: argument 'x' must be of mode character")
##-     return(NULL)
##-   }
##-   if (length(x)>1) {
##-     warning(":i.getopt: argument 'x' has length >1. Only 'x[1]' is used.")
##-     x <- x[1]
##-   }
##-   if (is.null(default))
##-     return(i.getoption(x))
##-   if (x %in% names(.i.getopt)&&!is.null(i.getoption(x)))
##-         i.getoption(x)
##-     else default
##- }
## -----------------------------------------------------------------------
regroptions <-
  function (x=NULL, list=NULL, default=NULL, regroptions = NULL,
            assign=TRUE, ...)
{ ##
  lrgdef <- get("regroptionsDefault", pos=1)
  lnewo <- loldo <-
    if (is.null(regroptions)) {
    if (exists(".regroptions", where=1)) get(".regroptions", pos=1)
    else  lrgdef
    } else regroptions
  largs <- c(list, list(...))
  lop <- c(names(largs))
  ##
  if (!is.null(x)) {  ## asking for options
    if(!is.character(x)) {
      warning(":regroptions: First argument 'x' must be of mode character")
      return(NULL)
    }
    return(if(length(x)==1) loldo[[x]] else loldo[x])
  }
  ## ---
  if (length(default) && u.notfalse(default)) { ## get default values
    if (u.true(default)) default <- "all"
    if (!is.character(default))
      stop("!regroptions! Argument 'default' must be of mode character or logical")
    if (default[1]=="all") return(regroptions(list=regroptionsDefault, assign=assign))
    ## resets all available components
    if (default[1]=="unset")
      return(regroptions(list=lrgdef[names(regroptionsDefault)%nin%names(loldo)],
                         assign=assign) )
    if (any(default!="")) {
      llopt <- regroptionsDefault[default[default%in%names(regroptionsDefault)]]
      return( regroptions(list=llopt) )
    }
  }
  ## set options
  ## check
  largs <- check.option(list=largs) ## fname="regroptions"
  if (length(largs))  lnewo[names(largs)] <- largs
  lo <- intersect(names(largs),names(loldo))
  if (length(lo)) attr(lnewo, "old") <- loldo[lo]
  if (assign) assign(".regroptions", lnewo, pos=1)
  ## assignInMyNamespace does not work
  invisible(lnewo)
}
##- ## -----------------------------------------------------
##- c.colors.ra <- c("gray3","gray2","blue","cyan","darkgreen","green",
##-                  "burlywood4","burlywood3","burlywood4")
##- c.colors <- c("black","red3","deepskyblue3","darkgreen",
##-               "orange2","purple2","blue2","green3","brown",
##-               "violet", "aquamarine3", "brown3") ## "darkgoldenrod3", "burlywood", 
.regroptions <- regroptionsDefault <- list(
  digits = 4, 
  ##-   colors = c.colors, colors.ra = c.colors.ra,
  ##-   mar=c(3,3,3,1), mgp=c(2,0.8,0), plext=0.05,
  regr.contrasts = c(unordered="contr.wsum", ordered="contr.wpoly"),
  factorNA = TRUE, testlevel = 0.05,  termtable = TRUE, vif = TRUE,
  show.termeffects = TRUE, termcolumns = c("coef",  "df", "ciLow","ciUp","R2.x",
    "signif", "p.value", "p.symbol"),
  termeffcolumns = "coefsymb",
  coefcolumns = c("coef", "ciLow","ciUp", "signif", "p.value", "p.symbol"),
  na.print = ".",
  notices = TRUE, 
  ## smoothFunction = "smoothRegr", smoothMinobs = 8,
  debug = 0
  )
## -----------------------------------------------------------------------
check.option <- function(optname, value, list = NULL) {
  if (is.null(list)) list <- setNames(list(value), optname)
  lnl <- length(list)
  loptnames <- names(list)
  for (lil in seq_len(lnl)) {
    lnm <- loptnames[lil]
    lvalue <- list[[lnm]]
    lcheck <- regroptionsCheck[[lnm]]
    if (length(lcheck)) {
##    if (is.list(lcheck)) {
      if (!is.list(lcheck[[1]])) lcheck <- list(lcheck)
      lnopt <- length(lcheck)
      lmsg <- rep("", lnopt)
      for (lj in seq_len(lnopt)) {
        lch <- lcheck[[lj]]
        lfn <- get(lch[[1]])
        lmsg[lj] <- lmsgj <-
          switch(paste("v",length(lch),sep=""), v0="", v1=lfn(lvalue),
                 v2=lfn(lvalue, lch[[2]]), v3=lfn(lvalue, lch[[2]], lch[[3]]),
                 "")
        if (lmsgj=="") break
      }
##    }
      if (all(lmsg!="")) {
        warning(":check.option: argument '", lnm,
                "' not suitable. It should\n    ",
                paste(lmsg, collapse=" -- or \n  "),
                "\n  instead of (str())\n    ", format(str(lvalue)))
        list[lnm] <- NULL
      }
    }
  }
  list
}
## -----------------------------------------------------------------------
regroptionsCheck <- list()
## ---------------------------------------------------------------------------------
i.getoption <- function(opt, plo=NULL) {
  ## opt is character, plo list or NULL
  lpldef <- get("regroptionsDefault", pos=1)
  if (is.null(plo))
    plo <- get("regroptions", envir=parent.frame()) ## list in calling fn
  if (is.function(plo)) plo <- NULL
  lopt <- plo[[opt]]  
  if (is.null(lopt)||(!is.function(lopt)&&all(is.na(lopt)))) ## NULL or NA
      lopt <- regroptions(opt)
  else {lopt <- check.option(opt, lopt)
    if (length(lopt)) lopt <- lopt[[1]]
  }
  if (length(lopt)==0) lopt <- lpldef[[opt]]
##  names(lopt) <- opt
  lopt
}
## -----------------------------------------------------------
i.getopt <- function(opt, plo = NULL) {
  ldef <- get("regroptionsDefault", pos=1)
  if (is.null(plo))
    plo <- get("regroptions", envir=parent.frame()) ## list in calling fn
  lnam <- as.character(substitute(opt))
  lopt <- opt 
  if (is.function(plo)) plo <- NULL
  if (is.null(lopt)||(is.atomic(lopt)&&all(is.na(lopt))))
    lopt <- plo[[lnam]]
  if (is.null(lopt)||(is.atomic(lopt)&&all(is.na(lopt)))) 
    lopt <- regroptions(lnam)
  else unlist(check.option(lnam, lopt))   ## check
  if (is.null(lopt)) lopt <- ldef[[lnam]]
##  names(lopt) <- opt
  lopt
}
## ==========================================================================
shift <- function(x, k = 1)
  if (k>0) c(rep(NA, k), last(x, -k)) else
  if (k==0) x else c(last(x, k), rep(NA, -k))
## --------------------------------------------------------------------------
notice <- function(..., printnotices = NULL)
  if (i.getopt(printnotices)) message("Notice in ",...)


##- plotregr <- function(x) {}
##- .onLoad <- function() {
##-   plotregr <- plgraphics::plregr
##- }
##-   require("knitr")
##-   tools::vignetteEngine("knitr")
##- ## , weave = vweave, tangle = vtangle, pattern = "[.]Rmd$", package = "knitr")

