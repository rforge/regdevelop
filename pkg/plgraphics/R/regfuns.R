## ===================================================================
residuals.regr <- function(object, type=NULL, standardized=FALSE, ...)
{
##-   if (!is.na(pmatch(type,"condquant"))) {
##-     ## this seems to apply only if residuals.regr is called explicitly
  lcall <- match.call()
  lcall$type <- if (is.null(type)||is.na(type)) NULL else type
  lff <- object$fitfun
  if (lff=="glm" && !is.null(type) && substr(type,1,4)=="cond") lff <- "polr"
  lcall[[1]] <-
    switch(as.character(lff),
           "polr" = quote(residuals.polr),
           "survreg" = quote(residuals.regrsurv),
           "coxph" = quote(residuals.regrcoxph),
           quote(residuals))
  class(object) <- setdiff(class(object), "regr")
  lcall$object <- object
  rr <- structure( eval(lcall, envir=parent.frame()), type=type)
  if (standardized) {
    lstr <- i.stresx(object, resid=rr)
    attr(rr,"stresiduals") <- lstr$stresiduals
    attr(rr,"leverage") <- lstr$leverage
    attr(rr,"strratio") <- lstr$strratio
  }
  rr
}
## ===================================================================
residuals.regrcoxph <- function(object, type=NULL, ...)
{
  ## Purpose:    conditional quantiles and random numbers for censored obs
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: Aug 2010
  if(u.nuna(type)) type <- "CoxSnellMod"
  if (type!="CoxSnellMod") {
    object$residuals <- object$resid.orig
    return(structure( survival:::residuals.coxph(object, type=type, ...),
                     type=type) )
  }
  lnaaction <- object$na.action
  lres <- object$residuals
  if (inherits(lres, "condquant")) return(lres)
  ly <- object$y
  lst <- ly[,2]  # status
  ltt <- attr(object$response, "type")
  ltl <- length(ltt)>0&&ltt=="left"
  ##  martingale --> coxsnell
  lres <- lst - lres
  lrs <- qnorm(exp(-lres)) # --> normalscores
  ## censoring
  lii <- lst!=1
  li <- which(lii)
  if (length(li)) {
      llim <- if(ltl) cbind(-Inf,lres[li]) else cbind(lres[li],Inf)
      lr <- condquant(llim, "normal", sig=1) ## ???
      lres[li] <- lr[,"median"]
      lr[,'index'] <- which(naresid(lnaaction, lii))
  } else lr <- NULL
  structure( naresid(lnaaction, lres), condquant=lr, type=type)
}
## ===================================================================
residuals.regrsurv <- function(object, type=NULL, ...)
{
  ## Purpose:    conditional quantiles and random numbers for censored obs
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 24 Oct 2008, 10:16
  if(u.nuna(type)) type <- "condquant"
  if (type!="condquant") 
    return(structure( survival:::residuals.survreg(object, type, ...),
                     type=type) )
  lnaaction <- object$na.action
  ly <- object$y  ## log for weibull
  lsig <- object$sigma
  if (length(lsig)==0) lsig <- summary(object)$scale
  if (length(lsig)==0) {
    warning(":residuals.regrsurv: no sigma found. Setting =1")
    lsig <- 1
  }
  lfit <- object$linear.predictors
##  lres <- ly[,1]-lfit
  lres <- ly[,1]-lfit
  ldist <- if (length(object$dist)>0) object$dist  else  "normal"
  li <- match(ldist, c("weibull","lognormal","loglogistic"))
  if (!is.na(li)) ldist <- c("revgumbel","normal","logistic")[li]
  ## for user-defined survreg.distributions with transformation,
  ##   this is not enough.
  ## censoring
  lst <- ly[,2]  # status
  ltt <- attr(object$response, "type")
  ltl <- length(ltt)>0&&ltt=="left"
  lii <- lst!=1
  li <- which(lii)
  if (length(li)) {
      llim <- if(ltl) cbind(-Inf,lres[li]) else cbind(lres[li],Inf)
      lr <- condquant(llim, ldist, lsig)
      lres[li] <- lr[,"median"]
      lr[,'index'] <- which(naresid(lnaaction, lii))
  } else lr <- NULL
  structure( naresid(lnaaction, lres), condquant=lr, type=type)
}
## ==============================================================
nobs.survreg <- function(object, use.fallback = TRUE) {
  lnobs <- length(object$linear.predictors)
  if (lnobs==0) lnobs <- NROW(residuals(object))
  lnobs
}
nobs.coxph <- function(object, use.fallback = TRUE) {
  object$n
}
## ===================================================================
residuals.polr <- function(object, ...) ## na.action=object, 
{
  ## Purpose:   residuals for cumulative logit regression
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  2 Oct 2007, 11:31
  if (!((lpolr <- inherits(object, "polr"))|
        (lbin <- inherits(object, "glm")&&object$family=="binomial")))
    stop ("!residuals.polr! unsuitable first argument")
  lnaaction <- object$na.action
  lmodel <- object$model
  if (is.null(lmodel)) lmodel <- model.frame(object)
  ly <- model.response(lmodel)
  if (length(ly)==0) stop ("!residuals.polr! bug: no response values found")
  if (length(dim(ly))) {
    warning(":residuals.polr: returning simple deviance residuals for non-binary (grouped) data")
    return(residuals(object, type="deviance"))
  }
  ly <- as.numeric(ordered(ly))
  ##  ly <- naresid(object$na.action, ly)
  object$na.action <- NULL
  lfit <- fitted(object, type="link")
  ##  if (length(lnaaction)) lfit <- lfit[-lnaaction]
  lthres <- c(-100, if (lpolr) object$zeta else 0, 100)
  llim <- structure(cbind(lthres[ly],lthres[ly+1])-lfit,
                    dimnames=list(names(lfit), c("low","high")))
  lr <- cbind(condquant(llim,"logis"),fit=lfit,y=ly)
  lr[,"index"] <- which(naresid(lnaaction, rep(TRUE,nrow(lr))))
  ##
  structure(naresid(lnaaction, lr[,"median"]), condquant=lr) 
}
## ===========================================================================
linear.predictors <- function(object) {
  llp <- object$linear.predictors
  if (is.null(llp)) llp <- object$fitted.values
  if (is.null(llp))
    stop("linear.predictors! no component linear predictor")
  naresid(object$na.action, llp)
}
linpred <- linear.predictors
## ===========================================================================
fitcomp <-
  function(object, data=NULL, vars=NULL, transformed=FALSE, 
           se=FALSE,  xm=NULL, xfromdata=FALSE, noexpand=NULL, nxcomp=51)
{
  ## Purpose:    components of a fit
  ## ----------------------------------------------------------------------
  ## !!! make nxcomp >= maximal number of factor levels !!!
  ## !!! why vars??? can possibly be much simpler!!!
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  6 Jan 2004, 22:24
  ltransf <- i.def(transformed, FALSE, TRUE, FALSE)
  if (ltransf) data <- model.frame(object)[,-1,drop=FALSE]
  else {
    if (length(data)==0) {
      data <- object$allvars
      if (length(data)==0)
        data <- eval(object$call$data)
    }
  }
  if (se & inherits(object, c("polr"))) {
    warning(":fitcomp: standard errors for fit not available for ",
            paste(class(object), collapse="  "), "models. 'se' set to FALSE")
    se <- FALSE
  }
  lform <- formula(object)[-2]
  lvars <- if (ltransf) names(data) else all.vars(lform)
  if (lIlns <- inherits(object, "nls")) ## c(eval(object$call$nonlinear),FALSE)[1]
    lvars <- lvars[match(lvars,names(coefficients(object)),nomatch=0)==0]
  lvmiss <- setdiff(lvars,names(data))
  if (length(lvmiss)) 
    stop(paste("!fitcomp! variable(s)", paste(lvmiss,collapse=", "),
               "not in data"))
  if (length(lvars)==0)
    stop("!fitcomp! no variables selected")
  ldata <- data[,lvars,drop=FALSE]
  for (lj in 1:length(ldata))
    if (is.character(ldata[[lj]])|is.factor(ldata[[lj]]))
      ldata[,lj] <- factor(ldata[,lj])
  lformfac <- NULL
  if (!lIlns) {
    if (length(c(grep("factor *\\(", format(lform)),
        grep("ordered *\\(", format(lform))))) {
      warning(
        ":fitcomp: Using 'factor(...)' or 'ordered(...)' in the formula ",
        "is hazardous for fitcomp.\n",
        "  : I try to continue.")
      lformfac <- i.findformfac(lform)
##          return(structure(list(comp=NULL), class="try-error") )
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
      lnhalf <- ceiling(sum(!is.na(lv))/2)
      xm[1,lj] <-
        if (is.factor(lv)) {
          levels(lv)[
          if (is.ordered(lv)) sort(as.numeric(lv))[lnhalf] else
                     which.max(table(as.numeric(lv)))
                   ]
        } else   ## median(as.numeric(lv),na.rm=TRUE)
        sort(lv)[lnhalf]
        ## median should be attained in >=1 cases
    }
  }
  if (is.null(attr(terms(object), "predvars"))) { # from  model.frame
    lmf <- if (ltransf) data
           else lm(formula(object), data=data, method="model.frame")
    lterms <- attr(lmf,"terms")
    attr(object$terms,"predvars") <- attr(lterms,"predvars")
  }
  ltype <- if (inherits(object,"coxph")) "lp"  else "response"
  if (inherits(object, "glm")) ltype <- "link"
  if (ltransf) {
    lobj <- structure(object, class="list")
    lobj[["terms"]] <- lterms <- delete.response(terms(object))
    lobj[["x"]] <- model.matrix(lobj, structure(xm, terms=lterms))
    lobj[["offset"]] <- NULL
    ## somehow, it needs the  terms  in both places...
    lprm <- c(predict(structure(lobj, class=class(object)), type=ltype))
  } else 
    lprm <- c(predict(object, newdata=xm, type=ltype)) # lf.
  lny <- length(lprm)
##  expand to matrix
  if (xfromdata) {
    lx <- ldata
  } else {
    lnxj <- sapply(ldata,
                   function(x) if (is.factor(x)) length(levels(x)) else 0)
    if(!is.null(noexpand) && is.numeric(noexpand))
      noexpand <- names(noexpand)[noexpand>0]
    noexpand <- c(noexpand, lformfac) ## 
    lvconv <- names(ldata) %in% noexpand
    names(lvconv) <- names(ldata)
    if (any(lvconv)) lnxj[lvconv] <-
      sapply(ldata[lvconv], function(x) length(unique(x)) )
    lnxc <- max(nxcomp, lnxj)
    lx <- ldata[1,,drop=FALSE][1:lnxc,,drop=FALSE]
    row.names(lx) <- 1:lnxc
  }
  ##  lxm: data.frame of suitable dimension filled with "median"
  lxm <- lx
  for (lv in names(lxm)) lxm[,lv] <- xm[,lv]
  ##
  lvcomp <- names(ldata)
  if (!is.null(vars)) lvcomp <- intersect(lvcomp, vars)
  if (is.null(lvcomp)) {
    warning(":fitcomp: no variables found. Result is NULL")
    return(NULL)
  }
##  components
  lcomp <- array(dim=c(nrow(lx), length(lvcomp), lny)) 
  dimnames(lcomp) <- list(dimnames(lx)[[1]], lvcomp, names(lprm))
  lcse <- if (se) lcomp  else NULL
  if (ltransf)   {
    lobj[["residuals"]] <- rep(NA,nrow(lx))
    ## needed since predict reads tne number of obs from $residuals
    lobj[["weights"]] <- rep(1,nrow(lx))
    ## forces  predict.lm  to invert X from scratch
    lsigma <- 
      if (is.null(scale <- object$scale)) {
        w <- object$weights
        r <- object$residuals
        rss <- sum(if (is.null(w)) r^2 else r^2 * w)
        df <- object$df.residual
        sqrt(rss/df)
      } else scale
  }
  ## ------------------------
  for (lv in lvcomp) {
    if (xfromdata) {
      ld <- lxm
      ld[,lv] <- ldata[,lv]
      lfc <- sapply(ld,is.factor) # eliminate extra levels of factors
      if (any(lfc)) ld[lfc] <- lapply(ld[lfc], factor)
    } else { # +++
      ldv <- ldata[,lv]
      if (lnxj[lv]) { # factor levels
        ldx <- if (lvconv[lv]) sort(unique(ldv)) else factor(levels(ldv))
        lnl <- length(ldx)
        ld <- lxm[1:lnl,,drop=FALSE]
        ld[,lv] <- ldx
##-         lx[,lv] <- factor(c(1:lnl,rep(NA,lnxc-lnl)),labels=levels(ldv))
        lx[,lv] <- c(ldx,rep(NA,lnxc-lnl))
        ##
        if (ltransf) {
          lobj[["x"]] <- model.matrix(lobj, structure(ld, terms=lterms))
          lpr <- predict(structure(lobj, class=class(object)), scale=lsigma, 
                     type=ltype, se.fit=se)
        } else 
          lpr <- try( predict(object, newdata=ld, se.fit = se),
                     silent=TRUE)
        if (class(lpr)=="try-error") {
          warning(":fitcomp: no fitcomp for variable  ", lv)
          ## predict finds new levels of formfac variables
          next
        }
        if (se) {
          lc <- lpr$fit
          lcse[1:lnl,lv,] <- lpr$se.fit
        } else lc <- lpr
        lcomp[1:lnl,lv,] <- lc
        next # end for loop
      } else { # continuous var
        ld <- lxm
        lx[,lv] <- ld[,lv] <-
          seq(min(ldv,na.rm=TRUE),max(ldv,na.rm=TRUE),length=lnxc)
      } # ---
    } # +++
    ## continuous variable or xfromdata
    if (ltransf) {
      lobj[["x"]] <- model.matrix(lobj, structure(ld, terms=lterms))
      lpr <- predict(structure(lobj, class=class(object)), scale=lsigma, 
                     type=ltype, se.fit=se)
    } else 
      lpr <- predict(object, newdata=ld, type=ltype, se = se) # lf.
    if (se) {
      lcomp[,lv,] <- lpr$fit
      lcse[,lv,] <- lpr$se.fit
    } else lcomp[,lv,] <- lpr
  }
  if (lny==1) {
    dim(lcomp) <- dim(lcomp)[1:2]
    dimnames(lcomp) <- dimnames(lx[,lvcomp,drop=FALSE])
    lcomp <- lcomp-lprm
    if (se) {
      dim(lcse) <- dim(lcse)[1:2]
      dimnames(lcse) <- dimnames(lcomp)
    }
  } else   lcomp <- sweep(lcomp,3,lprm)
  list(comp=lcomp, x=lx[,lvars,drop=FALSE], xm=xm[,lvars,drop=FALSE], se=lcse)
}
## ==========================================================================
i.findformfac <- function(formula) {
  ## find variable involved in explicit factor terms in formula
  lfo <- format(formula)
  lmf <- c(gregexpr("(factor *\\([^)]*\\))", lfo),
           gregexpr("(ordered *\\([^)]*\\))", lfo) )
  lf <- function(x) 
    if(x[1]!=-1) substring(lfo, x, x+attr(x,"match.length"))
  all.vars(as.formula(
        paste("~",paste(unlist(lapply(lmf, lf)), collapse="+"))))
}
## =======================================================================
leverage <- function(object)
{
  ## Purpose:  extract leverages
  ## Author: Werner Stahel, Date: 10 Nov 2010, 09:31
  lh <- object$leverage
  if(is.null(lh)) lh <- hatvalues(object)
  names(lh) <- names(object$resid)
  naresid(object$na.action, lh)
}
## ==========================================================================
i.stres <-
  function(x, residuals=x$residuals, sigma=x$sigma, weights=x$weights,
           leveragelim = c(0.99, 0.5))
{ ## sigma=x$sigma, 
  ## Purpose:  calculate  hat  and  standardized residuals
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  1 Mar 2018, 15:45
  leveragelim <- i.def(leveragelim, c(0.99, 0.5), valuefalse=c(0.99, 0.5))
  sigma <- i.def(sigma, x$sigma)
  if (is.null(sigma)||(!is.finite(sigma))||sigma<=0) {
    warning(":i.stres: sigma is missing or <=0. I use sigma=1")
    sigma <- 1
  }
  if (inherits(x, "nls")) {
    warning(":i.stres: no leverage and standardized residuals",
            " are available for a nonlinear model")
    return(list(leverage = NULL, stresiduals = NULL, strratio = NULL))
  }
  llev <- x$leverage
  if (length(llev)==0) llev <- hatvalues(x)
##-   {
##-     lmm <- x[["x"]]
##-     if (length(lmm)==0) lmm <- model.matrix(terms(x),x$model)
##-     llev <- hat(lmm)
  ##-   }
  llevlim <- leveragelim[1]
  if (mean(llev)>llevlim) {
    warning(":i.stres: leveragelim not suitable. I use 0.99")
    llevlim <- 0.99
  }
  lres <- residuals
  if (length(lres)==0) {
    warning(":regr/i.stres: no residuals round -> no standardized res.")
    return( list(leverage=llev) )
  }
  if (length(llev)!=NROW(lres)) {
    warning(":regr/i.stres: no leverages available, I set them 0")
    llev <- rep(0,NROW(lres))
  }
  names(llev) <- rownames(as.matrix(lres))
  ##
  lstrratio <- 1 / outer(sqrt(1-pmin(llevlim,llev)), sigma)
  lnwgt <- length(weights)
  lIwgt <- lnwgt==NROW(lres)
  if (lnwgt) {
    if (lIwgt) lstrratio <- lstrratio / sqrt(weights)
    else warning(":regr/i.stres: weights not suitable -> not used")
  }
  lstres <- as.matrix(lres)*lstrratio
  list(leverage = llev, stresiduals = structure(lstres, weighted=lIwgt),
       strratio = lstrratio)
}
## -----------------------------------------------------------------------
i.stresx <-
  function(x, resid=NULL, weights=NULL, leveragelim = c(0.99, 0.5))
{ ## Purpose:  calculate  hat  and  standardized residuals
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel
  sigma <- i.def(i.def(x$sigma, x$scale), 1)
  llev <- leverage(x)
  if (length(llev)==0) {
    lmm <- x[["x"]]
    if (length(lmm)==0) lmm <- model.matrix(terms(x),x$model)
    llev <- naresid(x$na.action, hat(lmm))
  }
  resid <- as.data.frame(i.def(resid, residuals(x)))
  lstrratio <- i.def(x$strratio,
                     outer(1/sqrt(1-pmin(leveragelim[1],llev)), sigma, "/") )
  lwgt <- weights
  if (lIwgt <- length(lwgt)==NROW(resid)) lstrratio <- lstrratio * sqrt(lwgt)
  lstres <- resid*lstrratio
  lmres <- ncol(lstres)
  for (lj in 1:lmres) {
    if (length(lcq <- attr(resid[,lj], "condquant")))
      attr(lstres[,lj], "condquant") <-
        cbind(lcq[,1:4]*lstrratio[lcq[,"index"]], lcq[,-(1:4)])
  }
  if (lmres==1) lstres <- lstres[[1]]
  list(leverage = llev, stresiduals = structure(lstres, weighted=lIwgt),
       strratio = lstrratio)
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
  if (length(x)==0) stop("!condquant! bug: no data")
  ## functions for calculating probab. and quantiles
  fp <- switch(dist, normal=pnorm, gaussian=pnorm, unif=function (x) x,
               logis=plogis, logistic=plogis, revgumbel=prevgumbel)
  if (is.null(fp)) stop(paste("!condquant! distribution ", dist, " not known"))
  fq <- switch(dist, normal=qnorm, gaussian=qnorm, unif=function (x) x,
               logis=qlogis, logistic=qlogis, revgumbel=qrevgumbel)
  ##
  x <- rbind(x)  ## rbind for a single observation -> vector x
  lii <- which(x[,1]!=x[,2])
  rr <- structure(matrix(lii,length(lii),6),
                  dimnames=
                    list(row.names(x)[lii],
                         c("median","lowq","uppq","random","prob","index")),
                  class="condquant")
  if (length(lii)==0) {
    message(".condquant. All intervals of length 0. ",
            "I return a matrix with 0 lines")
    return(rr)
  }
  lx <- t(apply(x[lii,], 1,sort))
  lp <- fp(lx/sig)
  lpp <- rbind(rbind(lp)%*%rbind(c(0.5,0.75,0.25),c(0.5,0.25,0.75)))
  lprand <- lp[,1]+(lp[,2]-lp[,1])*
    runif(nrow(lp),(1-randomrange)/2,(1+randomrange)/2)
  rr[,1:4] <- cbind(median=fq(lpp[,1]),lowq=fq(lpp[,2]),
                    uppq=fq(lpp[,3]), random=fq(lprand))*sig
  rr[,5] <- lp[,2]-lp[,1]
  if (any(lp0 <- lp[,2]<=0)) rr[lp0,1:4] <- matrix(lx[lp0,2],sum(lp0),4)
  if (any(lp1 <- lp[,1]>=1)) rr[lp1,1:4] <- matrix(lx[lp1,1],sum(lp1),4)
  if (any(c(lp0,lp1)))
    warning(":condquant: probabilities <0 or >1")
  ##
  rr
}
## ==========================================================================
fitted.regr <-
  function (object, type=NULL, ...)
{
  if (is.null(type)&&pmatch("fitted",names(object),nomatch=0))
    return( naresid(object$na.action, object$fitted) )
  lres <- object$residuals
  if (inherits(lres, "condquant"))
    structure(lres[,"fit"], names=row.names(lres))  else  {
      class(object) <- setdiff(class(object), "regr")
      predict(object, type=type, ...)
    }
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
    lfit <- object$fitted
  if (type=="class")
    lfit <- factor(max.col(lfit), levels = seq_along(object$lev),
            labels = object$lev)
  naresid(object$na.action,lfit)
}
## ==========================================================================
predict.regr <-
  function (object, newdata = NULL, type = NULL, se.fit=FALSE,
            scale = object$sigma, df=object$df.residual, ...)
{
  if (length(type)==0)
    type <- if (inherits(object,c("glm","polr"))) "link" else "response"
  if (se.fit & inherits(object, c("polr"))) {
    warning(":predict.regr: standard errors for fit not available for ",
            paste(class(object), collapse="  "), "models. Set to FALSE")
    se.fit <- FALSE
  }
  ## !!!
  if (object$fitfun=="rlm")
    if (!is.matrix(object[["x"]]))
      object$x <- model.matrix(formula(object), data=object$allvars) ## !!! was $model
  class(object) <- class(object)[-1]
  ldt <- newdata
  if (is.null(ldt)) return(
    predict(object, type=type, scale=object$sigma,
            dispersion=object$dispersion^2, ... ) )
  ## analyze variables
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
          df=df, dispersion=object$dispersion^2, se.fit=se.fit, ... )
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
    attr(Terms, "intercept") <- 1
    environment(Terms) <- environment() ## ! WSt
    m <- model.frame(Terms, newdata, na.action = function(x) x,
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts = object$contrasts)
    ## changed by WSt
    coef <- coefficients(object)
    eta <- drop(X[,-1] %*% coef) ## without the intercept
    n <- nrow(X)
    q <- length(object$zeta)
    ##  pgumbel <- function(q) exp(pweibull(log(q))) # ???
    pfun <- switch(object$method, logistic = plogis, probit = pnorm,
                   cloglog = prevgumbel, cauchit = pcauchy)
    cumpr <- matrix(pfun(matrix(object$zeta, n, q, byrow = TRUE) - eta), , q)
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
## ==========================================================================
predict.mlm <-
  function (object, newdata=NULL, scale = NULL, df = Inf,
            interval = c("none", "confidence", "prediction"),
            se.fit = FALSE, level = 0.95, type = c("response", "terms"),
            terms = NULL, na.action = na.pass,
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
simresiduals <- function(object, ...)  UseMethod("simresiduals")
  ## Purpose:   simulate residuals according to regression model
  ##            by permuting residuals of actual model or by random numbers
  ## ----------------------------------------------------------------------
  ## Arguments:  simfunction: how are residuals generated?
## ---------------------------------------------------------------------
simresiduals.glm <- function(object, nrep=19, simfunction=NULL,
                             glm.restype="working", ...)
{
  lcall <- object$call
  if ("weights"%in%names(lcall)) {
    warning(":simresiduals: I cannot simulate for weighted regression (yet)")
    ## get_all_vars  contains  weights  without parentheses -> danger!
    return(NULL)
  }
  loverd <- attr(object$scale, "fixed")
  if (length(loverd) && !loverd)
    warning(":simresiduals: Cannot simulate from overdispersed model.",
            " Using dispersion 1")
  lcall[[1]] <- as.name("glm")
  ## -------
  lnaaction <- object$na.action
  ldata <- object$allvars
  if (is.null(ldata)) {
      ldata <- if (u.debug())  eval(lcall$data)  else
      try(eval(lcall$data))
      if (class(ldata)=="try-error"||is.null(dim(ldata))) {
        warning(":simresiduals: data not found -> No simulated residuals")
        return(NULL)
      }
  }
  if (length(lnaaction)) ldata <- ldata[-lnaaction,]
  lfit <- object$fitted.values
  ## prepare call
  lform <- update(formula(object), .Y. ~.)
  lynm <- all.vars(lform[[2]])
  environment(lform) <- environment()
  lcl <- call("glm", formula=lform, data=as.name("ldata"),
              family=object$family, start=object$coef, model=FALSE,
              y=FALSE, na.action=lcall$na.action)
  lfam <- object$distrname
  if (length(lfam)==0) lfam <- object$family$family
  if (is.null(lfam)) lfam <- ""
  ly <- object$response
  ln <- nrow(ldata)
  lone <- rep(1,ln)
  if (!is.function(simfunction)) {
    if(lfam%in%c("binomial","quasibinomial")) {
      if (NCOL(ly)==1 && length(unique(ly))!=2) {
        warning(":simresiduals: binomial distribution with ",
                "unsuitable response.\n  No residuals simulated")
        return(list(simres=numeric(0)))
      }
      simfunction <-
        if(NCOL(ly)==1) function(n, fit, sig=NULL) rbinom(n, lone, fit)
      else {
        lnbin <- ly[,1]+ly[,2]
        function(n, fit, sig=NULL) {
          ly1 <- rbinom(n, lnbin, fit)
          cbind(N1=ly1,N2=lnbin-ly1)
        }
      }
    } else {
      if (lfam%in%c("poisson","quasipoisson")) 
        simfunction <- function(n, fit, sig) rpois(ln, fit)
      else {
        warning(":simresiduals: not (yet) available for this ",
                "type of model.\n  No residuals simulated")
        return(list(simres=numeric(0)))
      }
    }
  }
  ## ---
  lsimres <- matrix(NA, ln, nrep)
  for (lr in 1:nrep) {
    ldata$.Y. <- simfunction(ln, lfit)
    lrs <- eval(lcl, environment())
    lsimres[,lr] <-
      if (substr(glm.restype,1,4)=="cond")
          residuals.polr(lrs)[,"random"] else
          residuals(lrs, type=glm.restype)
  }
  naresid(lnaaction, lsimres)
}
## ======================================================================
simresiduals.default <-
  function(object, nrep=19, simfunction=NULL, stresiduals=NULL,
           sigma=object$sigma, ...)
  ## glm.restype="deviance")
{
  ## Purpose:   simulate residuals according to regression model
  ##            by permuting residuals of actual model or by random numbers
  ## ----------------------------------------------------------------------
  ## Arguments:  simfunction: how are residuals generated?
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 10 Aug 2008, 07:55
##-   if (!class(object)[1]%in%c("regr","lm","glm")) {
##-     warning(":simresiduals: ",
##-             "I can simulate only for `regr`, `lm`, or `glm` objects")
##-     return(NULL)
  ##-   }
  if (!inherits(object, c("lm", "lmrob")))
      stop("!simresiduals! I cannot simulate for model class  ",
           paste(class(object), collapse="  "))
  lcall <- object$call
  if ("weights"%in%names(lcall)) {
    warning(":simresiduals: I cannot simulate for weighted regression (yet)")
    ## get_all_vars  contains  weights  without parentheses -> danger!
    return(NULL)
  }
  ldata <- object$allvars
  if (is.null(ldata)) {
      ldata <- if (u.debug())  eval(lcall$data)  else
      try(eval(lcall$data))
      if (class(ldata)=="try-error"||is.null(dim(ldata))) {
        warning(":simresiduals: data not found -> No simulated residuals")
        return(NULL)
      }
  }
  if ("omit"==class(lnaaction <- object$na.action)) {
    lkeep <- rep(TRUE, nrow(ldata))
    lkeep[lnaaction] <- FALSE
    ldata <- ldata[lkeep,]
  }
##  lcall$na.action <- NULL
##-   lfit <- object$fitted.values
##-   if (is.null(lfit)) lfit <- object$linear.predictors
##-   if (is.null(lfit)) lfit <- 0
  lres <- stresiduals
  if (is.null(lres)) lres <- attr(object$stresiduals, "stresiduals")
  lnc <- NCOL(object$coefficients)
  ## -------
  if (lrgen <- length(simfunction)>0) {
    if (!is.function(simfunction)) simfunction <- rnorm
    ## ---
  ## weibull not yet implemented
    lsig <- sigma
    if (length(lsig)==0) lsig <- rep(1,lnc)  ## only standardized res useful!
    if (length(lsig)!=lnc) {
      warning(":simresiduals.default: inadequate length of 'sigma'")
      return(NULL)
    }
    if (length(lres)==0) lres <- matrix(0, nrow(ldata),length(lsig)) 
  } else { ## not yet for multivariate !!!
    lrs <- if (length(lres)==0) matrix(0,1,1) else as.matrix(lres)
    if(all(lrs[,1]==lrs[1,1])||all(!is.finite(lrs[,1]))) {
      lres <- object$resid
      if(is.null(lres)) {
        warning(":simresiduals: no (distinct) residuals found",
                "-> No simulated residuals")
        return(NULL)
      } else lres <- as.matrix(lres)
    }
    else
      lres <- sweep(as.matrix(lres), 2, object$sigma, "*") ## may still be d.f
    if (length(lcq <- attr(lres[,1], "condquant"))>0) { ##!!! mult!
      li <- lcq[,"index"]
      lres[li,1] <- lcq[,"random"] ##structure(lres[,"random"], names=row.names(lres))
    }
  }
  lresj <- lres[,1] ## wrong for multivariate
  if (nrow(ldata)!=length(lresj)) {
    li <- match(names(lresj),row.names(ldata))
    if (anyNA(li)) {
      warning(":simresiduals: data not suitable -> No simulated residuals")
      return(NULL)
    }
    ldata <- ldata[li,]
  }
    ##!!! weights
  lina <- is.na(lresj)
  if (any(lina)) {
    lresj <- lresj[!lina]
    ldata <- ldata[!lina,]
  ##  lfit <- rep(lfit,length=length(lresj))[!lina]
  }
  if (nrow(ldata)<=2) {
    warning(":simresiduals: <=2 residuals found -> No simulated residuals")
    return(NULL)
  }
## ---
  ## prepare call
  lcall$data <- as.name("ldata")
  lform <- formula(object)
  lynm <- all.vars(lform[[2]])
  environment(lform) <- environment()
  lcall$formula <- lform
  lcall <- lcall[names(lcall)%nin%c("yy","fname","family","vif",
                                    "calcdisp","suffmean","termtable")]
  lcall$model <- NULL
  lcall$termtable <- NULL
  lnrow <- nrow(ldata)
  lsimres <- matrix(NA,lnrow,nrep)
  lfam <- object$distrname
  if (length(lfam)==0) lfam <- object$family$family
##  if (lfam%in%c("gaussian","Gaussian")) 
  lcall$formula <- update(lform, paste(lynm,"~.")) ## needed for transformed y
  for (lr in 1:nrep) {
    ldata[,lynm] <-
      ##-    if (lrgen) simfunction(lnrow,lfit,lsig) else lfit + sample(lresj)
      if (lrgen) simfunction(lnrow,0,lsig) else sample(lresj)
    ## this would not work with polr or other matrix residuals
    lenv <- environment()
    lrs <- eval(lcall, envir=lenv) ## update(x, formula=lfo, data=ldata)
    lrsr <- residuals(lrs)
    if (inherits(lrsr, "condquant"))
      lrsr <- attr(lrsr, "condquant")[,"random"]
    lsimres[,lr] <- lrsr
  }
  ##naresid(lnaaction, lsimres)
  lsimres
}
## ==========================================================================

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
## =========================================================================
Tobit <- function(data, limit=0, limhigh=NULL, transform=NULL, log=FALSE, ...)
{
  ## Purpose:   create a Surv object for tobit regression
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  1 Jan 2010, 21:49
##  require(survival)  ## !?!
  ltrs <- as.character(substitute(transform))
  data <- pmax(data,limit)
  lright <- !is.null(limhigh)
  if (lright) data <- pmin(data, limhigh)
  if (log[1]) { ## model.frame  evaluates  log  in  data ! Whence [1]
      transform <- logst
      ltrs <- "logst"
  }
  if (!is.null(transform)) {
    if (is.character(transform)) transform <- get(transform)
    if (!is.function(transform))
      stop("!Tobit! argument 'transform' does not yield a function")
    ldt <- transform(c(limit,limhigh,data), ...)
    data <- ldt[-1]
    limit <- ldt[1]
    if (lright) {
      limhigh <- ldt[2]
      data <- data[-1]
    }
  }
  if (sum(data<=limit,na.rm=TRUE)==0)
    warning(":Tobit: no observation <= `limit`")
  if (lright&&sum(data>=limhigh,na.rm=TRUE)<=1)
    warning(":Tobit: no observation >= `limhigh`")
  if (lright) { 
    rr <- survival::Surv(time = data, time2=data,
                         event = (data<=limit) + (data<limhigh),
                         type="interval")
    rr[,2] <- rr[,1]
    } else
      rr <- survival::Surv(data, event = data>limit, type="left")
  structure(rr, distribution="gaussian", transform=ltrs,
            limit=c(limit,limhigh), class=c(class(rr), "matrix"))
}

## =========================================================================
