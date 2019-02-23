##- Functions in  regfuns.R
##-
##- residuals.regr
##- residuals.regrcoxph
##- residuals.regrsurvreg
##- nobs.survreg
##- nobs.coxph
##- residuals.regrpolr
##- fitted.regr
##- predict.regrpolr
##- fitted.regrpolr
##- 
##- linear.predictors
##- fitcomp
##- condquant
##- Tobit
##- 
##- i.findformfac
##- i.stres
##- 
##- simresiduals
##- simresiduals.glm
##- simresiduals.default
##- 
##- drevgumbel
##- prevgumbel
##- qrevgumbel
##- rrevgumbel

## ===================================================================
residuals.regrcoxph <- function (object, type="CoxSnellMod", ...)
{
  ## Purpose:    conditional quantiles and random numbers for censored obs
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: Aug 2010
  type <- i.def(type, "CoxSnellMod", valuefalse="")
  if (type=="") return()
  if (type%nin%c("CoxSnellMod","condquant")) {
    if (!inherits(object, "coxph")) {
      warning(":residuals.regrcoxph: 'object' does not inherit from 'coxph'")
      return()
    } ## $residuals <- object$resid.orig
    return(structure( residuals(object, type=type, ...), type=type) )
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
    lr[,"index"] <- li
    lres[li] <- lr[,"median"]
    lr[,'index'] <- which(naresid(lnaaction, lii))
  } else lr <- NULL
  structure( naresid(lnaaction, lres), condquant=lr, type=type)
}
## ===================================================================
residuals.regrsurvreg <- function (object, type="condquant", ...)
{
  ## Purpose:    conditional quantiles and random numbers for censored obs
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 24 Oct 2008, 10:16
  type <- i.def(type, "condquant", valuefalse="")
  if (type=="") return()
  if (type!="condquant")
    return( residuals(object, type=type, ...) )
  if (!inherits(object, "survreg")) {
    warning(":residuals.regrsurvreg: 'object' does not inherit from 'survreg'.",
            " I use usual residuals")
    return(residuals(object, ...) )
  }
  lnaaction <- object$na.action
  ldist <- if (length(object$dist)>0) object$dist  else  "normal"
  ly <- naresid(lnaaction, object$y) 
  if (ldist %in% c("weibull", "exponential", "lognormal", "loglogistic"))
    ly[,1] <- log(ly[,1])
  lsig <- object$sigma
  if (length(lsig)==0) lsig <- summary(object)$scale
  if (length(lsig)==0) {
    warning(":residuals.regrsurvreg: no sigma found. Setting =1")
    lsig <- 1
  }
  lfit <- naresid(lnaaction, object$linear.predictors)
##  lres <- ly[,1]-lfit
  lres <- ly[,1]-lfit
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
    lr[,"index"] <- li
    lres[li] <- lr[,"median"]
    ## lr[,'index'] <- which(naresid(lnaaction, lii)) ## !! not needed (?)
  } else lr <- NULL
  structure(lres, condquant=lr, type=type, family=ldist)
}
## ==============================================================
nobs.survreg <- function (object, ...) { ##use.fallback = TRUE
  lnobs <- length(object$linear.predictors)
  if (lnobs==0) lnobs <- NROW(residuals(object))
  lnobs
}
nobs.coxph <- function (object, ...) {
  object$n
}
## ===================================================================
residuals.regrpolr <- function (object, type="condquant", ...) ## na.action=object, 
{
  ## Purpose:   residuals for cumulative logit regression
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  2 Oct 2007, 11:31
  lbin <- inherits(object, "glm") && object$family$family=="binomial" &&
    length(unique(object$y))==2
  if (!((lpolr <- inherits(object, "polr")) | lbin))
    stop ("!residuals.regrpolr! unsuitable first argument")
  type <- i.def(type, "condquant")[1]
  if (type!="condquant") return(residuals(object, type=type))
  ## ---
  if (inherits(object, "polr")) class(object) <- "regrpolr"
  lnaaction <- object$na.action
  lmodel <- object$model
  if (is.null(lmodel)) lmodel <- model.frame(object)
  ly <- model.response(lmodel)
  if (length(ly)==0) stop ("!residuals.regrpolr! bug: no response values found")
  if (length(dim(ly))) {
    warning(":residuals.regrpolr: returning simple deviance residuals for non-binary (grouped) data")
    return(residuals(object, type="deviance"))
  }
  ## ---
  ly <- as.numeric(ordered(ly))
  ##  ly <- naresid(object$na.action, ly)
  object$na.action <- NULL
  lfit <- fitted(object, type="link")
  ##  if (length(lnaaction)) lfit <- lfit[-lnaaction]
  lthres <- c(-Inf, if (lpolr) object$zeta else 0, Inf)
  llim <- structure(cbind(lthres[ly],lthres[ly+1])-lfit,
                    dimnames=list(names(lfit), c("low","high")))
  lr <- cbind(condquant(llim,"logis"),fit=lfit,y=ly) ## all obs are 'censored'
  lr[,"index"] <- which(naresid(lnaaction, rep(TRUE,nrow(lr))))
  ##
  structure(naresid(lnaaction, lr[,"median"]), condquant=lr) 
}
residuals.polr <- residuals.regrpolr
## ===========================================================================
predict.regrpolr <-
  function (object, newdata=NULL, type = c("class", "probs", "link"), ...)
  ## type link added by WSt, newdata=NULL
{
  if (!inherits(object, c("polr","regrpolr")))
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
## ===================================================================
fitted.regrpolr <- function (object, type= c("class", "probs", "link"), ...) {
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
  napredict(object$na.action,lfit)
}
## ==========================================================================
linear.predictors <- function (object) {
  llp <- object$linear.predictors
  if (is.null(llp)) llp <- object$fitted.values
  if (is.null(llp))
    stop("linear.predictors! no component linear predictor")
  naresid(object$na.action, llp)
}
linpred <- linear.predictors
## ===========================================================================
fitcomp <-
  function (object, data=NULL, vars=NULL, transformed=FALSE, 
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
##  llog <- inherits(object, c("survreg")) &&
##    object$dist %in% c("weibull","lognormal","loglogistic") ## log(predictions)
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
  ##
  ## if (inherits(object, "polr")) class(object) <- "lm"
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
  ltype <- if (inherits(object,c("survreg","coxph"))) "lp"  else "response"
  ## lp = linear predictor
  if (inherits(object, c("glm", "polr"))) ltype <- "link"
  if (inherits(object, "polr")) class(object) <- "regrpolr"
  if (ltransf) {
    lobj <- structure(object, class="list")
    lobj[["terms"]] <- lterms <- delete.response(terms(object))
    lobj[["x"]] <- model.matrix(lobj, structure(xm, terms=lterms))
    lobj[["offset"]] <- NULL
    lobj[["na.action"]] <- NULL
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
      ld <- lxm  ## suitable dimension
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
          lpr <- try( predict(object, newdata=ld, type=ltype, se.fit = se),
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
##-   if (llog) {
##-     lcomp <- log(lcomp)
##-     lprm <- log(lprm)
##-   }
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
condquant <- function (x, dist="normal", mu=0, sig=1, randomrange=0.9)
{
  ## Purpose:   conditional quantiles and random numbers
  ## works only for centered scale families
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 24 Oct 2008, 09:30
  if (length(x)==0) stop("!condquant! no data")
  x <- rbind(x)  ## rbind for a single observation -> vector x
  if (NCOL(x)!=2) stop("!condquant! 'x' must be a matrix with 2 columns")
  ## functions for calculating probab. and quantiles
  fp <- switch(dist, normal=pnorm, gaussian=pnorm, Gaussian=pnorm,
               unif=function(x) x,
               logis=plogis, logistic=plogis, revgumbel=prevgumbel)
  if (is.null(fp)) stop(paste("!condquant! distribution ", dist, " not known"))
  fq <- switch(dist, normal=qnorm, gaussian=qnorm, unif=function(x) x,
               logis=qlogis, logistic=qlogis, revgumbel=qrevgumbel)
  ##
  lii <- which(x[,1]!=x[,2])
  rr <- structure(
    matrix(lii,length(lii),8),
    dimnames=
      list(row.names(x)[lii],
           c("median","lowq","uppq","random","limlow","limup","prob","index")),
    class="condquant")
  if (length(lii)==0) {
    message(".condquant. All intervals of length 0. ",
            "I return a matrix with 0 lines")
    return(rr)
  }
  lx <- t(apply(x[lii,], 1,sort))
  lp <- fp((lx-mu)/sig)
  lpp <- rbind(rbind(lp)%*%rbind(c(0.5,0.75,0.25),c(0.5,0.25,0.75)))
  lprand <- lp[,1]+(lp[,2]-lp[,1])*
    runif(nrow(lp),(1-randomrange)/2,(1+randomrange)/2)
  rr[,1:4] <- cbind(median=fq(lpp[,1]),lowq=fq(lpp[,2]),
                    uppq=fq(lpp[,3]), random=fq(lprand))*sig
  rr[,5:6] <- lx
  rr[,7] <- lp[,2]-lp[,1]
  if (any(lp0 <- lp[,2]<=0)) rr[lp0,1:4] <- matrix(lx[lp0,2],sum(lp0),4)
  if (any(lp1 <- lp[,1]>=1)) rr[lp1,1:4] <- matrix(lx[lp1,1],sum(lp1),4)
  if (any(c(lp0,lp1))) {
    if (any(lp[,2]<0|lp[,1]>1))
      warning(":condquant: BUG: probabilities <0 or >1")
    else
      message(".condquant. probabilities <=0 or >=1")
  }
  attr(rr,"distribution") <- list(distribution = dist, mu = mu, sigma = sig)
  ##
  rr
}
## =========================================================================
Tobit <- function (data, limit=0, limhigh=NULL, transform=NULL, log=FALSE, ...)
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
## ==========================================================================
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
## =======================================================================
leverage <- function (object)
{
  ## Purpose:  extract leverages
  ## Author: Werner Stahel, Date: 10 Nov 2010, 09:31
  lnaaction <- object$na.action
  lh <- object$leverage
  lnm <- names(i.def(object$residuals, object$resid))
  if (is.null(object$rank)) object$rank <- length(coefficients(object))
  if (is.null(lh)) {
    lqr <- object$qr
    if (is.null(lqr)) {
      if (length(lx <- object[["x"]])==0) {
        if (inherits(object, "regrMer"))
          lx <- model.matrix(terms(object$model),object$model)
        else {
          if (inherits(object, "lm"))
            lx <- model.matrix(object)
          else if(length(ld <- object$model))
            lx <- model.matrix(formula(object), data=object$model)
        }
      }
      if (length(lx))
        lqr <- qr(if (length(lwgt <- object$weights)) sqrt(lwgt)*lx else lx)
    }
    if(length(lqr)) {
      lh <- setNames(hat(lqr, intercept=FALSE), lnm)
      ##  lh <- naresid(object$na.action, lh)
    }
  }
  if (length(lh)==0) {
    lmf <- object$model
    if (length(lmf)) {
      lx <- try(model.matrix(formula(object), data=lmf))
    } else {
      lcl <- object$call
      lcl <- lcl[setdiff(names(lcl),"dist")]
      names(lcl)[2] <- "object"
      lcl[1] <- list(quote(model.matrix))
      lx <- try(eval(lcl))
    }
    if (class(lx)!="try.error") 
      lh <- setNames(hat(lx, intercept=FALSE), lnm)
  }
  if (length(lh)==0) message("*leverage* no model.matrix found -> no leverages",
              ". call fitting function with 'x=TRUE'")
  
  if (length(lh)) lh <- naresid(lnaaction, lh)
  lh
}
## ==========================================================================
i.stres <-
  function (x, residuals=NULL, sigma=x$sigma, weights=NULL,
           leveragelim = c(0.99, 0.5))
{ ## sigma=x$sigma, 
  ## Purpose:  calculate  hat  and  standardized residuals
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  1 Mar 2018, 15:45
  leveragelim <- i.def(leveragelim, plgraphics::ploptions("leveragelim"))
  lnaaction <- x$na.action
  if (is.null(residuals)& length(lrs <- x$resid))
    residuals <- naresid(lnaaction, lrs)
  if (length(residuals)==0)  residuals <- residuals(x)
  lres <- as.data.frame(residuals)
  lmres <- ncol(lres)
  if (is.null(weights)) weights <- naresid(lnaaction, x$weights)
  sigma <- i.def(i.def(sigma, x$sigma), x$scale)
  if (length(sigma)==0 && inherits(x, "lm") && !inherits(x, "glm"))
    sigma <-
      if (inherits(x, "mlm")) sapply(summary(x), function(z) z$sigma)
      else summary(x)$sigma
  if (is.null(sigma)||(!all(is.finite(sigma)))||any(sigma<=0)) {
    if (!(inherits(x, "glm") && x$family$family%in%c("binomial","poisson") ||
          inherits(x, "polr")))
      warning(":i.stres: sigma is missing or <=0. I use sigma=1")
    sigma <- 1
  }
  sigma <- rep(sigma, length=lmres)
  ## --- leverage
  llev <- x$leverage
  llev <- if (length(llev)) naresid(lnaaction, llev) else leverage(x)
  if (length(llev)==0)
    return(list(leverage = NULL,
                stresiduals = sweep(lres,2,sigma,"/"), strratio = 1/sigma))
  ##
  llevlim <- leveragelim[1]
  if (mean(llev,na.rm=TRUE)>llevlim) {
    warning(":i.stres: leveragelim not suitable. I use 0.99")
    llevlim <- 0.99
  }
  ## --- standardized residuals
  if (length(lres)==0) {
    warning(":regr/i.stres: no residuals found -> no standardized res.")
    return( list(leverage=llev) )
  }
  if (length(llev)!=nrow(lres)) {
    warning(":regr/i.stres: no leverages available, I set them 0")
    llev <- rep(0,nrow(lres))
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
  lstres <- as.data.frame(as.matrix(lres)*lstrratio)
  for (lj in 1:lmres) {
    if (length(lcq <- attr(lres[,lj], "condquant")))
      attr(lstres[,lj], "condquant") <-
        cbind(lcq[,1:4]*lstrratio[lcq[,"index"]], lcq[,5:6])
  }
  list(leverage = llev, stresiduals = structure(lstres, weighted=lIwgt),
       strratio = lstrratio)
}
## ===================================================================
simresiduals <- function (object, ...)  UseMethod("simresiduals")
  ## Purpose:   simulate residuals according to regression model
  ##            by permuting residuals of actual model or by random numbers
  ## ----------------------------------------------------------------------
  ## Arguments:  simfunction: how are residuals generated?
## ---------------------------------------------------------------------
simresiduals.glm <- function (object, nrep=19, simfunction=NULL,
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
      if (length(lnaaction)) ldata <- ldata[-lnaaction,]
  }
  lfit <- object$fitted.values
  ## prepare call
  lform <- update(formula(object), .Y. ~.)
  lynm <- all.vars(lform[[2]])
  environment(lform) <- environment()
  lcl <- call("glm", formula=lform, data=as.name("ldata"),
              family=object$family, start=object$coef, model=FALSE,
              y=FALSE) ## na.action
  lfam <- object$distrname
  if (length(lfam)==0) lfam <- object$family$family
  if (is.null(lfam)) lfam <- ""
  ly <- object$y
  if (is.null(ly)) ly <- object$response
  if (is.null(ly)) {
    warning(":simresiduals: response not found. No simulated residuals")
    return(NULL)
  }
  ln <- nrow(ldata)
  lone <- rep(1,ln)
  if (!is.function(simfunction)) {
    if(lfam%in%c("binomial","quasibinomial")) {
      if (NCOL(ly)==1 && length(unique(ly))!=2) {
        warning(":simresiduals: binomial distribution with ",
                "unsuitable response.\n  No residuals simulated")
        return(NULL) ## list(simres=numeric(0))
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
        return(NULL) ## list(simres=numeric(0))
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
          residuals.regrpolr(lrs)[,"random"] else
          residuals(lrs, type=glm.restype)
  }
  naresid(lnaaction, lsimres)
}
## ======================================================================
simresiduals.default <-
  function (object, nrep=19, simfunction=NULL, stresiduals=NULL,
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
  lnaaction <- object$na.action
  ldata <- object$allvars ## NAs dropped
  if (is.null(ldata)) {
      ldata <- if (u.debug())  eval(lcall$data)  else
      try(eval(lcall$data))
      if (class(ldata)=="try-error"||is.null(dim(ldata))) {
        warning(":simresiduals: data not found -> No simulated residuals")
        return(NULL)
      }
      if (length(lnaaction)) ldata <- ldata[-lnaaction,]
  }
##-   if ("omit"==class(lnaaction <- object$na.action)) {
##-     lkeep <- rep(TRUE, nrow(ldata))
##-     lkeep[lnaaction] <- FALSE
##-     ldata <- ldata[lkeep,]
##-   }
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
  naresid(lnaaction, lsimres)
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
    location + scale * log(-log(1-p))
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
## pseudoreplicate variability
xdistResdiff <-
  function (object, perc=c(3,10,80), trim=0.1, nmax=100, out="aggregate") ##nsim=100 
{
  ## Purpose:   distance in x space and absolute residual difference
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 13 Oct 2011, 09:40
  if (!inherits(object,"lm")) stop("only suitable for lm like objects")
  lnaaction <- object$na.action
  if (length(lnaaction)) {
    class(lnaaction) <- "omit"
    object$na.action <- lnaaction
  }
  ## this function works with 'short' vectors (without NA elements)
  lstres <- object$stres 
  if (is.null(lstres)) {
    lsr <- try(i.stres(object), silent=TRUE)
    if (class(lsr)=="try-error") {
      warning(":xdistResdiff: no standardized residuals. I use raw residuals")
      lstres <- residuals(object)
    } else  lstres <- lsr$stresiduals
  }
  lqr <- object$qr
  ln <- nrow(lqr$qr)
  lq <- qr.qy(lqr, diag(1, nrow = ln, ncol = lqr$rank)) # like in hat
  li <- which(!is.na(lstres))
  if (ln>nmax) {
    li <- sample(1:ln, nmax, replace=FALSE) # [replace=FALSE]
    lq <- lq[li,]
    lstres <- lstres[li,, drop=FALSE]
    ln <- nmax
  }
  ldist <- dist(cbind(lq))
  lio <- order(ldist)
  ## id's
  lm <- diag(ln)
  lnm <- row.names(lstres)
  lid <- cbind(id1=lnm[rep(1:(ln-1),(ln-1):1)],
               id2=lnm[row(lm)[row(lm)>col(lm)]])
  ## ---
  lrd <- apply(lstres, 2, function(r) as.dist(abs(outer(r,r,"-"))) )
##-   for (lj in 1:ncol(lstres)) {
##-     lrs <- lstres[,lj]
##-     lrd <- abs(outer(lrs,lrs,"-"))
##-     if (nsim) {
##-       lrsim <- matrix(NA,length(ldist),nsim)
##-       for (ls in 1:nsim) {
##-         li <- sample(ln)
##-         lrsim[,ls] <- as.dist(lrd[li,li])
##-       }
##-     }
##-     lrd <- as.dist(lrd)
  rr <- data.frame(lid[lio,], xdist=ldist[lio], resdiff=lrd[lio])
  ##-  if (nsim) rr <- cbind(rr, rdsim=lrsim[lio,])
  class(rr) <- c("xdistResdiff", "data.frame")
  if (out=="aggregate")  xdistResscale(rr, perc=perc)  else rr
}
## ====================================================================
xdistResscale <- function (x, perc=c(3,10,90), trim=1/6)
{
  ## Purpose:  aggregate  xdistResdiff  data
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 18 Oct 2011, 08:40
  if (!inherits(x,"xdistResdiff"))
    stop ("only programmed for  xdistResdiff  objects")
  lxd <- x$xdist
  llim <- c(-0.001*max(lxd), lxd[perc*nrow(x)/100], max(lxd))
  lgrp <- list(cut(lxd, llim))
  lrd <- x[,-(1:3), drop=FALSE]
  lrda <- aggregate(sqrt(lrd), lgrp, mean, trim=trim)[,-1]
  lmn <- apply(sqrt(lrd), 2, mean, trim=trim)
##-  ljsim <- FALSE
##-   if (ncol(x)>5) {
##-     ljsim <- TRUE
##-   lrdsim <-
##-     aggregate(sqrt(as.matrix(x[,-(1:4),drop=FALSE])), lgrp, mean, trim=trim)[,-1]
##-     lnsim <- ncol(lrdsim)
##-     lrdmn <- apply(lrdsim,1,mean)
##-     lrdse <- apply(lrdsim,1,sd)
##-     lrss <- apply(sweep(lrdsim,1,lrdse,"/")^2,2,sum)
##-     lrss0 <- sum((lrd/lrdse)^2)
##- ##-     ltestall <- sum(((lrd-lmn)/lrdse)^2)
##- ##-     lpv <- pchisq(ltestall, df=length(lrd)-1, lower=FALSE)
##- ##-     lpv1 <- pnorm((lrd[1]-lmn)/lrdse[1]) # one-sided is good
##-     lpv <- apply(lrdsim<lrd,1,sum)/lnsim
##-     names(lpv) <- paste("pv",1:length(lpv),sep="")
##-     lpv <- c(lpv, pv.rssq=sum(lrss0<lrss)/lnsim)
##-   }
  rr <- cbind(xdist=aggregate(lxd, lgrp, mean)[,2], rdMean=lrda)
##-   if (ljsim) rr <- cbind(rr, resd.simmean=lrdmn, resd.se=lrdse)
  attr(rr,"limits") <- llim
  attr(rr,"resdMean") <- lmn
  attr(rr,"trim") <- trim
  attr(rr,"perc") <- perc
  ## !!! calculate se and se from resample
##-   if (ljsim) attr(rr,"pvalues") <- lpv #c(shortdist=lpv1, overall=lpv)
  class(rr) <- c("xdistResscale", "matrix")
  rr
}
## =======================================================================
plot.xdistResscale <- function (x, lwd=2, cex=2, xlab="distance in x space",
       ylab="average abs. residual difference", col.aux="grey30", ...)
{
  ## Purpose:   plot average residual difference^2 vs. x distance
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 14 Oct 2011, 08:46
  lxdist <- sqrt(x[,"xdist"])
  lrd <- x[,-1,drop=FALSE]
  lse <- attr(x,"se")
  lIse <- length(lse)>0
  if (lIse) ly2 <- lrd+2*lse
  lymax <- if (lIse) max(ly2) else max(lrd)
  llim <- attr(x,"limits")
  ##
  plot(lxdist, lrd, xlab=xlab, ylab=ylab,type="b", cex=cex,
       xlim=c(0,sqrt(max(llim))), xaxs="i", yaxs="i", ylim=c(0,1.02*lymax),
       axes=FALSE, ...)
  box()
  axis(2)
  lxat <- pretty(c(0,0.8*last(llim)))
  axis(1, at=sqrt(lxat), labels=format(lxat))
  if (lIse) {
    segments(lxdist,lrd-2*lse,lxdist,lrd-2*lse, lwd=lwd, col=col.aux)
    lysim <- attr(x,"resd.simmean")
    if (length(lysim)) {
      lxd <- sqrt(max(llim))/(nrow(x)*10)
      segments(lxdist-lxd, lysim, lxdist+lxd, lysim, col=col.aux)
    }
  }
##-   axis(3,at=c(0,sqrt(llim[-1])),labels=rep("",length(llim)), col=col.aux)
  abline(h=attr(x,"resdMean"), lty=3, col=col.aux)
  abline(v=sqrt(llim[2:(length(llim)-1)]), lty=3, col=col.aux)
  invisible(NULL)
##  "plot.xdistResscale done"
}
## ============================================================================
