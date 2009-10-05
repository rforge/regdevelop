#  regr.R  Functions that are useful for regression, W. Stahel,
##  April 2008
##  source("/u/stahel/R/regr/regr.R")
##  d.blast <- read.table(blast.dat",header=T)
##  d.blast$location <- factor(d.blast$location)
## ==========================================================================
regr <- function(formula, data, tit=NULL, method="lm", family=NA,
  init.reg="f.ltsreg", subset=NULL, weights=NULL, na.action=nainf.exclude,
  model=TRUE, calcdisp=NULL, termtable=TRUE, vif=TRUE, ...)
{
  ## !!! dispersion: allow to be set.
  ## Purpose:    fit all kinds of regression models
  ## -------------------------------------------------------------------------
  ## Arguments:
  ##   formula, data, ...  as with lm
  ##   tit        title (becomes tit attribute of result)
  ##   method     method to be called for fitting
  ##   family     family argument of  glm  (only used if method == "glm")
  ##   init.reg   function used to initialize robust methods
  ##   calcdisp   should dispersion be calculated for
  ##              family=binomial and family=poisson
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  Jan 02
  ## -------------------------------------------------------------------------
## preparations
  lcontrol <- list(calcdisp=calcdisp, suffmean=3)
  if (length(method)==0||is.na(method)) method <- "lm"## data
  dataname <- deparse(substitute(data))
  lcall <- match.call()
  lform <- as.formula(formula)
  environment(lform) <- environment()
##-   lcall$formula <- lform
  ldata <- data
    lnobs <- nrow(ldata)
  if (lnobs==0) stop("!regr! no observations in data")
## variables used
  lvars <- all.vars(lform)
  if (method=="nls") {
    lvars <- lvars[match(lvars,names(ldata),nomatch=0)>0]
    if (length(lvars)==0)
      stop("!regr (method='nls')! no variables of  data  used in formula")
  } else {
    lj <- match(".",lvars,nomatch=0)
    if (lj>0) lvars <- unique(c(lvars[-lj],names(ldata)))
  ## checking variables -- potentially not needed    
    if (any(lvj <- match(lvars,names(ldata),nomatch=0)==0)) {
    lvwrong <- lvconst <- NULL
    for (lv in lvars[lvj]) {
      lvv <- try(get(lv),silent=TRUE)
      if (class(lvv)[1]=="try-error") lvwrong <- c(lvwrong, lv) else {
        if (NROW(lvv)!=lnobs) {
          if (NROW(lvv)!=1) lvwrong <- c(lvwrong, lv) else
          lvconst <- c(lvconst, lv) } else {
            if (is.matrix(lvv))
      stop("!regr! Matrices not allowed in formula. Use cbind(name1,name2,...)")
            ldata[[lv]] <- lvv
          }
      }
    }
    if (length(lvconst)>0) {
      warning (paste(":regr: scalar constants used in formula:",
          paste(lvconst,collapse=", ")))
      lvars <- lvars[match(lvars,lvconst,nomatch=0)==0]
    }
    if (length(lvwrong)>0) 
      stop (paste("!regr! invalid variables in formula:",
          paste(lvwrong,collapse=", ")))
  }
}
  ldt <- ldata[,lvars,drop=FALSE]
  lvy <- all.vars(lform[[2]])
  ly <- cbind(ldt[,lvy])   ## not ,drop=F !
  ## subset
  if (length(subset)) {
    lsubset <- eval(substitute(subset), ldata)
    lsubset <- nna((1:lnobs)[lsubset]) ## logicals and negative integers
    if (length(lsubset)==0) stop("!regr! no data left after subsetting")
    ldrop <- (1:lnobs)[-lsubset]
    if (length(ldrop)) ly[ldrop,] <- NA
    ldt[,lvy] <- ly
  }
## weights
  lwgts <- NULL
  lwgt <- length(weights)>0
  if (lwgt) {
    lwgts <- eval(weights, data)
    if (length(weights)!=nrow(data))
      stop("!regr! weights and data have incompatible dimensions")
    ldt[[".Weights"]] <- lwgts ## needed?
  }
## missing values
  ldta <- na.action(ldt)
  if (NROW(ldta)==0) stop("!regr! no observations without NAs")
  .Weights <- ldta[[".Weights"]]
  ## as.character(lform[2])==lvars[1] ## left side is simple
## strange variables
  l1v <- sapply(ldta, function(x) all(x==c(x[!is.na(x)],0)[1],na.rm=TRUE) )
                                ## covers case of several or all NAs
  if (any(l1v)) {
    warning(paste(":regr: variable(s)", paste(lvars[l1v],collapse=", "),
                  "has (have) no distinct values")) #  -> dropped.
##-     lvars <- lvars[!l1v]
##-     stop("bug: elimination of const.variables not yet programmed")
##-                                         # drop from lform!
  }
## type of variables, convert binary factors to numeric
  ldta <- datatype(ldta)
  lvartype <- attr(ldta,"vartype")
  lbinlevels <- attr(ldta,"binlevels")
  ## family 
  lfam <- as.character(substitute(family))[1]
  if (is.na(lfam)||lfam=="") 
    lfam <- c(c="gaussian",f="multinomial",o="ordered",
              b="binomial")[lvartype[1]]
  if (substring(lfam,1,7)=="multinom") lfam <- "multinomial"
  lfitname <- c(gaussian="i.lm", multinomial="i.multinomial",
                ordered="i.polr")[lfam]
  if (is.na(lfitname)) lfitname <- "i.glm"
  if (!exists(lfitname))
    stop (paste("!regr! fitting function",lfitname, "not found"))
  lfitfun <- get(lfitname)
  ## ---------------------
  lreg <- lfitfun(formula=lform, data=ldta, wgts=.Weights,
                  family=family, fname=lfam, method=method, 
                  control=lcontrol, model=model, vif=vif,
                  termtable=termtable, 
                  ...)
  ## ---------------------
  lreg$funcall <- lreg$call
  lreg$call <- lcall
    ltit <- if (length(tit)==0) attr(data,"tit") else tit
    ldoc <- attr(data,"doc") 
    tit(lreg) <- ltit
    doc(lreg) <- ldoc
  lterms <- if (method=="nls") NULL else terms(lreg)
  if (is.null(attr(lterms, "predvars"))&method!="nls") { ## needed for survreg
    lterms <- attr(lm(lform, data=ldta, method='model.frame'),'terms')
    attr(lreg$terms,'predvars') <- attr(lterms,'predvars')
  }
##-   ## expand for subset and NAs
##-   if (length(lreg$na.action))
##-     if(length(ldrop)&class(lreg$na.action)=="exclude") {
##-     li <- naresid(ldrop,naresid(lreg$na.action,1:length(lreg$resid)))
##-     lex <- which(is.na(li))
##-     class(lex) <- "exclude"
##-     lreg$na.action <- lex
##-   }
#  lnaac <- lreg$na.action
  lresnm <- colnames(lreg$residuals)
  ## leverage
  lhat <- pmax(0,hat(lreg$qr))
##-   if (length(lnaac)) lhat <- naresid(lnaac,lhat)
##-     lhat <- u.merge(lhat,NA,lnnahat, names=lrl)
  lreg$h <- lhat
  ## standardized res
  if (length(lhat)==NROW(lreg$stres))
    lreg$stres <- lreg$stres/ifelse(lhat>=1,1,sqrt(pmax(1e-10,1-lhat)))  else {
      warning(":regr: bug: leverages not available")
      lhat <- NULL
    }
  ## cov of estimates
  lsig <- lreg$sigma
  if (length(lsig)==0)
    lsig <- if (!is.null(lreg$dispersion)) sqrt(lreg$dispersion) else 1
  lcov <- lreg$cov.unscaled
  if (length(lcov)>0) {
    ldfr <- lreg$df.residual
    if (is.null(ldfr))  ldfr <- lreg$df[2]
    if (is.na(ldfr)|ldfr==0) {
      warning(':regr: residual df are 0')
      lsig <- 0
    }
    if (length(lsig)==1) {
      lcov <- lcov*lsig
      lreg$covariance <- lcov
    }
    lse <- sqrt(diag(lcov))
    lreg$correlation <- lcov/outer(lse, lse)
  }
  ## misc
  lreg$response <- ly
  lreg$vartype <- lvartype
  lreg$binlevels <- lbinlevels
  if (is.null(lreg$familyname)) lreg$familyname <- lfam
  lreg$fitfun <- class(lreg)[1]
  lreg$na.action <- attr(ldta,"na.action")
##-   class(lreg) <- if (class(lreg)[1]=="try-error") class(lreg)[-1] else
  class(lreg) <- c("regr",class(lreg))
## result
  lreg  
}
## -----------------------------------------------------------------------
datatype <- function(data)
## type of variables, convert binary factors to numeric
{
## convert character to factor
  lvars <- names(data)
  for (lvn in lvars) {
    lv <- data[[lvn]]
    if (is.character(lv)) data[[lvn]] <- factor(lv)
  }
 ## analyze variables
  ltype <- rep("c",length(lvars))
  names(ltype) <- lvars
  lbinlevels <- list()
  for (lvn in lvars) {
    lv <- data[[lvn]]
    if (length(unique(lv))==2) {
      ## binary variable
      ltype[lvn] <- "b"
      lvf <- factor(lv)
##-       if (is.factor(lv)) {
##-         llv <- levels(data[lvn])
      data[[lvn]] <- as.numeric(lvf)-1  # else
##-       if (is.logical(lv)) data[[lvn]] <- as.numeric(lv)
      lbinlevels <- c(lbinlevels, list(levels(lvf)))
        names(lbinlevels)[length(lbinlevels)] <- lvn
    } else  if (is.factor(lv)) {
      ltype[lvn] <- if (is.ordered(lv)) "o" else "f"
      ## repair factor bug, since otherwise, rlm does not work
      ## ... unless model.matrix is called with ...
      data[[lvn]] <- factor(lv)  
    }
  }
  attr(data,"vartype") <- ltype
  attr(data,"binlevels") <- lbinlevels
  data
}
## -----------------------------------------------------------------------
i.lm <- function(formula, data, wgts=data[[".Weights"]],
                 family, fname="gaussian", method="ls", control=NULL,
                 model=TRUE, vif=TRUE, termtable=TRUE, ...)
{
  ## Purpose:  internal: fit lm
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
  .Weights <- wgts 
  lreg <- if (method[1]=="rlm") {
    require(MASS)
    lmeth <- c(method,'MM')[2]
    rlm(formula, data=data, weights=.Weights, model=model, method=lmeth,...)
  } else if (method[1]=="nls") {
    nls(formula, data=data, weights=.Weights, model=model, ...)  
  } else lm(formula, data=data, weights=.Weights, model=model, ...)
  lreg$call$formula <- formula
  if (class(lreg)[1]=="mlm")
    return(i.mlmsum(lreg, wgts, model, vif, termtable))
  lreg1 <- summary(lreg)
  ly <- lreg$model[[1]]
  lsdy <- sqrt(var(ly))
  lsig <- lreg1$sigma
  lreg$sigma <- lsig
  ## ---
  if (is.finite(lsig)&&lsig>0) {
    lreg$stres <- lreg$residuals/lsig
    if (length(.Weights)) lreg$stres <- lreg$stres*sqrt(.Weights)
  }
  ## leverage not taken into account here
  ldfr <- lreg$df.residual
  if (length(ldfr)==0||is.na(ldfr)) lreg$df.residual <- ldfr <- lreg1$df[2] # rlm
  lcomp <- c("r.squared","fstatistic","colregelation","aliased",
             "df","cov.unscaled")
  lreg[lcomp] <- lreg1[lcomp]
  lreg$adj.r.squared <- 1-(1-lreg$r.squared)*(length(lreg$residuals)-1)/ldfr
  ## --- deviances
  if (termtable&method!="nls") {
    ltt <- i.termtable(lreg, lreg1$coef, "F", lsdy=lsdy, vif=vif)
    lreg$testcoef <- ltt$test
    lreg$allcoef <- ltt$allcoef
  }
  ## result
  lreg
}
##- ## -----------------------------------------------------------------------
##- i.mlm <- function(formula, data, wgts=data[[".Weights"]],
##-                  family, fname="gaussian", method="ls", control=NULL,
##-                  model=TRUE, vif=TRUE, termtable=TRUE, ...)
##- {
##-   ## Purpose:  internal: fit multivariate lm
##-   ## ----------------------------------------------------------------------
##-   ## Arguments:
##-   ## ----------------------------------------------------------------------
##-   ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
##- ##-   .Weights <- wgts 
##- ##-   lreg <- lm(formula, data=data, weights=.Weights, model=model, ...)
##- ##-   lreg$call$formula <- eval(lreg$call$formula) # patch
##- ##-   lreg1 <- summary(lreg)
##-   lts <- ltp <- NULL
##-   for (ly in 1:length(lreg1)) {
##-     lrg <- lreg1[[ly]]
##-     lts <- cbind(lts,c(lrg[["sigma"]],lrg[["r.squared"]],
##-                        lrg[["fstatistic"]]))
##-     ltp <- cbind(ltp,lrg[["coefficients"]][,4])
##-   }
##-   lmodel <- nrow(lts)>=5  # non-trivial model
##-   if (lmodel) {
##-     lts[4,] <- 1-pf(lts[3,],lts[4,],lts[5,])
##-     lts <- lts[1:4,]
##-   }
##-   dimnames(lreg$coefficients)[[2]] <- as.character(formula[[2]])[-1]
##-   dimnames(ltp) <- dimnames(lreg$coefficients)
##-   dimnames(lts) <- list(rep(c("sigma","r.squared","fstatistic","p-value"),
##-                             length=nrow(lts)), dimnames(ltp)[[2]])
##-   lreg$pvalues <- ltp
##-   lreg$stats <- lts
##-   lreg$sigma <- lsig <- lts["sigma",]
##-   lres <- residuals(lreg)
##-   if (all(lsig>0)) {
##-     lreg$stres <- sweep(lres,2,lsig,"/")
##-     if (length(.Weights)) lreg$stres <- lreg$stres*sqrt(.Weights)
##-   }
##-   lreg$resmd <- mahalanobis(lres,0,var(lres))
##-   ldfr <- lreg$df.residual
##-   lreg$r.squared <- lr2 <- lts["r.squared",]
##-   lreg$adj.r.squared <- 1-(1-lr2)*(nrow(lreg$resid)-1)/ldfr
##-   lcomp <- c("aliased","df","cov.unscaled")
##-   lreg[lcomp] <- lreg1[[1]][lcomp]
##-   lreg$drop1 <- if (lmodel) drop1.mlm(lreg)
##- #  class(lreg) <- c("mregr","mlm","lm")
##-   lreg
##- }
## -----------------------------------------------------------------------
i.mlmsum <- function(object, weights=NULL, model=TRUE, vif=TRUE, termtable=TRUE)
{
  ## Purpose:  internal: fit multivariate lm
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
    lts[4,] <- 1-pf(lts[3,],lts[4,],lts[5,])
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
    if (length(weights)) object$stres <- object$stres*sqrt(weights)
  }
  object$resmd <- mahalanobis(lres,0,var(lres))
  ldfr <- object$df.residual
  object$r.squared <- lr2 <- lts["r.squared",]
  object$adj.r.squared <- 1-(1-lr2)*(nrow(object$resid)-1)/ldfr
  lcomp <- c("aliased","df","cov.unscaled")
  object[lcomp] <- lreg1[[1]][lcomp]
  object$drop1 <- if (lmodel) drop1.mlm(object)
#  class(lreg) <- c("mregr","mlm","lm")
  object
}
## -----------------------------------------------------------------------
i.glm <- function(formula, data, wgts=data[[".Weights"]],
                 family, fname="gaussian", method="ls", control=NULL,
                 model=TRUE, vif=TRUE, termtable=TRUE, ...)
{
  ## Purpose:  internal: fit glm
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
  lfamily <- get(fname)
##-   environment(formula) <- environment()
  .Weights <- wgts 
  lreg <- glm(formula, family=lfamily, data=data, weights=.Weights, model=model, ...)
  lreg1 <- summary(lreg)
  lcoeftab <- lreg1$coef
  ly <- lreg$model[,1]
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
        lreg$familyname <- paste("quasi",fname,sep="")
        lcoeftab[,2] <- lcoeftab[,2]*sqrt(ldisp)
        lcoeftab[,3] <- lcoeftab[,3]/sqrt(ldisp)
##-         lcoeftab[,4] <- 2*pnorm(lcoeftab[,3],lower.tail=FALSE)
      }
      else ldisp <- 1
    }
  }  # else calcdisp <- FALSE
  lreg$dispersion <- ldisp
  lreg$sigma <- sqrt(ldisp)
  ## ---
  if (ldisp>0) lreg$stres <- residuals(lreg,"pearson")/sqrt(ldisp)
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
  ltbd[1:2,3] <- 1-pchisq(ltbd[1:2,1],ltbd[1:2,2])
  if (!lsuffmean) ltbd[2,3] <- NA
  lreg$devtable <- ltbd
  ## ---
  if (termtable) {
    ltt <- i.termtable(lreg, lcoeftab, ltesttype, lsdy=1, vif=vif)
    lreg$testcoef <- ltt$test
    lreg$allcoef <- ltt$allcoef
  }
  ## result
  lreg
}
## -----------------------------------------------------------------------
i.multinomial <- function(formula, data, wgts=data[[".Weights"]],
                 family, fname="gaussian", method="ls", control=NULL,
                 model=TRUE, vif=TRUE, termtable=TRUE, ...)
{
  ## Purpose:  internal: fit multinom
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
#  ltr <- control$trace
#  if (length(ltr)==0) ltr <- trace
  require(nnet)
  .Weights <- wgts 
  lreg <- multinom(formula, data=data, weights=.Weights, model=TRUE,
                   trace=FALSE, ...)
  lreg$call$formula <- formula
  lreg1 <- summary(lreg)
  lreg$dispersion <- lreg$sigma <- 1
  lres <- lreg1$residuals
  lreg$residuals <- lres
  lcf <- lreg1$coefficients
  lreg$coefficients <- lcf
  lreg$aic <- lreg1$AIC
  ldfm <- lreg1$edf-nrow(lcf)
  lreg$df <- c(ldfm,prod(dim(lres)-1)-ldfm,ldfm)
##-   environment(lreg$call$formula) <- environment()
  ldr1 <- try(drop1(lreg, test="Chisq", trace=FALSE),
              silent=TRUE)
  if (class(ldr1)[1]=="try-error") {
    warning(paste(":regr/i.multinom: drop1 did not work.",
                  "I return the multinom object"))
    return(lreg)
  } else {
  ldr1 <- ldr1[-1,]}
  names(ldr1) <- c("df", "AIC", "Chisq","p.value") #if(calcdisp) "F" else
  lreg$testcoef <- lreg$drop1 <- ldr1
  class(lreg) <- c(class(lreg1),class(lreg))
  lreg
}
## -----------------------------------------------------------------------
i.polr <- function(formula, data, wgts=data[[".Weights"]],
                 family, fname="gaussian", method="ls", control=NULL,
                 model=TRUE, vif=TRUE, termtable=TRUE, ...)
{
  ## Purpose:  internal: fit ordered y
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 11:18
  require(MASS)
  .Weights <- wgts 
  lfo <- formula
  environment(lfo) <- environment()
  lreg <- # if (length(data[,".Weights"])>0)
    polrregr(lfo, data=data, weights=.Weights, Hess=TRUE, ...)
#      else  polr(lfo, data=data, Hess=TRUE, ...)
  lreg$call$formula <- lfo
  lreg1 <- try(summary(lreg))
  if (class(lreg1)[1]=="try-error") {
    warning(paste(":regr/i.polr: summary did not work.",
                  "I return the polr object"))
#    lreg$call$data <- call$data
    class(lreg) <- c("try-error","polr")
    return(lreg)
  }   ## ---
  lcf <- lreg1$coefficients
  lreg$intercepts <- lcf[(lreg1$pc+1):nrow(lcf),1:2]
  lreg$stres <- NULL
  ## --- deviances
  if (termtable) {
    ltt <- i.termtable(lreg, lreg1$coef, ltesttype="Chisq", lsdy=1, vif=vif)
    lreg$testcoef <- ltt$test
    lreg$allcoef <- ltt$allcoef
  }
  lreg$dispersion <- 1
#  lreg$call$data <- call$data
  ## result
  lreg
}
## -----------------------------------------------------------------------
i.termtable <- function(lreg, lcoeftab, ltesttype="F", lsdy, vif=TRUE) 
{
  ## Purpose:  generate term table for various models
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  4 Aug 2004, 15:37
  if(length(attr(terms(lreg),"term.labels"))==0) 
    return(list(test=data.frame(coef=lreg$coef,stcoef=NA,signif=NA,
                  R2.x=NA,df=length(lreg$coef),p.value=NA)))
  ldr1 <-if (is.R())
    try(drop1(lreg, test=ltesttype, scope=terms(lreg),trace=FALSE),
        silent=TRUE) else
    try(drop1(lreg, test=ltesttype, scope=terms(lreg),trace=FALSE))
  if (class(ldr1)[1]=="try-error") {
    warning(paste(":regr: drop1 did not work. I return the (g)lm table"))
    lsum <- summary(lreg)
    lcft <- lsum$coef
    if (length(lcft)==0) lcft <- lsum$parameters ## nls
    return(list(test=lcft)) # !!! noch reparieren
  }
#  ldr1 <- ldr1[-1,]
  ldr1$RSS <- NULL # same ncol for lm and glm
##-   if (lglm&&!is.R()) {
##-     ldr1 <- data.frame(ldr1,NA,1-pchisq(ldr1[[2]],ldr1[[1]]))
##-     if (nrow(ldr1)==1) row.names(ldr1) <- attr(terms(lreg),"term.labels")
##-   }
  if (inherits(lreg,"rlm"))  ldr1[,4] <- ldr1[,2]/ldr1[,1]
  if (inherits(lreg,"mlm")||inherits(lreg,"manova"))
    return(list(test=ldr1))  ## !!! needs much more
##-   ldfr <- lreg$df[2] # !!! df.residual
  ldfr <- lreg$df.residual
  if (length(ldfr)==0) ldfr <- lreg$df[2]
##  ldfy <- sum(lreg$df[1:2])
  lcoef <- lreg$coefficients
  if (ncol(lcoeftab)==3)
    lcoeftab <- cbind(lcoeftab,pt(-abs(lcoeftab[,3]),ldfr)*2)
  ldr1 <- ldr1[-1,]
  lvif <- rep(NA,nrow(ldr1)-1)
  if (vif) {
    lvift <- try(vif.regr(lreg))
    if (class(lvift)[1]=="try-error") {
      warning(":regr/i.termtable: error in the calculation of R2.x's")
    } else lvif <- lvift[,3]
  }
  ## intercept
  if ("(Intercept)"==names(lcoef)[1]) {
    ldr1 <- rbind("(Intercept)"=c(1,NA,NA,NA,lcoeftab[1,4]),ldr1)
    lvif <- c(NA,lvif)
  }
  ltq <- qt(0.975,ldfr) # not quite right for rlm
  lfq <- if (ltesttype=="F") qf(0.95,ldr1[,1],ldfr) else {
    if (ltesttype=="Chisq") qchisq(0.95,ldr1[,1]) else NA }
  ldr1t <- sqrt(max(0,ldr1[,4])/lfq) 
  ltlb <- dimnames(ldr1)[[1]]
  lcont <- ldr1[,1]==1
  lprcol <- pmatch("Pr(",names(ldr1))
  if (!is.na(lprcol))
    ldr1t[!lcont] <- -qnorm(ldr1[!lcont,lprcol]) ## z score for factors
  lnt <- nrow(ldr1)
  ltb <- c(list(coef=rep(NA,lnt), stcoef=rep(NA,lnt), signif=ldr1t), ldr1)
  if (any(lcont)) {
    lclb <- ltlb[lcont]
    licf <- pmatch(lclb,dimnames(lcoeftab)[[1]])
    if (any(is.na(licf))) warning(":regr: bug: coefficients not found")
    ## !!! bug if some coefs have 0 df
    lcf <- lcoef[lclb]
    ltst <- lcoeftab[lclb,3]
    ## standardized coefficients
    lmmt <- model.matrix(lreg)
    lmmcol <- match(lclb,dimnames(lmmt)[[2]],nomatch=0)>0 ## ginge einfacher
    lsd <- if (any(lmmcol)) 
      ifelse(lmmcol,sqrt(apply(lmmt[,lclb[lmmcol],drop=FALSE],2,var)),
             NA) else NA
    lstcf <- lcf*lsd/lsdy
    ## fill in
    ltb$coef[lcont] <- lcf
    ltb$stcoef[lcont] <- lstcf
    ltb[["signif"]][lcont] <- ltst/ltq #sign(lcf)*sqrt(ltb[["F value"]][lcont])
  }
  lr2 <- if (length(lvif)==length(ltb[[1]])) 1-1/lvif else NA
  ltb <- data.frame(ltb[1:3],R2.x=lr2,ltb[c(4,length(ltb))],
                    row.names=ltlb)
  names(ltb)[5:6] <- c("df","p.value")
  ## dummy coef
  lallcf <- if (is.R()) try(dummy.coef(lreg),silent=TRUE)  else
    try(dummy.coef(lreg)) 
  if (class(lallcf)=="try-error") {
    warning("dummy.coef did not work")
    lallcf <- NULL
  }
  list(test=ltb, allcoeff=lallcf)
}
## =======================================================================
##- 
##- 
##- 
##-   if (method=="robust") {
##-     lireg <- init.reg
##-     if (is.character(lireg)) lireg <- get(lireg)
##-     if (!is.function(lireg)) {
##-     lireg <- f.ltsreg
##-     warning("unsuitable argument  init.reg . Will use  f.ltsreg") }
##- ##-     lreg <- mmreg(formula, na.action=na.omit, data=ldata, weights=.Weights,
##- ##-                   init.reg=lireg, ...)
##-     lreg <- mmreg(formula, data=ldt, weights=.Weights, init.reg=lireg, ...)
##-     lreg$qr <- lqr }
##-   if (method=="rlm") {
##-     lqr <- lreg$qr
##-     lreg <- if (length(.Weights)==0) rlm(formula, qr=T, data=ldt,...) #
##-             else rlm(formula, qr=T, weights=.Weights, data=ldt,...) 
##-     lreg$qr <- lqr }
##-   attr(lr0$terms, "formula") <- formula
##- ## ---------------------
##- 
## ==========================================================================
print.regr <- function (x, correlation = FALSE,
    dummycoef = NULL,
    digits = max(3, options("digits")[[1]] - 3), 
    symbolic.cor = p > 4, signif.stars = options("show.signif.stars")[[1]],
    doc = options("doc")[[1]], 
    residuals=FALSE, niterations=FALSE, ...) 
{
##
  if (length(doc)==0) doc <- 0
  if (doc>=1) if (length(tit(x)))
    cat("\n",tit(x),"\n")
  if (doc>=2) if (length(doc(x)))
    cat(paste(doc(x),"\n"),"\n")
  if (inherits(x,"mlm")) { print.mregr(x, ...)
    return() }
  if (length(dummycoef)==0)
    dummycoef <- c(options("show.dummy.coef")[[1]],TRUE)[1]
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),"\n", sep = "")
  cat("Fitting function ",x$fitfun,"\n")
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
    nsingular <- df[3] - df[1]
    if ((!is.na(nsingular))&&nsingular>0) 
        cat("\nCoefficients: (", nsingular,
            " not defined because of singularities)\n", sep = "")
    # else {
    if (!is.null(x$sigma))
      if((!is.finite(x$sigma))||x$sigma<=0)
        cat("\n!!! Error variance is 0 !!!")
  ltc <- x$testcoef
  if (length(ltc)>0) {
    cat("\nTerms:\n")
    lip <- nna(pmatch(c("R2","p."),colnames(ltc)))
    if (length(lip)) ltc[,lip] <- round(ltc[,lip],max(3,digits))
    print(ltc)
      # }
    if (x$familyname=="multinomial") {
      cat("\nCoefficients:\n")
      print(t(x$coefficients))
    } else {
      lidf <- match("df",colnames(x$testcoef))
      if (is.na(lidf)) warning(":print.regr: df of coef not available")
      else {
        mterms <- row.names(x$testcoef)[x$testcoef[,"df"]>1]
        if (length(mterms)>0&dummycoef & length(x$allcoef)>0) {
          imt <- mterms%in%names(x$allcoef)
          mt <- x$allcoef[mterms[imt]]
          if (length(mt)>0) {
            cat("\nCoefficients for factors:\n")
##-             lia <-
##-               apply(outer(row.names(x$testcoef)[mterms],names(x$allcoef),
##-                           function(X,Y) pmatch(X,table=Y)),2,any,na.rm=TRUE)
##-             lia <- rep(FALSE,length(x$allcoef))
##-             for (lt in row.names(x$testcoef)[mterms])
##-               lia[pmatch(lt,names(x$allcoef))] <- TRUE
            print(mt) }
        } else  cat("\n")
      }}
  } else {
    lcf <- x$coef
    if (length(lcf)) {
      cat("\nCoefficients:\n")
      print(x$coef)
    }
  }
  if (length(x$intercepts)) { # polr
    cat("Intercepts:\n")
    print(x$intercepts)
  }
  if (length(x$fstatistic)>0) {
    cat("St.dev.error: ", formatC(x$sigma, digits = digits),
        "  on", rdf, "degrees of freedom\n")
    cat("Multiple R^2: ", formatC(x$r.squared, digits = digits))
    cat("   Adjusted R-squared:", 
        formatC(x$adj.r.squared, digits = digits),
          "\nF-statistic:  ", formatC(x$fstatistic[1], 
            digits = digits), "  on", x$fstatistic[2], "and", x$fstatistic[3], 
            "d.f.,  p.value:", formatC(1 - pf(x$fstatistic[1], 
            x$fstatistic[2], x$fstatistic[3]), dig = digits), 
         "\n")
    }
  if (length(x$deviance)>0) {
    if (length(x$devtable)>0) {
      cat("\n")
      print(x$devtable)
    }
    cat("\nFamily is ",x$familyname, ".  Dispersion parameter ", 
        if (x$dispersion==1) "taken to be 1" else
        paste("estimated to be ", format(x$dispersion)),
        ".\n",sep="")
    cat("AIC: ", format(x$aic, digits = max(4, digits + 1)), "\n", sep = "")
    if (niterations&&length(x$iter)>0)
      cat("Number of iterations:", x$iter, "\n")
  }
  correl <- x$correlation
  if (correlation && length(correl)>0) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (symbolic.cor) 
        print(symnum(correl)[-1, -p])
      else {
        correl[!lower.tri(correl)] <- NA
        print(correl[-1, -p, drop = FALSE], digits = digits, 
              na = "")
      }
    }
  }
  if (length(x$binlevels)>0) {
    cat("\nFactor(s) with two levels converted to 0-1 variable(s):\n")
    print(as.matrix(data.frame(x$binlevels,row.names=0:1)))
  }
#    cat("\n")
  invisible(x)
}
## ==========================================================================
summary.regr <- function(object, dispersion=NULL, ...)  object
## ==========================================================================
drop1.regr <-
  function (object, scope=NULL, scale = 0, test = NULL, k = 2,
            sorted = FALSE, add=FALSE, ...) 
{
  ## Purpose:    drop1/add1 for regr objects
  ## ----------------------------------------------------------------------
  lfam <- object$familyname
  lres <- object$residuals
  if (is.null(test)) test <- if (is.null(lfam)) "none" else {
    if (lfam=="gaussian"|((lfam=="binomial"|lfam=="poisson")&&
          object$dispersion>1)) {
          if (inherits(object,"mlm")) "Wilks" else "F" }
            else "Chisq"
  }
  if (length(scope)==0) {
    scope <- if (add) attr(terms(update.formula(object, ~(.)^2-.)), 
                            "term.labels") else
        drop.scope(object)
  } else 
    if (!is.character(scope)) 
            scope <- attr(terms(update.formula(object, scope)), 
                "term.labels")
  if (length(scope)==0) {
    warning(":drop1/add1.regr! no valid scope")
    ldr1 <- data.frame(Df = NA, "Sum of Sq" = NA, RSS =NA, AIC = NA,
                      row.names = "<none>")
    return(ldr1)
  }
  class(object) <- setdiff(class(object), "regr")
  dr1 <- if (add)
    add1(object, scope=scope, scale=scale, test=test, k=k, ...) else 
    drop1(object, scope=scope, scale=scale, test=test, k=k, ...)
##-   rnm <- row.names(dr1) 
##-   row.names(dr1) <- paste(ifelse(substring(rnm,1,1)=="<","",
##-                                  if (add) "+ " else "- "),rnm,sep="")
  attr(dr1,'drop') <- !add
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
  drop1.regr(object, scope=scope, scale=scale, test=test, k=k,
             sorted=sorted, add=TRUE, ...)
}
## ==========================================================================
predict.regr <- 
function (object, newdata = NULL, se.fit = FALSE, scale = object$sigma, 
    interval = c("none", "confidence", "prediction"), level = 0.95, 
    type = NULL, terms = NULL, na.action = na.pass, 
    ...)
  ## bug: if used with NULL newdata, predictions will be produced
  ## for obs in model.matrix, which excludes those with missing values
{
##-   lglm <- inherits(object,"glm")
##-   lmeth <- object$call$method
##-   lnls <- length(lmeth)>0 && lmeth=="nls"
  if (length(type)==0)
    type <- if (inherits(object,"glm")) "link" else
  if (inherits(object, "polr")) "link" else "response"
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
##-     lpred <- if (lglm)
##-       predict.glm(object, newdata=ldt, type=type, se.fit=se.fit, 
##-                   dispersion=lscale^2, terms=terms, na.action = na.action )
##-     else predict.lm(object, newdata=ldt, type=type, se.fit=se.fit, 
##-                   terms=terms, na.action = na.action, scale=object$sigma )
##-   lpred <- predict(object, newdata=ldt, type=type, se.fit=se.fit, 
##-                   dispersion=lscale^2, terms=terms, na.action = na.action )
##-   lpred
  predict(object, newdata=ldt, type=type, se.fit=se.fit, scale=object$sigma,
          dispersion=object$dispersion^2, terms=terms, na.action = na.action )
}
## ==========================================================================
extractAIC.regr <- function (fit, scale = 0, k = 2, ...) 
{  #  AIC, divided by n to allow for comparing models with different n
  fit$residuals <- nna(fit$residuals)
  class(fit) <- setdiff(class(fit),"regr")
  rr <- extractAIC(fit, scale = scale, k = k, ...) 
  rr
}
## ==========================================================================
vif.regr <- function(mod)
{
  ## Purpose:   vif.lm  of library  car
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 11 Mar 2005, 09:18
    v <- mod$sigma^2* mod$cov.unscaled
    cls <- dimnames(model.matrix(mod))[[2]]%in%dimnames(v)[[2]]
                                        # needed for singular cases
    assign <- attributes(model.matrix(mod))$assign[cls]
    terms <- labels(terms(mod))
    n.terms <- length(terms)
    if (n.terms < 2) {
##-         stop("model contains fewer than 2 terms")
      return(matrix(1,1,3))
    }
    if (length(v)==0) { # ||n.terms!=nrow(v)|nrow(v)!=ncol(v)
      warning(":vif.regr: mod$cov.unscaled  is inappropriate. no vifs")
      return(matrix(NA,n.terms,3))
    }
    if (names(coefficients(mod)[1]) == "(Intercept)") {
        v <- v[-1, -1]
        assign <- assign[-1]
    }
    else warning("No intercept: vifs may not be sensible.")
    R <- cov2cor(v)
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
#        cat("trying -", tt, "\n")
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
summary.mreg <- function(object, ...)
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
  lts[4,] <- 1-pf(lts[3,],lts[4,],lts[5,])
  lts <- lts[1:4,]
  dimnames(ltp) <- dimnames(object$coefficients)
  dimnames(lts) <- list(c("sigma","r.squared","fstatistic","p-value"),
                        dimnames(ltp)[[2]])
  list(coefficients=object$coefficients, pvalues=ltp, stats=lts)
}
## =================================================================
drop1.mlm <- 
function (object, scope = NULL,
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
    if (!(inherits(object, "maov")||inherits(object, "mlm"))) 
        stop("object must be of class \"manova\", \"maov\", or \"mlm\"")
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
    ldata <- get(as.character(object$call$data),
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
    if (rss.qr$rank < ncol(res)) 
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
      eigs <- Re(eigen(qr.coef(rss.qr, bss),symmetric = FALSE)$values)
      stats <- rbind(stats, "<total>"=tstfn(eigs, object$df[1], rdf))
    }
    data.frame(stats, p.value = pf(stats[,2],stats[,3],stats[,4],
                           lower.tail = FALSE))
#    attr(stats,"df") <- ldf
#    stats
}
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
  ##  determine rows with NA's in model.frame for expanded model
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
        pgumbel <- function(q) exp(pweibull(log(q))) # ???
        pfun <- switch(object$method, logistic = plogis, probit = pnorm, 
            cloglog = pgumbel, cauchit = pcauchy)
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
condquant <- function(x, dist='normal', sig=1)
{
  ## Purpose:   conditional quantiles and random numbers
  ## ----------------------------------------------------------------------
  ## Arguments:  
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 24 Oct 2008, 09:30
  fp <- switch(dist, normal='pnorm(x,0,sig)', logis='plogis(x)')
  fq <- switch(dist, normal=qnorm, logis=qlogis)
  if (c(dim(x),0,0)[2]!=2) stop('!condquant! x must have 2 columns')
  lp <- eval(parse(text=fp))
  lpp <- lp%*%rbind(c(0.5,0.25,0.75),c(0.5,0.75,0.25))
  lprand <- lp[,1]+(lp[,2]-lp[,1])*runif(nrow(lp))
  lr <- cbind(cbind(median=fq(lpp[,1]),lowq=fq(lpp[,2]),
              uppq=fq(lpp[,3]), random=fq(lprand))*sig,
              prob=abs(lp[,2]-lp[,1]))
#  dimnames(lr)[[1]] <- row.names(x)
  class(lr) <- "condquant"
  lr
}
## ===================================================================
residuals.survreg <- function(object, type='response', ...)
{
  ## Purpose:    conditional quantiles and random numbers for censored obs
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 24 Oct 2008, 10:16
  ly <- object$y
  lsig <- object$sigma
  ldist <- 'normal'
#  if (length(lsig)==0) lsig <- summary(object)$sigma
  if (length(lsig)==0) lsig <- summary(object)$scale
  if (length(lsig)==0) {
    warning('!residuals.survreg! no sigma found. Setting =1')
    lsig <- 1
  }
  lfit <- object$linear.predictors
  lres <- ly[,1]-lfit
  ## fill matrix with values adequate for non-censored obs
  lrr <- matrix(lres,length(lres),5)
  dimnames(lrr) <- list(row.names(ly),c('median','lowq','uppq','random','prob'))
  lrr[,'prob'] <- 0
  ## censoring
  lst <- ly[,2]  # status
  li <- lst!=1
  if (any(li)) {
    llim <- cbind(lres[li],c(Inf,NA,-Inf)[lst[li]+1]) # 
    lr <- condquant(llim, ldist, lsig)
    lrr[li,] <- lr
  }
  lrr <- cbind(lrr, fit=lfit)
##-   lna <- object$na.action
##-   if (length(lna)) naresid(lna, lrr) else  lrr
  class(lrr) <- 'condquant'
  naresid(object$na.action, lrr)
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
  ly <- object$model[,lyr]
  if (length(dim(ly))) {
    warning(":residuals.polr: returning simple deviance residuals for non-binary (grouped) data")
    return(residuals(object, type="deviance"))
  }
  #if (lpolr)
  ly <- as.numeric(ordered(ly))
  lfit <- fitted(object, type="link", na.action=list(na.action=NULL))     
  lthres <- c(-100, if (lpolr) object$zeta else 0, 100)
  llim <- cbind(lthres[ly],lthres[ly+1])-lfit
  lr <- cbind(condquant(llim,'logis'),fit=lfit,y=ly)
##-   lp <- cbind(plogis(),plogis(lthres[ly+1]-lfit))
##-   lpp <- lp%*%rbind(c(0.5,0.25,0.75),c(0.5,0.75,0.25))
##-   lprand <- lp[,1]+(lp[,2]-lp[,1])*runif(length(lfit))
##-   lr <- cbind(median=qlogis(lpp[,1]),lowq=qlogis(lpp[,2]),
##-               uppq=qlogis(lpp[,3]), random=qlogis(lprand), fit=lfit,y=ly)
  dimnames(lr)[[1]] <- row.names(object$model)
  class(lr) <- "condquant"
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
##-     function(object, newdata, type='linear.predictor', se.fit = FALSE) { 
##-     if (type=='linear.predictor') {
##-       lmod <- model.matrix(formula(object)[-2], newdata)
##-       if (inherits(object,'polr')) lmod <- lmod[,colnames(lmod)!='(Intercept)']
##-       rr <- lmod%*%object$coefficients[colnames(lmod)]
##-     } else rr <- predict(object, newdata, type=type, se.fit=se.fit)
##-     rr
##-   }
  if (length(data)==0) {
    data <- eval(object$call$data)
##-     datanm <- as.character(object$call)[3]
##-     if (is.na(datanm)) stop("no data found")
##-     data <- eval(parse(text=datanm))
  }
  lmeth <- object$call$method
  lnls <- length(lmeth)>0 && lmeth=="nls"
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
      xm[1,lj] <- if (is.factor(lv))
        levels(lv)[which.max(table(as.numeric(lv)))[1]]
      else   ## median(as.numeric(lv),na.rm=TRUE)
        sort(lv)[ceiling(length(lv)/2)] 
        ## median should be attained in some cases
    }
  }
  if (is.null(attr(terms(object), "predvars"))) { # from  model.frame
    lterms <- attr(lm(formula(object), data=data, method='model.frame'),'terms')
    attr(object$terms,'predvars') <- attr(lterms,'predvars')
  }
  lprm <- c(predict(object, newdata=xm)) # lf.
  lny <- length(lprm)
##  expand to matrix
  if (xfromdata) {
    lx <- ldata
  } else {
    lx <- ldata[1,,drop=FALSE][1:nxcomp,,drop=FALSE]
    row.names(lx) <- 1:nxcomp
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
        ld <- lxm[1:lnl,]
        ld[,lv] <- ldx
        lx[,lv] <- NA
        lx[1:lnl,lv] <- ldx
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
          seq(min(ldv,na.rm=TRUE),max(ldv,na.rm=TRUE),length=nxcomp)
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
  function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf, 
    interval = c("none", "confidence", "prediction"), level = 0.95, 
    type = c("response", "terms"), terms = NULL, na.action = na.pass, 
##-     pred.var = res.var/wgts, wgts = 1, ...) 
    pred.var = NULL, weights = 1, ...)
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
                stop("'weights' as formula should be one-sided")
            d <- if (missing(newdata) || is.null(newdata)) 
                model.frame(object)
            else newdata
            weights <- eval(weights[[2L]], d, environment(weights))
        }
    }
    type <- match.arg(type)
    if (se.fit || interval != "none") {
        res.var <- if (is.null(scale)) {
#            r <- object$residuals
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
modelTable <- 
  function(models, data=NULL, seq=NULL)
{
  ## Purpose:   collect several models into a table
  ## -------------------------------------------------------------------------
  ## Arguments:
  ##   models   character vector of names of model objects
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  9 Aug 2000, 11:50
##-   stopifnot(is.atomic(models),is.character(models))
  t.islist <- is.list(models)
  if ((!t.islist)&!(is.atomic(models)&is.character(models)))
    stop('!modelTable! Inadequate argument  model')
  t.nm <- length(models)
  t.ls <- vector("list",t.nm)
  t.sig <- rep(NA,t.nm)
  t.mnm <- if (t.islist) names(models) else models
  if (length(t.mnm)==0) t.mnm <- paste('model',1:t.nm,sep='.')
  names(t.ls) <- names(t.sig) <- t.mnm
  t.terms <- NULL
  t.dt <- as.character(substitute(data))
  for (li in 1:t.nm) {
    t.r <- if (t.islist) models[[li]] else get(models[li],envir=parent.frame())
    if (!(inherits(t.r,'lm')|inherits(t.r,'multinom')|inherits(t.r,'polr')))
      stop(paste('!modelTable! Model',li,'is not an adequate model'))
    t.dr <- drop1(t.r,test="F")
    t.tnm <- row.names(t.dr)
#    t.tnm <- substring(t.nm,3,max(nchar(t.nm)))
    t.dr$cf <- t.r$coef[match(t.tnm,names(t.r$coef))]
    t.ls[[li]] <- t.dr
    t.terms <- c(t.terms,t.tnm[-1])
    t.sig[li] <- c(summary(t.r)$sigma,NA)[1]
    t.dt <- c(t.dt, as.character(t.r$call[3]))
  }
  t.tr <- unique(t.terms)
  # --- standard deviations
  t.sd <- rep(1,length(t.tr))
  names(t.sd) <- t.tr
  t.dt <- unique(t.dt)
  if (length(t.dt)>1)
    warning(':modelTable: data arguments of model call are different: ',
            paste(t.dt,collapse=', '),'. Using the first one.')
  t.data <- get(t.dt[1])
  t.form <- as.formula(paste("~",paste(t.tr,collapse="+")))
  t.stnd <- all(all.vars(t.form)%in%names(t.data))
  if (t.stnd) {
    t.x <- model.matrix(t.form,t.data)
    t.vr <- apply(t.x,2,var)
    t.sd[t.tr] <- sqrt(t.vr)[t.tr]
  } else  {
    warning(':modelTable: cannot get the following terms from data:',
            paste(t.tr[!t.tr%in%names(t.data)], collapse=', '),
            '\n  Coefficients are not standardized.')
  }
  t.nt <- length(t.tr)
  t.cf <- matrix(NA,t.nt,t.nm)
  dimnames(t.cf) <- list(t.tr,t.mnm)
  t.pr <- t.cf
  for (li in 1:t.nm) {
    t.dr <- t.ls[[li]][-1,]
##-     t.t <- substring(row.names(t.dr),3,max(nchar(row.names(t.dr))))
    t.t <- row.names(t.dr)
    t.cf[t.t,li] <- t.dr[,"cf"]*t.sd[t.t]
    t.pr[t.t,li] <- t.dr[,"Pr(F)"]
  }
  if (length(seq)>0) {
    t.i <- match(seq, t.tr)
    t.t <- c(t.tr[t.i[!is.na(t.i)]],t.tr[!t.tr%in%seq])
    if (any(!t.tr%in%t.t|!t.t%in%t.tr)) warning("bug. not all terms. ask wst")
    t.cf <- t.cf[t.t,]
    t.pr <- t.pr[t.t,]
    t.sd <- t.sd[t.t]
  }
  attr(t.cf,'standardized') <- t.stnd
  if (all(is.na(t.sig))) t.sig <- NULL
  t.r <- list(coef=t.cf, p=t.pr, sd.terms=t.sd, sigma=t.sig)
  class(t.r) <- 'modelTable'
  t.r
}
## ==========================================================================
print.modelTable <-
  function(x, tex = FALSE, digits=3,
           stars = c("***","** ","*  ",":  ",".  "), ...)
{
  ## Purpose:   
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   x        a modelTable object
  ##   tex      if TRUE, the output will be suitable for pasting into
  ##            (la)tex source
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 23 Dec 2008, 10:09
  sep <- mc1 <- mc2 <- end <- headend <- tabstart <- tabend <- ''
  if(tex) {
    sep <- '&'
    mc1 <- '&\\mc{2}{c}{'
    mc2 <- '}'
    end <- '\\\\'
    headend <- '\\\\ \\hline'
    tabstart <- '  \\providecommand{\\begin{tabular}'
    tabend <- '  \\end{tabular}'
  }
  t.sd <- x$sd.terms
  t.st <- stars[cut(x$p,c(0,0.001,0.01,0.05,0.1,1))]
  t.st[is.na(t.st)] <- '   '
  dim(t.st) <- dim(x$p)
  if (!attr(x$coef,'standardized')) 
    warning(':: Coefficients are not standardized',call.=FALSE)
  t.cf <- x$coef
  if (length(x$sigma)) {
    t.cf <- rbind(t.cf, sigma=x$sigma)
    t.st <- rbind(t.st, '   ')
    t.sd <- c(t.sd, sigma=NA)
  }
  t.o <- paste(sep, format(t.cf, digits=digits), sep, t.st)
  t.out <-
    cbind(sep, format(t.sd), matrix(t.o, nrow=nrow(t.cf)))
  colnames(t.out) <- c(sep,'sd',paste(mc1,dimnames(x$p)[[2]],mc2))
  if (tex) {
    t.end <- c(rep(end,nrow(x$p)), headend, '')
    cat(tabstart, '{', rep('rl',ncol(t.out)),'}\n  ',
        colnames(t.out),headend,'\n')
    for (li in 1:nrow(t.out)) cat('  ',t.out[li,],t.end[li],'\n')
    cat(tabend,'\n')
  } else  print(t.out,quote=FALSE)
  invisible(t.out)
}
## ==========================================================================
## plotting functions
## ==========================================================================
plot.regr <-
function(x, data=NULL, markprop=NULL, lab=NULL, cex.lab=0.7, 
         mf = NULL, mfcol=FALSE, mar=c(3,3,2,1), mgp=c(2,0.7,0),
         oma = 2*(prod(mf)>1), cex=par("cex"), ask = NULL,
         multnrows = 0, multncols=0, 
         lty = c(1,2,5,3,4,6,1,1), lwd = c(1,1,2,1,1.5,1,1,1),
         colors = options("colors.ra")[[1]], pch=NULL, col=NULL, 
         main = NULL, cex.title = NULL, wsymbols=NULL, symbol.size=NULL, 
         smooth = TRUE, smooth.par=NA, smooth.iter=NA, 
         smooth.sim=19, nxsmooth=51, 
         plotselect = NULL, 
         weights = NULL, hat.cooklim = 1:2, reslim = TRUE, ylim=TRUE,
         ylimfac=3.0, ylimext=0.1,
         glm.restype = "deviance", condprobrange=c(0.05,0.8), 
	 sequence=NA, xplot = TRUE, x.se=FALSE, x.smooth = smooth,
         addcomp=FALSE, ...)
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
##  main       main title, defaults to the formula of  x
##             if it starts by : it will be appended to the formula
##  cex.title  character expansion for the main title
##  wsymbols   plot points by circles according to weights (of x)
##  symbol.size  ... to be used for these symbols
##  plotselect which plots should be shown?
##             0 : do not show,  1: show without smooth,  2: show with smooth
##             The default is
##             c( yfit=0, ta=3, tascale = NA, weights = NA, qq = NA, hat = 2)
##             modify this vector to change the selection and the sequence in
##             which the plots appear
##    yfit     response against fitted
##    ta       Tukey-Anscombe plot (residuals against fit)
##    tascale  absolute standardized residuals against fit
##             (smooth is calculated from sqrt(abs(stres)) )
##    weights  standardized residuals against weights to be given in
##             argument weights or as x$weights
##    qq       normal plot of standardized residuals 
##    hat      residuals against diagonal of hat matrix
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
## family
  lfam <- x$familyname
  if (length(lfam)==0) lfam <- x$family$family
  if (is.null(lfam) || lfam=="" || is.na(lfam)) lfam <- "gaussian"
  lfgauss <- lfam == "gaussian"
  lglm <- !lfgauss
#  lpolr <- inherits(x,"polr")
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
              hat = 2, resmatrix = 1, qqmult = 3)
  if (length(plotselect)>0) {
    lplnm <- names(plotselect)
    if (length(lplnm)==0) {
      if (length(plotselect)==length(lplsel))
        lplnm <- names(lplsel)
        else {
        warning(":plot.regr: Inadequate argument plotselect")
        plotselect <- NULL }
    } else {
      lplsel[] <- if ("default"%in%lplnm)
        pmin(lplsel,plotselect["default"]) else 0
      lina <- is.na(match(lplnm,names(lplsel)))
      if (any(lina)) {
        warning(paste(":plot.regr: Inadequate elements plotselect:",
                      paste(names(plotselect)[lina]),collapse=", "))
        lplnm <- lplnm[!lina] }
      lplsel[lplnm] <- plotselect[lplnm]
    }
  }
  lseq <- NA
## -----------------------------------
## prepare objects needed for plotting
  lnaaction <- x$na.action
  x$na.action <- NULL
  rtype <- "response" 
  if (lglm) rtype <- glm.restype
##-   if (lsurv) rtype <- "condquant"
##-   lcondq <- !is.na(pmatch(rtype,"condquant"))
##-   if(!lcondq)
##-     lres <- residuals(x, type=rtype) else {
##-       lres <- if (lpolr) residuals.polr(x) else residuals.survreg(x)
##-       lcondq <- NCOL(lres)>1
##-     }
  lres <- residuals(x, type=rtype)
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
  lnna <- !is.na(lres[,1])
  ln <- nrow(lres)
  lresa <- lres[lnna,,drop=FALSE]
  class(lresa) <- class(lres)
  lna <- nrow(lresa)
  lrname <- paste("res(", lyname, ")")
  ## plot range
  if (is.logical(reslim)&&!is.na(reslim)) reslim <- 
    if (reslim) {
      if(lcondq) robrange(c(lresa)) else apply(lresa,2,robrange)
    } else  NULL
  if (!is.null(reslim))
    if (any(dim(cbind(reslim))!=c(2,lmres))) {
      warning(":plot.regr: unsuitable argument  reslim ")
      reslim <- NULL
    }
##-   ldfres <- x$df.residual
  ldfres <- x$df
  if(length(ldfres)>1) ldfres <- ldfres[2] 
  ## fit
  lf <- x$linear.predictors  # predict(x)  does not work for nls
  lfname <- "Linear Predictor"
  if(is.null(lf)) {
    lf <- fitted(x)
    lfname <- "Fitted Values"
  }
  lf <- cbind(lf)[lnna,,drop=FALSE]
  lsigma <- x$sigma
  if (length(lsigma)==0) lsigma <- if (lfcount) 0 else
    sqrt(apply(lresa^2,2,sum)/ldfres)
#  lsigma <- x$sigma
## weights
  if (lwgt) {
    lwgts <- lweights[lnna]
    if (length(lwgts)!=lna) {
      warning(c(":plot.regr: there is something wrong with the weights.",
                "They are not used."))
      lwgt <- lwsymbols <- FALSE }
  } else lwgts <- NULL
## hat
  lhat <- x$h[lnna]
  if(length(lhat)==0&&length(x$qr)>0) lhat <- hat(x$qr) # [lnna]
##-   if (lwgt) lhat <- lhat/lwgts use almost-unweighted h for plotting
## standardized residuals
  lstres <- x$stres
  lstres <- if (length(lstres)>0)
    cbind(lstres)[lnna,,drop=FALSE] else {
      if (all(lsigma>0)) {
        if (length(lhat)==0||any(is.na(lhat))) sweep(lresa,2,lsigma,"/") else
        lresa/outer(ifelse(lhat>=1,1,sqrt(pmax(0,1-lhat))),lsigma)
      } else  lresa
    }
#  lstres <- lstres[lnna,,drop=FALSE]
##-     if (lwgt) lstres <- lstres*sqrt(lwgts) ## !!! check
  lrabs <- abs(lstres)
  lstrname <- paste("st", lrname, sep = ".")
  if (lmult) {
    lresmd <- x$resmd
    if (is.null(lresmd)) lresmd <- mahalanobis(lresa,0,var(lresa))
    lresmd <- lresmd[lnna]
  }
## color vector
  if (length(col)) col <- rep(col,ln)[lnna]
## smooth
  if (is.logical(smooth)) smooth <- if (smooth) 
    function(x,y,weights,par,iter)
      fitted(loess(y~x, weights=weights, span=par, iter=iter,
                   na.action=na.exclude))
    else NULL
  lsmpar <- if (is.na(smooth.par)) 3*lna^log10(1/2)*(1+lglm) else
               smooth.par# 2
  if (length(smooth.iter)==0||is.na(smooth.iter))
    smooth.iter <- if (lfgauss) 3 else 0
  lnsims <- smooth.sim
  if (!lfgauss) lnsims <- 0
  if (length(lnsims)==0) lnsims <- 0
  lnsims <- if (is.logical(lnsims)&&lnsims) 19 else as.numeric(lnsims)
  if (!lfgauss) lnsims <- 0
  if (lmult) lnsims <- 0 # not yet programmed for mlm
  lsimres <- NULL
  if (lnsims>0) {
    lsimr <- simresiduals(x, lnsims)
    lsimres <- lsimr$simres
    lsimstres <- lsimr$simstres
    if (length(lsimres)==0) lnsims <- 0
  }
## labels
  ## priorities:  lab , pch , row.names
  lpch <- pch
  if (length(lpch)==ln) {
    if (length(lab)==0) lab <- lpch
    lpch <- NULL }
  if (length(lpch)==0) lpch <- ifelse(ln>200,".",3)
  lrown <- names(lres)
    if (length(lrown)==0) lrown <- as.character(1:ln)
  if (length(lab)>1&length(lab)!=ln) {
    warning(":plot.regr: argument  lab  has unsuitable length")
    lab <- NULL
  } 
  llabels <- if (length(lab)==0) lrown else
  rep(if (is.character(lab)) lab else as.numeric(lab), length=ln)[lnna]
                                        # factors: as.numeric
  if (length(markprop)==0||is.na(markprop))
    markprop <- if (length(lab)>1) 1 else ceiling(sqrt(lna)/2)/lna
  if (markprop==0) {
    llab <- rep(lpch, length=lna)
    llabna <- NA
  } else  {
    llab <- llabels[lnna]
    if (markprop<1) {
      li <- if (lmult)  order(lresmd)[1:(lna*(1-markprop))]  else 
      order(abs(lstres[,1]))[1:(lna*(1-markprop))]  # [,1] for condq
      llab[li] <- if (is.numeric(llab)) {
        if(is.numeric(lpch)) lpch[1] else 1
      } else  {
        if (is.numeric(lpch)) "+" else lpch[1]
      }
      llabna <- llab
      llabna[li] <- NA
    }
  }
  if (!lwgt) llabna <- llab
## weights to be used for symbols
## plot symbols
  if (length(symbol.size)==0||is.na(symbol.size)) symbol.size <- 3/ln^0.3
  if (lwsymbols) {
    liwgt <- is.na(llabna)
    if (any(liwgt)) {
      lwgts <- lwgts/mean(lwgts,na.rm=TRUE)
      lsyinches <- 0.02*symbol.size*par("pin")[1] # *max(lwgts,na.rm=TRUE)
      llab[liwgt] <- ifelse(is.character(llab),"",0)
    } else lwgt <- FALSE
  }
  ltxt <- is.character(llab) # &&any(nchar(llab)>1)
  lpty <- ifelse(ltxt,"n","p")
## -----------------------------------
## prepare graphical elements
  lty <- rep(c(lty,1:6),length=6)
  lwd <- rep(c(lwd,1),length=6)
  if (length(ask)==0) ask <- !last(c("",names(dev.list())))=="postscript"
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
    c("black","gray","blue","cyan","red","magenta","darkgreen",
      "green","lightgreen")  else
    rep(colors,length=9)
  lftext <- paste(as.character(lform)[c(2,1,3)],collapse="")
  if (length(main)==0) main <- lftext
  if (is.logical(main)) main <- if (main) lftext else ""
  if (is.character(main)&&substring(main,1,1)==":")
    main <- paste(lftext,substring(main,2,30))
  main <- as.character(main)
  if (length(main))  tit(main) <- tit(x)
  if (is.null(cex.title)) cex.title <- max(0.5, min(1.2,
      par("mfg")[4]*par("pin")[1]/(par("cin")[1]*nchar(main))))
  outer.margin <- par("oma")[3]>0
## --------------------------------------------------------------------------
## start plots
  if (!lmult) lplsel[c("resmatrix","qqmult")] <- 0
  lplsel <- lplsel[is.na(lplsel)|lplsel>0]
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
      if (length(ly)==0) ly <- lf + lresa
      lsimr <- c(lf)+lsimres
      i.plotlws(lf,ly, lfname,lyname,main,outer.margin,cex.title,
                colors,lty,lwd,lpch, col,
                lpty, llabna,cex.lab,ltxt, lwsymbols,lwgt,lwgts,liwgt,lsyinches,
                lpllevel>1,smooth,lsmpar,smooth.iter, ylim=ylim,
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
      i.plotlws(lf,lresa, lfname,lrname,main,outer.margin,cex.title,
                colors,lty,lwd,lpch, col, lpty, llabna,cex.lab,ltxt, 
                lwsymbols,lwgt,lwgts,liwgt,lsyinches,
                lpllevel>1,smooth,lsmpar,smooth.iter, ylim=reslim,
                reflinex=lrx,refliney=lry,
                lnsims=lnsims, simres=lsimres,
                condprobrange=condprobrange)
    }
  }
## --- 
  if(lpls=="tascale")  {
    if (is.na(lpllevel)) lpllevel <- 3*(!lfcount)
    if (lpllevel>0) 
      i.plotlws(lf,lrabs, lfname, paste("|",lstrname,"|"),
        main, outer.margin, cex.title,
        colors,lty,lwd,lpch, col, lpty, llabna,cex.lab,ltxt,
        FALSE,FALSE, # lplwgt=F,wgt=F,
        lwgts,liwgt,lsyinches, lpllevel>1, smooth,lsmpar,smooth.iter, 
        smooth.power=0.5,
        ylim=c(0,1.05*max(lrabs,na.rm=TRUE)),yaxs="i",
                condprobrange=c(0.01,0))
    if (lpllevel>2&lnsims>0&&length(lsimstres)>0) {
      lio <- order(lf)
      for (lr in 1:lnsims) {
        ys <- smooth(lf, sqrt(abs(lsimstres[,lr])),
                     weights=lwgts, par=lsmpar, iter=smooth.iter)
        lines(lf[lio], ys[lio]^2, lty=lty[4], lwd=lwd[4], col=colors[4])
      }
    }
  }
## --- plot abs. res vs. weights
  if(lpls=="weights") {  
  if (is.na(lpllevel)) lpllevel <- 3*(lfgauss&lwgt)
#  lplweights <- (is.logical(weights) && weights) | (length(weights)>1)
  if (lpllevel>0) { # plot on weights
    lweights <- if (length(weights)>1) weights[lnna] else lwgts
    if (length(lweights)!=length(lrabs)) 
      warning(":plot.regr: no weights foound. cannot plot absres on weights")
    else {
      i.plotlws(lweights,lrabs, "weight (log)",
        paste("|",lstrname,"|"),main,outer.margin, cex.title,
        colors,lty,lwd,lpch, col, lpty, llab,cex.lab,ltxt, FALSE,FALSE, # lplwgt=F,wgt=F,
        lwgts,liwgt,lsyinches,
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
      lxy <- qqnorm(lstres[,lj], ylab = lstrname[lj], main="", type="n")
      abline(0,1,lty = lty[2], col=colors[2])
      if (lpllevel>2&lnsims>0) {
        lxx <- qnorm(ppoints(lna))
        lsimstr <- simresiduals(x, nrep=lnsims, resgen=rnorm)$simstres
        for (lr in 1:lnsims) {
          lines(lxx,sort(lsimstr[,lr]),lty=lty[4], lwd=lwd[4], col=colors[4])
        }
#        lines(lxx, sort(lstres[,lj]), col=colors[1])
      }
      li <- order(lxy$x)
      lxx <- lxy$x[li]
      lyy <- lxy$y[li]
      lines(lxx,lyy, col=colors[1])
      if (is.character(llab)) text(lxx,lyy,llab[li], col=colors[1]) else
        points(lxx,lyy, pch=lpch, col=colors[1])
##-     lquart <- quantile(lstresa,c(0.25,0.75))
##-     abline(0, diff(lquart)/(2*qnorm(0.75)), lty = lty[2], col=colors[2])
      i.main(main, outer.margin=outer.margin)
      stamp(sure=FALSE) }
    
  }
## --- leverage plot. If weight are present, use "almost unweighted h"
  if(lpls=="hat") 
  if ((!is.na(lpllevel))&&lpllevel>0) {
    if (diff(range(lhat,na.rm=TRUE))<0.001)
      warning(":plot.regr: all hat elements equal, no leverage plot")
    else {
      llabh <- llab
      if (markprop>0 & markprop<1) {
        li <- order(lhat, decreasing=TRUE)[1:(lna*markprop/2)]
        llabh[li] <- llabels[lnna][li]
      }
      lhattit <- paste("hat diagonal", if(lwgt) "(unweighted)")
      for (lj in 1:lmres) {
        i.plotlws(lhat, lstres[,lj], lhattit, lstrname[lj],
              main,outer.margin,cex.title,
              colors,lty,lwd,lpch, col, lpty, llabh,cex.lab,ltxt,
              lwsymbols,lwgt, lwgts,liwgt,lsyinches,
              FALSE, ylim=reslim,
                condprobrange=condprobrange) 
  ##  line with constant Cook distance
      if (length(hat.cooklim)>0) {
        llx <- seq(min(c(hat.cooklim,4),na.rm=TRUE)^2*(1-ldfres/lna)/6,
                   max(lhat,na.rm=TRUE),length=50)
        llr <- outer(sqrt((1-llx)^2*(lna-ldfres)/(llx*ldfres)),
                     c(hat.cooklim,-hat.cooklim))
        matlines(llx, llr, lty=lty[2], lwd=lwd[2], col=colors[2])
      }
      }
    }
  }
  if(lpls=="resmatrix") { ## residual matrix for multivariate regr
    if (lmult) plmatrix(lresa, pch=lab, main=main)
  }
  if(lpls=="qqmult") { ## qq plot of Mahalanobis lenghts for multivariate regr
  if ((!is.na(lpllevel))&&lpllevel>0) {
    lxx <- sqrt(qchisq(ppoints(lresmd),ncol(lresa)))
    lor <- order(lresmd)
    lyy <- sqrt(lresmd[lor])
    lop <- par(mfrow=c(1,1))
    plot(lxx,lyy, xlab="sqrt(Chisq.quantiles)",type="n",
             ylab = "Mahal.oulyingness", main="", col=colors[1])
    if (ltxt) text(lxx,lyy,llab[lor]) else points(lxx,lyy,pch=llab[lor])
    axis(1)
    axis(2)
    abline(0,1,lty = lty[2], col=colors[2])
    i.main(main)
    stamp(sure=FALSE)
    par(lop)
  }}
## --- 
##-   if (lpls=="seq")
##-   if ((!is.na(lpllevel))&&lpllevel>0) {
##-     lsq <- (1:length(lres))[lnna]
##-     i.plotlws(lsq,lresa, "sequence",lrname,main,outer.margin,cex.title,
##-               colors,lty,lpch, col,
##-             lpty, llabna,cex.lab,ltxt, lwsymbols,lwgt,lwgts,liwgt,lsyinches,
##-             lpllevel>1,smooth,lsmpar,smooth.iter, ylim=reslim)
##-   }
##-     lseq <- lpllevel
} ## end lplsel
## ----------------------------------------------------------------
## plot residuals vs. explanatory variables by calling plresx
  lseq <- sequence
  if (length(lseq)==0) lseq <- NA
  if (is.logical(xplot)) xplot <- if(xplot) lform else NULL
  if(lxpl <- length(xplot) > 0) {
    if (is.character(xplot))
      xplot <- formula(paste("y ~",paste(xplot,collapse="+")))
    xplot <- update(lform, xplot)
    xvars <- all.vars(xplot[-2])
  } else xvars <- FALSE
#  llabs <- rep(NA,ln)
#  llabs[lnna] <- llabna
  lylim <- if (addcomp) ylim else reslim
  if (lxpl | is.na(lseq)|lseq) {
    plresx(x, resid=lres, partial.resid=TRUE, glm.restype = glm.restype,
           lab=llab, cex.lab=cex.lab,
           vars=xvars, sequence=lseq, se=x.se, addcomp = addcomp,
           smooth=x.smooth, smooth.par=lsmpar, 
           smooth.sim=lsimres, lty=lty, lwd=lwd, colors = colors, col=col,
           main=main, cex.title=cex.title, ylim = lylim, 
           wsymbols=lwsymbols, symbol.size=symbol.size, pch=lpch, ...)
  }
## --- end
  par(loldpar)
  "plot.regr done"
}
## ==========================================================================
i.plotlws <- function(x,y, xlab="",ylab="",main="", outer.margin=FALSE,
  cex.title=1.2, colors=1:8, lty=1:6, lwd=1, pch=1, col=1,
  pty="p", lab="+", cex.lab=1,
  txt=FALSE,  plwgt=FALSE, wgt=1, weights=NULL, iwgt=NULL, syinches=0.1, 
  do.smooth=TRUE, smooth=i.smooth, smooth.par=NULL, smooth.iter=10,
  smooth.power=1, ylim=NULL, ylimfac=3.0, ylimext=0.1,
  reflinex=NULL, refliney=NULL, reflineyw=NULL, lnsims=0, simres=NULL,
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
  if (lcq) lnsims <- 0
  lnc <- if (lcq) 1 else ncol(ly)
  lnr <- nrow(lx)
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
  reflineyw <- cbind(reflineyw)
  ldjrflx <- ncol(reflinex)==lnc
  ldjrfly <- ncol(refliney)==lnc
##-   lrfyw <- FALSE
##-   if (!is.null(reflineyw)) {
##-     if (nrow(reflineyw)<2)
##-       warning(":plot.regr/i.plotlws: reference line width not suitable") else
##-     lrfyw <- TRUE
##-   }
  lrfyw <- length(reflineyw)>0
## do plot(s)
  ljx <- ljrflx <- ljrfly <- 1
  for (lj in 1:lnc) {
    lxj <- lx[,ljx]
    if (lcq) {
      lyj <- y
      if (llimy) {
        lylimj <- lylim[,lj]
        lyplj <- lypl <- plcoord(y, range=lylimj, limext=ylimext)
      } else lyplj <- lypl <- y
    } else {
      lyj <- ly[,lj]
      if (llimy) {
        lylimj <- lylim[,lj]
        lyplj <- plcoord(lyj, range=lylimj, limext=ylimext)
      } else lyplj <- lyj
    }
    if (new) {
      lyrgj <- range(lyplj,na.rm=TRUE)
      plot(range(lxj,na.rm=TRUE), lyrgj, xlab = lxlab[lj], ylab = lylab[lj],
           type="n", bty="n", axes=FALSE, ...)
        if (1%in%axes) axis(1)
      if (llimy & length(attr(lyplj,"nmod"))) {
        lusr <- par("usr")
        box(lty=3)
        lines(lusr[c(1,1,2,2,1)],
              c(max(lylimj[1],lusr[3]),min(lylimj[2],lusr[4]))[c(1,2,2,1,1)],
              xpd=TRUE)
        if (2%in%axes) {
          lat <- pretty(lylimj, min.n=5)
          lat <- lat[lat>=lylimj[1]&lat<=lylimj[2]]
          axis(2,at=lat)
        }
      } else {
        box()
        if (2%in%axes) axis(2)
      }
      abline(0, 0, lty = lty[2], col=colors[2])
    }
    lusr <- par("usr")
    ## conditional quantiles
    if(lcq) {
      li <- lypl[,"prob"]>=condprobrange[1] & lypl[,"prob"]<=condprobrange[2]
      if (any(li))
      segments(lxj[li],lypl[li,"lowq"],lxj[li],lypl[li,"uppq"],col=colors[9])
      ldx <- diff(lusr[1:2])*0.02
      segments(lxj-ldx,lypl[,"median"],lxj+ldx,lypl[,"median"],col=colors[8])
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
##-     }
##-     if (lrf) {
      if (llimy) {
        if (length(lrfxj)==2) {
          lrfxj <- seq(lusr[1],lusr[2],length=50)
          lrfyj <- reflinex[,ljrflx]+refliney[,ljrfly]*lrfxj
          if (lrfyw) {
            lrfylj <- lrfyj-reflineyw[,ljrfly]*lrfxj
            lrfyuj <- lrfyj+reflineyw[,ljrfly]*lrfxj
          }
        }
##-         lrfyj <- plcoord(lrfyj, range=lylimj, limext=ylimext)
        lrfyj[lrfyj<lylimj[1]|lrfyj>lylimj[2]] <- NA
        if (lrfyw) {
          lrfylj[lrfylj<lylimj[1]|lrfylj>lylimj[2]] <- NA
          lrfyuj[lrfyuj<lylimj[1]|lrfyuj>lylimj[2]] <- NA
        }
      }
##-     }
##-     if (lrf) {
      lines(lrfxj,lrfyj,lty=lty[5],col=colors[5],lwd=lwd[5])
      if (lrfyw) {
        lines(lrfxj,lrfylj,lty=lty[6],col=colors[6],lwd=lwd[6])
        lines(lrfxj,lrfyuj,lty=lty[6],col=colors[6],lwd=lwd[6])
      }
    }
    if(do.smooth) {
      if (lcq) lyj <- ly[,"median"]
      lna <- (is.na(lxj)|is.na(lyj))
      lxj[lna] <- NA
      lio <- order(lxj)[1:sum(!lna)]
      lxjo <- lxj[lio] # sorted without NA
      lwgo <- weights[lio]
      if (lnsims>0) {
        simres <- simres[lio,]
        for (lr in 1:lnsims) {
          lsms <- smooth(lxjo, simres[,lr]^smooth.power, weights=lwgo,
                       par=smooth.par, iter=smooth.iter)^(1/smooth.power)
          if (llimy)  lsms[lsms<lylimj[1]|lsms>lylimj[2]] <- NA
          lines(lxjo, lsms,
                lty=lty[4], col=colors[4],lwd=lwd[4])
        }
      }
      lsm <- smooth(lxjo,lyj[lio]^smooth.power, weights=lwgo,
                    par=smooth.par, iter=smooth.iter)^(1/smooth.power)
      if (llimy)  lsm[lsm<lylimj[1]|lsm>lylimj[2]] <- NA
      lines(lxjo, lsm,
            lty = lty[3], col=colors[3], lwd=lwd[3])
    }
##- points
    lpclr <- if (length(col)) rep(col,length=lnr) else colors[1]
    if (lcq) {
      lyplj <- lypl[,"random"]
      lpclr <- if (length(col)) rep(col,length=lnr) else colors[7]
    }
    if (length(lyj)) {
    if (txt) {
      text(lxj,lyplj, lab, cex=cex.lab, col=lpclr)
      if (plwgt) {
        if (is.null(iwgt)) iwgt <- rep(iwgt,length(x))
        symbols(lxj[iwgt], lyplj[iwgt], circles=sqrt(weights[iwgt]),
                     inches=syinches, fg=lpclr, add=TRUE)
      }
      else  if (any(li <- is.na(lab)|lab==""))
        points(lxj[li], lyplj[li], pch=pch, cex=0.7*cex.lab, col=lpclr)
    } else  points(lxj, lyplj, pch=lab, cex=0.7*cex.lab, col=lpclr)
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
i.smooth <-function(x,y,weights,par=3*length(x)^log10(1/2),iter=50)
  fitted(loess(y~x, weights=weights, span=par, iter=iter,
               na.action=na.exclude))
## ==========================================================================
simresiduals <- function(object, nrep, resgen=NULL)
{
  ## Purpose:   simulate residuals according to regression model
  ##            by permuting residuals of actual model
  ## ----------------------------------------------------------------------
  ## Arguments:  resgen: how are residuals generated?
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 10 Aug 2008, 07:55
  lcall <- object$call
  lerror <- NULL
  ldata <- try(eval(lcall$data))
  if (class(ldata)=="try-error") lerror <- "data not found" else {    
    lres <- object$stres*object$sigma
    if ((length(lres)==0)||all(!is.finite(lres))) lres <- object$resid
    if (length(lres)==0||all(lres==lres[1]))
      lerror <- "no (distinct) residuals found"  else {
      if (nrow(ldata)!=length(lres)) {
        li <- match(names(lres),row.names(ldata))
        ldata <- ldata[li,] 
        if (any(is.na(li)))
          lerror <- "data not suitable" else {
        ##!!! weights
          lina <- is.na(lres)
          lres <- lres[!lina]
          ldata <- ldata[!lina,]
          if (nrow(ldata)<=2) lerror <- "<=2 non-missing residuals"
      }}}}
  if (length(lerror)) {
    warning(paste(":plot.regr/simresiduals:",lerror, "-> no simulated smooths"))
    return(list(simres=NULL, simstres=NULL))
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
  lform <- update(lform, paste(lynm,"~.")) ## needed for transformed y
  environment(lform) <- environment()
  lcall$formula <- lform
  lcall$model <- NULL
  lcall$termtable <- NULL
  lnrow <- nrow(ldata)
  lsimstres <- lsimres <- matrix(NA,lnrow,nrep)
  for (lr in 1:nrep) {
    ldata[,lynm] <- if (lrgen) resgen(lnrow)*lsig else sample(lres)
    lrs <- eval(lcall) ## update(x, formula=lfo, data=ldata)
    lsimres[,lr] <- lrs$resid
    lsimstres[,lr] <- lrs$stres
  }
  list(simres=lsimres, simstres=lsimstres)  
}
## ==========================================================================
plresx <-
  function (x, data = NULL, resid=NULL, partial.resid = TRUE,
  glm.restype = "deviance", weights = NULL, lab = NULL, cex.lab = 1,
  vars = NULL, sequence=FALSE, se = FALSE,
  addcomp = FALSE, rug = FALSE, jitter=NULL,
  smooth = TRUE, smooth.par=NA, smooth.iter=NA, smooth.sim=19, 
  nxsmooth=51, lty = c(1,2,5,3,6,4,1,1), lwd=c(1,1,2,1,1.5,1,1,1),
  colors = options("colors.ra")[[1]], pch=NULL, col=NULL, 
  xlabs = NULL, ylabs = NULL, main = NULL, cex.title = NULL, 
  ylim=TRUE, ylimfac=3.0, ylimext=0.1,
  cex = par("cex"), wsymbols=NULL, symbol.size=NULL, condprobrange=c(0.05,0.8),
  ask = NULL, multnrows = 0, multncols = 0, ...) 
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
## family
  lfam <- x$familyname
  if (is.null(lfam)) lfam <- x$family$family
  if (is.null(lfam) || lfam=="" || is.na(lfam)) lfam <- "gaussian"
  lfgauss <- lfam == "gaussian"
  lglm <- !lfgauss
  lpolr <- inherits(x,"polr")
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
  if (lmult) {
    if (is.null(lcn <- colnames(lres))) lcn <- 1:ncol(lres)
    colnames(lres) <- paste("res",lcn,sep=".")
  }
  lnna <- !is.na(lres[,1])
  ln <- nrow(lres)
  lresa <- lres[lnna,,drop=FALSE]
  class(lresa) <- class(lres)
  lna <- nrow(lresa)
  llab <- if (length(lab)==0) rep(NA,length=lna)  else
    rep(lab,length=length(lres))[lnna]
  ltxt <- is.character(llab) # &&any(nchar(llab)>1)
  lpty <- ifelse(ltxt,"n","p")
  liwgt <- is.na(llab)
  if (is.character(llab)) {
    liwgt <- liwgt| llab=="NA"|llab==""
    lpch <- if (length(pch)) pch else ifelse(lna>200,".","+")
  } else lpch <- if (length(pch)) pch else 3
## weights
  lwgt <- if (length(weights)==1) as.logical(weights) else NA
  lwgts <- if (length(weights)<=1) x$weights else weights
  if (is.na(lwgt)) lwgt <- length(lwgts)>1
  if (lwgt&& length(lwgts)!=length(lres)) {
      warning(":plresx: no suitable weights found")
      lwgt <- FALSE
  }
  lwsymbols <- lwgt&any(liwgt)
## color vector
  if (length(col)) col <- rep(col,ln)[lnna]
## data
  ldata <- if (length(data)==0) eval(x$call$data) else data
  if (nrow(ldata)!=nrow(lres)) {
    li <- match(row.names(lres),row.names(ldata))
    if (length(li)==0||any(is.na(li))) 
      stop("!plresx! cannot match residuals and x's")
    else {
      ldata <- ldata[li,]
      if (lwgt) lwgts <- lwgts[li]
    }
  }
## weight symbols
  if (lwgt) {
    lwgts <- lwgts[lnna]/mean(lwgts,na.rm=TRUE)
    if (length(symbol.size)==0||is.na(symbol.size)) symbol.size <- 3/lna^0.3
    lsyinches <- 0.02*symbol.size*par("pin")[1] # *max(lwgts,na.rm=TRUE)
    if (lwsymbols) llab[liwgt] <- ifelse(is.character(llab),"",0)
  } else  {
    lwgts <- NULL
    llab[liwgt] <- lpch
  }
  lipts <- !(lwsymbols&liwgt) ## points to be shown by  lab
## Prepare vars
  lvmod <- all.vars(formula(x)[[3]])
  if (match(".",lvmod,nomatch=0)>0) lvmod <- names(x$model)[-1]
  if ((is.logical(vars)&&!vars)) vars <- NULL else {
    if (length(vars)==0) {
      vars <- lvmod
      if (length(data)==0) data <- x$model
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
  mv <- match(vars,dvars,nomatch=0)
  if (any(mv==0)) {
    warning (paste(":plresx: Variable(s)",
      paste(vars[match(vars,dvars,nomatch=0)==0],collapse=", "),
                "not found"))
    vars <- vars[mv!=0]
  }
## convert character to factor
  for (lvn in dvars) {
    lv <- ldata[[lvn]]
    if (is.character(lv)) ldata[[lvn]] <- factor(lv)
  }
  if (length(sequence)==0||is.na(sequence)) {
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
  }
  if (sequence>0) {
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
## smooth
  if (is.logical(smooth)) smooth <- if (smooth) { 
    function(x,y,weights,par,iter)
      fitted(loess(y~x, weights=weights, span=par, iter=iter,
                   na.action=na.exclude))
    } else NULL
  ldosm <- length(smooth)>0
  lsmpar <- if (is.na(smooth.par)) 3*lna^log10(1/2)*(1+lglm) else
               smooth.par# 2
  if (length(smooth.iter)==0||is.na(smooth.iter))
    smooth.iter <- if (lfgauss) 3 else 0
  lnsims <- smooth.sim
  if (length(lnsims)==0) lnsims <- 0
  lnsims <- if (is.logical(lnsims)&&lnsims) 19 else as.numeric(lnsims)
  if (!lfgauss|lcondq) lnsims <- 0
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
  if (length(jitter)==0) jitter <- 0.3*(1-10^(-0.01*pmax(0,lna-10)))
  if (is.na(smooth.par)) smooth.par <- 3*sum(lnna)^log10(1/2)*(1+lglm) # 2
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
      if(lcondq) robrange(c(lresa[,1:3])) else apply(lresa,2,robrange)
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
  ldt <- ldata[lnna,,drop=FALSE]
##-   if (sum(terminmodel)>0) {
##-     comps <- fitcomp(x, ldt, vars=tin, xfromdata=addcomp, se=se)
##-     if (addcomp) {
##-       if (length(ldt)==0&&(sum(lnna)!=nrow(ldt))) # length(res)
##-         stop("plresx: BUG: number of residuals != nrow(comps)")
##-       xcomp <- ldt[,tin]
##-       lcmp <- comps$comp # [lnna,]
##-     } else {
##-       xcomp <- comps$x
##-       lcmp <- -comps$comp
##-     }
##-   }
##-   environment(x$call$formula) <- environment()
##-   x$call$data <- as.name("ldt")
  if (sum(terminmodel)>0) {
    lcmp <- try(fitcomp(x, ldt, xfromdata=FALSE, se=se))
    if (class(lcmp)=='try-error') {
      warning(':plresx: fitcomp did not work. no reference lines')
      terminmodel[] <- FALSE
    } else {
      lcompx <- lcmp$x
      lcompy <- if (addcomp) lcmp$comp else -lcmp$comp
      lcompse <- lcmp$se
      if (addcomp) {
        if (length(ldt)==0&&(sum(lnna)!=nrow(ldt))) # length(res)
          stop("!plresx! BUG: number of residuals != nrow(comps)")
        lcompdt <- fitcomp(x, ldt, vars=vars[terminmodel], xfromdata=TRUE)$comp
      }
    }
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
                lpty, llab,cex.lab,ltxt, lwsymbols,lwgt,lwgts,liwgt,lsyinches,
                ldosm&lcnt,smooth,lsmpar,smooth.iter, 
##-                 ylim=lylim, 
                lnsims=lnsims, simres=lsimres, new=FALSE,
                reflinex=lcmpx, refliney=lcmpy)
    }
    plmatrix(lres,ldata[,vars,drop=FALSE], panel=lpanel, pch=llab, range=lylim,
             reference=FALSE, 
             nrows=multnrows, ncols=multncols, main=main) # clrsmooth=colors[3]
    return()
  }
## --- loop ---
  for (lj in 1:nvars) {
    lv <- vars[lj]
    lci <- if (terminmodel[lj]) lcompy[, lv] else 0 ## !!!
    rr <- lresa
    if (partial.resid) 
      if (addcomp&terminmodel[lj]) rr <- rr+lcompdt[, lv]
#!    class(rr) <- class(lresa)
## - plot range
##-     if (llimy) rr <- plcoord(rr, range=lylim, limext=ylimext)  else
##-       if (llimyrob) rr <- plcoord(rr, limfac=ylimfac, limext=ylimext)
##-     ylims <- if (partial.resid) range(rr, na.rm = TRUE) else 
##-                 range(lci, na.rm = TRUE)
##-     if (rug)  ylims[1] <- ylims[1] - 0.07 * diff(ylims)
## ---
    if (is.fac[lv]) { # ---
## factors
      ff <- factor(ldt[, lv])
      ll <- levels(ff)
      lnl <- length(ll)
      xx <- as.numeric(ff)+runif(lna,-jitter,jitter)
      xlims <- c(0.5,lnl+0.5) 
      i.plotlws(xx, rr, xlab = xlabs[lv], ylab = ylabs[lv],
        "", outer.margin, cex.title,
        colors, lty, lwd, lpch, col, lpty, llab, cex.lab, ltxt,
        lwsymbols, lwgt, lwgts, liwgt, lsyinches, 
        FALSE, smooth, smooth.par, smooth.iter, 1,
        lylim, ylimfac, ylimext,
        reflinex=NULL, lnsims=0, axes=2,
                condprobrange=condprobrange)
      axis(1,labels=ll,at=1:lnl)
      axis(2)
      box(lty=3)
      ## - 
      if (terminmodel[lj]) {
        lx <- seq(along = ll)
##-         ww <- if (addcomp) match(ll,as.character(ff)) else 1:lnl
        lcil <- lci[1:lnl]
        if (llimy) lcil[lcil<lylim[1]|lcil>lylim[2]] <- NA
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
      if (terminmodel[lj]) {
        lrefx <- lcompx[,lv]
        if (se) lrefyw <- qnt*lcompse[,lv] 
      }
      i.plotlws(as.numeric(ldt[, lv]), rr, xlab = xlabs[lv], ylab = ylabs[lv],
        "", outer.margin, cex.title,
        colors, lty, lwd, lpch, col, lpty, llab, cex.lab, ltxt,
        lwsymbols, lwgt, lwgts, liwgt, lsyinches, 
        ldosm, smooth, smooth.par, smooth.iter, smooth.power=1,
        ylim=if(llimy) lylim else NULL, ylimfac, ylimext,
        lrefx, lci, lrefyw, lnsims, lsimres,
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
plmatrix <- 
function(data1, data2=NULL, data=NULL, panel=l.panel, 
         nrows=0, ncols=0, save=TRUE, robrange.=FALSE, range.=NULL,
         pch=NULL, clr=1, reference=0, ltyref=3,
         log="", xaxs="r", yaxs="r", 
         vnames=NULL, main='', cex=NA, cexlab=1.3, cextext=1, cex.title=1,
         bty="o", oma=NULL, ...)
{
## Purpose:    pairs  with different plotting characters, marks and/or colors
##             showing submatrices of the full scatterplot matrix
##             possibly on several pages
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: 23 Jul 93; minor bug-fix+comments:
  ## M.Maechler
  l.panel <- function(x,y,indx,indy,pch=1,clr=1,...) {
    if (is.character(pch)) text(x,y,pch,col=clr) else
    points(x,y,pch=pch,col=clr,...)
  }
  oldpar <- par(c("mfrow","mar","cex","oma","mgp"))
  on.exit(par(oldpar))
##---------------------- preparations --------------------------
## data
  if (is.formula(data1))  data1 <- model.frame(data1,data)
  if (is.data.frame(data1)) {
    for (jj in 1:length(data1)) data1[[jj]] <- as.numeric(data1[[jj]])
    data1 <- as.matrix(data1)
  } else data1 <- cbind(data1)
#  stop("!plmatrix! first argument must either be a formula or a data.frame or matrix")
  nv1 <- dim(data1)[2]
  lv1 <- 0
  if (is.null(data2)) {
    ldata <- data1
    if (save) { nv1 <- nv1-1; lv1 <- 1 }
    nv2 <- nv1; lv2 <- 0
  } else { # cbind data2 to data for easier preparations
    save <- FALSE
    if (is.formula(data2))  data2 <- model.frame(data2,data)
    if (is.data.frame(data2)) {
      for (jj in 1:length(data2)) data2[[jj]] <- as.numeric(data2[[jj]])
      data2 <- as.matrix(data2)
    }
    ldata <- cbind(data1, as.matrix(data2))
    nv2 <- ncol(ldata)-nv1 ; lv2 <- nv1 }
  nvv <- ncol(ldata)
  tnr <- nrow(ldata)
## variable labels
  if (missing(vnames)) vnames <- dimnames(ldata)[[2]]
  if (is.null(vnames)) vnames <- paste("V",1:nvv)
## plotting characters
  if (length(pch)==0) pch <- 1
## range
  rg <- matrix(nrow=2,ncol=nvv,dimnames=list(c("min","max"),vnames))
  if(is.matrix(range.)) {
    if (is.null(colnames(range.))) {
      if (ncol(range)==ncol(rg)) rg[,] <- range.  else 
      warning('argument  range.  not suitable. ignored')
    } else {
      lj <- match(colnames(range.),vnames)
      if (any(is.na(lj))) {
        warning('variables', colnames(range.)[is.na(lj)],'not found')
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
  colnames(rg) <- vnames
## reference lines
  tjref <- (length(reference)>0)&&!(is.logical(reference)&&!reference)
  if (tjref) {
    if(length(reference)==1) lref <- rep(reference,length=nvv) else {
      lref <- rep(NA,nvv)
      lref[match(names(reference),vnames)] <- reference
    }
    names(lref) <- vnames
  }
## plot
  lpin <- par("pin")
  lnm <- if (lpin[1]>lpin[2]) c(5,6) else c(8,5)
  if (is.na(nrows)||nrows<1) nrows <- ceiling(nv1/((nv1-1)%/%lnm[1]+1))
  if (is.na(ncols)||ncols<1) ncols <- ceiling(nv2/((nv2-1)%/%lnm[2]+1))
  if (length(oma)!=4) oma <- 2 + c(0,0,!is.null(main)&&main!="",1)
  par(mfrow=c(nrows,ncols),oma=oma, mgp=c(1,0.5,0))
  ## log
  if (length(grep("y",log))>0) ldata[ldata[,1:nv1]<=0,1:nv1] <- NA
  if (length(grep("x",log))>0) ldata[ldata[,lv2+1:nv2]<=0,lv2+1:nv2] <- NA
  if (!is.na(cex)) par(cex=cex)
  cex <- par("cex")
  cexl <- cex*cexlab
  cext <- cexlab*cextext
  par(oma=oma*cexlab, mar=rep(0.2,4))
##
  npgr <- ceiling(nv1/nrows)
  npgc <- ceiling(nv2/ncols)
##----------------- plots ----------------------------
  for (ipgr in 1:npgr) {
    lr <- (ipgr-1)*nrows
  for (ipgc in 1:npgc) {
    lc <- (ipgc-1)*ncols
    if (save&&((lr+nrows)<=lc)) break
  for (jr in 1:nrows) { #-- plot row [j]
    jd1 <- lr+jr
    j1 <- lv1 + jd1
    if (jd1<=nv1)  v1 <- ldata[,j1]
    for (jc in 1:ncols) { #-- plot column  [j2-lv2] = 1:nv2
      jd2 <- lc+jc
      j2 <- lv2 + jd2
    if (jd1<=nv1 & jd2<=nv2) {
      v2 <- ldata[,j2]
      plot(v2,v1, type="n", xlab="", ylab="", axes=FALSE, 
           xlim <- rg[,j2], ylim <- rg[,j1],
           xaxs=xaxs, yaxs=yaxs, log=log)
      usr <- par("usr")
      if (jr==nrows||jd1==nv1)   mtext(vnames[j2], side=1, line=0.5, cex=cexl,
            at=mean(usr[1:2]))
      if (jc==1) mtext(vnames[j1], side=2, line=0.5, cex=cexl,
            at=mean(usr[3:4]))
      if (jr==1) axis(3,xpd=TRUE)
      if (jc==ncols||jd2==nv2) axis(4,xpd=TRUE)
      box(bty=bty)
      if (any(v1!=v2,na.rm=TRUE)) { # not diagonal
        panel(v2,v1,jd2,jd1, pch, clr, ...)
        if (tjref) abline(h=lref[j2],v=lref[j1],lty=ltyref)
      }
      else { uu <- par("usr") # diagonal: print variable name
             text(mean(uu[1:2]),mean(uu[3:4]), vnames[j1], cex=cext) }
    }
      else frame()
    }
  mtext(main,3,oma[3]-0.5,outer=TRUE,cex=cex.title)
  stamp(sure=FALSE,line=par("mgp")[1]+0.5)
  }}}
  on.exit(par(oldpar))
  "plmatrix: done"
}

## ===========================================================================
plres2x <-
  function(formula=NULL, reg=NULL, data=reg, restricted=NULL, size = 0,
  slwd = 1, scol = 2, xlab = NULL, ylab= NULL, xlim=NULL, ylim=NULL,
  main = NULL, cex.title= NULL, ...)
{
## Purpose:  plot residuals vs. two x's
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
##   size       the symbols are scaled so that 'size' is the size of
##              the largest symbol in cm.
##   main       main title, defaults to the formula
##   ...        additional arguments for the S-function 'plot'
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
#  if (length(pch)<lny) pch <- 1:lny
  lmx <- max(lpr)
  l.panel <- function(x,y,indx,indy,ly,clr, ssize) {
    lix <- indx==ly
    liy <- indy==ly
    x[!(lix|liy)] <- NA
    segments(x-ssize*lix,y-ssize*liy,x+ssize*lix,y+ssize*liy,col=clr)
    abline(1,-1,lty=3)
  }
  plmatrix(lpr, panel=l.panel, pch=ly, range=c(0,lmx), main=main, ssize=ssize)
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
dropdata <- function(data, rowid, incol="row.names")
{
  ## Purpose:   drop observations from a data frame
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 31 Mar 2004, 15:32
  if (incol=="row.names") 
    li <- match(as.character(rowid),row.names(data),nomatch=0)
  else {
    incol <- if (is.numeric(incol)) (1:ncol(data))[incol] else
    match(incol, colnames(data))
    if (is.na(incol)) stop("misspecified incolumn")
    li <- match(rowid,data[,incol],nomatch=0)
  }
  li <- li[li>0]
  ldata <- data[-li,]
  lattr <- attributes(data)
  lattr <- lattr[is.na(match(names(lattr),c("dim","dimnames","row.names")))]
  attributes(ldata) <- c(attributes(ldata),lattr)
  if (length(attr(data,"na.action")))
    li <- c(attr(data,"na.action"),
            which(naresid(data,(1:nrow(data)%in%li))))
  if (length(li)==0) warning(":dropobs: no observations dropped")
  if (length(li)==nrow(data)) warning(":dropobs: no observations left")
  class(li) <- "omit"
  attr(ldata,"na.action") <- li
  ldata
}
## ======================================================================
showd <- function(data, first=3, nrow.=4, ncol.=NULL,
                  doc=options("doc")[[1]])
{
## print some rows (and columns) of a matrix or data.frame
  lnl <- ""
  if (length(tit(data))>0) {
    cat(tit(data),":   ")
    lnl <- "\n"
  }
  if (length(dim(data))>0) {
    cat("dim: ",dim(data))
    lnl <- "\n"
  }
  cat(lnl)
  ldata <- cbind(data)
  l.nr <- nrow(ldata)
  l.nc <- ncol(ldata)
  l.ic <- if (length(ncol.)==0) 1:l.nc  else {
    if (length(ncol.)==1) {
      if (l.nc>ncol.) 
        c(seq(1,by=l.nc%/%ncol.,length=ncol.-1),l.nc) else 1:l.nc
    } else  {
      lic <- ncol.[ncol.>0&ncol<=l.nc]
      if (length(lic)>0) lic else 1:l.nc
    }
  }
  if (l.nr<=nrow.+first)  l.dc <- format(ldata[,l.ic])  else {
    l.ir <- c(1:first,round(seq(first,l.nr,length=nrow.+1))[-1])
    l.ir <- unique(c(last(l.ir,-1),l.nr))
    l.dc <- data.frame(u.merge(format(ldata[l.ir,l.ic]),"",after=first), 
                       stringsAsFactors=FALSE)
    names(l.dc) <- colnames(ldata)[l.ic]
    if (length(lrn <- row.names(ldata))>0)
      row.names(l.dc) <- c(lrn[1:first],"...", lrn[l.ir[-(1:first)]])
  }
  if (l.nc==1) {
    row.names(l.dc) <-
      format(rbind(row.names(l.dc),l.dc[,1]),justify="right")[1,]
    l.dc <- t(l.dc)
  }
  print(l.dc,quote=FALSE)
  if (length(doc)&&doc&&length(doc(data)))
    cat("\ndoc:  ",paste(doc(data),collapse="\n  "),"\n")
  invisible(l.dc)
}
## -------------------------------------------------------------------------
mframe <-
function(mfrow=NULL, mfcol=NULL, mft=NULL, row=TRUE, oma=c(0,0,2,1),
                 mar=options("mar")[[1]], mgp=options("mgp")[[1]], ...)
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
  if (length(mar)==0) mar <- c(3,3,1,1)+0.5
  if (length(mgp)==0) mgp <- c(2,0.8,0)
  invisible(if(row)
            par(mfrow=c(mfrow,mfcol), oma=oma, mar=mar, mgp=mgp, ...) else
            par(mfcol=c(mfrow,mfcol), oma=oma, mar=mar, mgp=mgp, ...) )
}
## ==========================================================================
robrange <-
function(data, trim=0.2, fac=3)
{
  ldat <- data[is.finite(data)]
  if (is.character(ldat)|length(ldat)==0) stop("!robrange! invalid data")
  trim <- c(trim, 0.2)[1]
  if (!is.finite(trim)) trim <- 0.2
  lmn <- mean(ldat,trim=trim)
  lds <- sort(abs(ldat-lmn))
  ln <- ceiling((1-trim)*length(ldat))
  if (ln<3) {
    warning(":robrange: not enough valid data. returning ordinary range")
    lsd <- Inf } else {
    lsd <- fac*sum(lds[1:ln]/(ln-1))
    if (lsd==0) {
      warning(":robrange: robust range has width 0. returning ordinary range")
      lsd <- Inf }
  }
  c(max(lmn-lsd,min(ldat)), min(lmn+lsd,max(ldat)))
}
## ==========================================================================
stamp <- function(sure=TRUE, outer.margin = NULL,
                  project=options("project")[[1]], 
                  step=options("step")[[1]], stamp=options("stamp")[[1]], ...)
{
## Purpose:   plot date and project information
## -------------------------------------------------------------------------
## Arguments:
##   sure     if F, the function only plots its thing if  options("stamp")[[1]]>0
##   outer    if T, the date is written in the outer margin
##   project  project title
##   step     title of step of data analysis
##   ...      arguments to  mtext , e.g., line=3
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: 13 Aug 96, 09:00
  if (length(stamp)==0) {
    warning(":stamp: setting options(stamp=1)")
    options(stamp=1)
    stamp <- 1
  }
  if (length(outer.margin)==0) outer.margin <- par('oma')[4]>0 
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
## ===========================================================================
last <-
function(data,n = 1)
{
  ldt <- length(data)
  data[sign(n)*((ldt-abs(n)+1):ldt)]
}
## ==============================================================
nainf.exclude <- function (object, ...)
  ## na.omit, modified to omit also Inf and NaN values
{
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
  ff <- if(inf) function(x) !is.finite(x) else is.na
  if (is.matrix(object)) apply(ff(object),2,sum)  else {
    if (is.list(object)) sapply(object,function(x) sum(ff(x)) )
    else if(is.atomic(object)) sum(ff(object))
  }
}
## ==========================================================================
logst <- function(data, mult=1)
{
  ## Purpose:   logs of data, zeros and small values treated well
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  3 Nov 2001, 08:22
  data <- cbind(data)
  lncol <- ncol(data)
  ljdt <- rep(TRUE, lncol)
  lthr <- rep(NA, lncol)
  for (lj in 1:lncol) {
    ldt <- data[,lj]
    ldp <- ldt[ldt>0&!is.na(ldt)]
    if(length(ldp)==0) ljdt[lj] <- FALSE else {
      lq <- quantile(ldp,probs=c(0.25,0.75),na.rm=TRUE)
      lthr[lj] <- lc <- lq[1]^(1+mult)/lq[2]^mult
      li <- which(ldt<lc)
      if (length(li)) 
        ldt[li] <- 10^(log10(lc) + (ldt[li]-lc)/(lc*log(10)))
      data[,lj] <- log10(ldt)
    }
  }
  if (length(names(data))) 
    lnmpd <- names(ljdt) <- names(lthr) <- names(data)  else
    lnmpd <- as.character(1:lncol)
  attr(data,"threshold") <- c(lthr)
  if (any(!ljdt)) {
    warning(':logst: no positive data for variables',lnmpd[!ljdt],
            '. These are not transformed')
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
"tit<-" <- function(x, value) ## ! argument must be 'value'. demanded by attr
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
                   doc=options("doc")[[1]])
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
getmeth <- function(fn,mt)  getS3method(as.character(substitute(fn)),
                                        as.character(substitute(mt))) 
warn <- function()   table(names(warnings()))
BR <- browser
DB <- function(on=TRUE) options(error=if(on) recover else NULL)
options(show.dummy.coef=TRUE)
## ===========================================================================
if (length(options("colors"))==0)
  options(colors) <- c("black","firebrick3","deepskyblue3","springgreen3",
                "darkgoldenrod3","olivedrab3","purple3","orange3","palegreen3")
if (length(options("colors.ra"))==0)
  options(colors.ra) <-
  c("black","gray","blue","cyan","red","magenta","darkgreen","gray30")
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
#     newform <- update.formula(object,
#                               paste(". ~ . +", paste(scope, collapse="+")))
#     data <- model.frame(update(object, newform)) # remove NAs
#     object <- update(object, data = data)
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2,
                  dimnames = list(c("<none>", scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
    n0 <- length(object$residuals)
    env <- environment(formula(object))
    for(i in seq(ns)) {
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying +", tt, "\n", sep='')
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
#    data <- model.frame(object) # remove NAs
#    object <- update(object, data = data)
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2,
                  dimnames =  list(c("<none>", scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
    n0 <- length(object$residuals)
    env <- environment(formula(object))
    for(i in seq(ns)) {
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying -", tt, "\n", sep='')
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

drop1.lm <-
  function (object, scope, scale = 0, all.cols = TRUE, test = c("none", 
    "Chisq", "F"), k = 2, ...) 
{
    ## robust regr
    if (inherits(object, 'rlm'))  {
      if (any(object$weights!=1))
        warning(':regr/drop1.lm: external weights ignored')
      object$weights <- object$w
    }
    x <- model.matrix(object)
    offset <- model.offset(model.frame(object))
    iswt <- !is.null(wt <- object$weights)
    n <- nrow(x)
    asgn <- attr(x, "assign")
    tl <- attr(object$terms, "term.labels")
    if (missing(scope)) 
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
    chisq <- deviance.lm(object)
    dfs <- numeric(ns)
    RSS <- numeric(ns)
##-     y <- object$residuals + predict(object)
    y <- object$residuals + object$fitted.values
    na.coef <- (1:length(object$coefficients))[!is.na(object$coefficients)]
    for (i in 1:ns) {
        ii <- seq_along(asgn)[asgn == ndrop[i]]
        if (all.cols) 
            jj <- setdiff(seq(ncol(x)), ii)
        else jj <- setdiff(na.coef, ii)
        z <- if (iswt) 
            lm.wfit(x[, jj, drop = FALSE], y, wt, offset = offset)
        else lm.fit(x[, jj, drop = FALSE], y, offset = offset)
        dfs[i] <- z$rank
        oldClass(z) <- "lm"
        RSS[i] <- deviance(z)
    }
    scope <- c("<none>", scope)
    dfs <- c(object$rank, dfs)
    RSS <- c(chisq, RSS)
    if (scale > 0) 
        aic <- RSS/scale - n + k * dfs
    else aic <- n * log(RSS/n) + k * dfs
    dfs <- dfs[1] - dfs
    dfs[1] <- NA
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
      } else dev <- dev/scale
      df <- aod$Df
      nas <- !is.na(df)
      dev[nas] <- pchisq(dev[nas], df[nas], lower.tail = FALSE)
      aod[, "Pr(Chi)"] <- dev
    } else  if (test == "F") {
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
    head <- c("Single term deletions", "\nModel:",
              deparse(as.vector(formula(object))), 
              if (scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
# <environment: namespace:stats>
deviance.lm <- function (object, ...) 
sum(weighted.residuals(object)^2, na.rm = TRUE)
# }## only for R version <= 2.7.1
## -----------------------------------
polrregr <- 
function (formula, data, weights, start, ..., subset, na.action, 
    contrasts = NULL, Hess = FALSE, model = TRUE, method = c("logistic", 
        "probit", "cloglog", "cauchit"))
  ## copy of polr from MASS. 1 line added by WSt to keep eta in the result
{
    logit <- function(p) log(p/(1 - p))
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
        fit <- switch(method, logistic = glm.fit(X, y1, wt, family = binomial(), 
            offset = offset), probit = glm.fit(X, y1, wt, family = binomial("probit"), 
            offset = offset), cloglog = glm.fit(X, y1, wt, family = binomial("probit"), 
            offset = offset), cauchit = glm.fit(X, y1, wt, family = binomial("cauchit"), 
            offset = offset))
        if (!fit$converged) 
            stop("attempt to find suitable starting values failed")
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
        stop("'start' is not of the correct length")
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
    if (model) 
        fit$model <- m
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    class(fit) <- "polr"
    fit
}
# <environment: namespace:MASS>
