## =========================================================
extract.lassogrp <-
  function(object, i=NULL, lambda=NULL, data=NULL, fitfun='lm', ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## does not make sense yet for more than one i
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 21 Aug 2009, 08:55
  lform <- object$formula
  if (is.null(lform))
    stop("extract.lassogrp needs an object with $formula")
  if (is.null(i)) {
    if (is.null(lambda))
      stop('!extract.lassogrp! Either arg.  i  or  lambda  must be given')
    if (length(lambda)>1 || lambda<0 | lambda>max(object$lambda))
      stop('!extract.lassogrp! argument "lambda" not suitable')
#    lamlam <- outer(lambda,object$lambda,'/')
    i <- which.min(abs(sqrt(lambda)-sqrt(object$lambda)))
##-       unique(apply( lamlam>0.99&lamlam<1.01, 1, which))
##-     if (is.null(i))
##-       stop('!extract.lassogrp!  lambda  not equal to any lambdas in lassogrp object')
  }
  ni <- length(i)
  if (ni!=1) {
    warning(':extract.lassogrp: not programmed for extracting more than 1 model')
    i <- i[1]
  }
  lpen <- object$norms.pen[,i]
  ltrms <- c(object$lasso.terms[object$lasso.terms<0],names(lpen[lpen>0]))
  lform <- update(eval(lform), paste('~',paste(ltrms, collapse=' + ')))
  result <- NULL
  lmod <- object$model
  if (!is.null(lmod)) {
    if (is.character(lmod))
      lmod <- pmatch(lmod,c('gaussian','binomial','poisson'))
    else 
      lmod <- pmatch(substring(lmod@name,1,3),c('Linear','Logistic','Poisson'))
  }
  lmeth <- c("lm","binomial","glm")[lmod]
  lcall <- object$call
  ldt <- if (is.null(data)) eval(lcall$data) else data
  if (is.null(ldt))
    stop ("!extract.lassogrp! I do not find the data. Specify the 'data' argument")
  BR()
  call <- if (fitfun=='regr')
    call(fitfun, formula=lform, data=as.name("ldt"), method = lmeth,
               family=c('gaussian','binomial','poisson')[lmod],
               subset=lcall$subset, weights=lcall$weights,
               na.action=lcall$na.action) # , '...'=...
  else {
    if (fitfun=='lm')
    call(fitfun, formula=lform, data=ldt, 
               subset=lcall$subset, weights=lcall$weights,
               na.action=lcall$na.action) # , '...'=...
    else 
    call(fitfun, formula=lform, data=as.name("ldt"), 
               family=c('gaussian','binomial','poisson')[lmod],
               subset=lcall$subset, weights=lcall$weights,
               na.action=lcall$na.action) # , '...'=...
  }
  lreg <- eval(call, parent.frame())
  lreg <- eval(call)
  for (li in seq_along(i)) { ## for loop not in use!
    lr <- lreg
    lk <- i[li]
    lcoef <- lr$coefficients <- object$coefficients[,lk]
    lr$fitted.values <- object$fitted[,lk]
    lrsd <- lr$residuals <- object$y - lr$fitted.values
    lr$df.residual <- lreg$df.residual - sum(lcoef==0)
    ## only good for raw residuals
    lcl <- lreg$call
    lcl[1] <- paste(fitfun,'lassogrp',sep='.')
    attr(lcl,'comment') <- 'call not R usable'
    lr$call <- lcl
    lrss <- sum(lrsd^2, na.rm=TRUE)
    lr$sigma <- sqrt(lrss/lr$df.residual)
    lr$r.squared <- 1-lrss/sum(object$y^2,na.rm=TRUE)
    lr$stres <- lr$testcoef <- lr$adj.r.squared <- lr$fstatistic <-
      lr$covariance <- lr$correlation <- NULL
    lr$fitfun <- 'lassogrp'
    lallcoef <- try(dummy.coef(lr), silent=TRUE)
    if (is.list(lallcoef)) lr$allcoef <- lallcoef
    
    result <- c(result,lr)
  }
  result <- if (length(result)==1) result[[1]] else result
  class(result) <- c('lassofit', class(lreg))
  result
}
## =======================================================================
cv.block <-
  function(model, blocks, comp = 'nloglik', data=NULL,
           seed=0, robust=TRUE, trace=FALSE)
{
  ## Purpose:   block cross validation and bootstrap of a model
  ##            with predictive error variability
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   model      a model object that has methods for  update  and  predict
  ##   blocks     name of the variable in  data  that defines the blocks
  ##   robust     if TRUE, means will be 20% trimmed, standard deviations
  ##              will be replaced by the robust Qn estimator
  ## Value:
  ##   ...        cv means and standard errors of model component comp
  ##   psig       predictive error standard deviation, mean of ...
  ##   psigblocks predictive error standard deviation for each block
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 16 Feb 2007, 09:30
#  on.exit(browser())
  if (is.null(data)) data <- eval(model$call$data)
  if (is.null(data)) stop("!cv.block! no data")
  if (!is.character(blocks)||!blocks%in%names(data))
    stop('cv.block! argument "blocks" unsuitable')
  lblocks <- data[,blocks] 
  if (length(lblocks)!=nrow(data))
    stop("!cv.block! argument 'blocks' unsuitable")
  lyy <- model$y
  if (is.null(lyy)) lyy <- model$residuals+model$fitted.values
  if (is.null(lyy))
    stop("!cv.block! response not found")
  if (length(lyy)!=nrow(data))
    stop("!cv.block! response unsuitable")
  if (!comp%in%names(model))
    stop("!cv.block! argument 'comp' unsuitable")
  ## ---
  lbli <- split(1:nrow(data),lblocks)
  lnbl <- length(lbli)
  ## ---
  f.scl <- if (robust) function(x) Qn(x[!is.na(x)]) else
                       function(x) sqrt(var(x,na.rm=TRUE))
  f.loc <- if (robust) function(x) mean(x,trim=0.2,na.rm=TRUE) else
                       function(x) mean(x,na.rm=TRUE)
  rcomp <- matrix(NA,lnbl,length(model[[comp]]))
  rpress <- NULL
  rfit <- NULL
  ## -----------------------
  for (lrep in 1:lnbl) {
    if (trace) cat('  run',lrep,'\n')
    li <- lbli[[lrep]]
    if (length(li)==0) { ## no data to be dropped
      rpress <- rbind(rpress,rep(NA,3))
      next
    }
    ldt <- data[-li,]
    lr <- update(model, data=ldt, weights=NULL)
    rcomp[lrep,] <- lr[[comp]]
    ## --- out of sample prediction
    ldto <- data[li,]
    lblo <- factor(ldto[,blocks])
    ldto[,blocks] <- ldt[1,blocks]
    ## this leads to "wrong" constant, which cancels for scale estimation
    lfit <- predict(lr, newdata=ldto)
    lpres <- lyy[li] - lfit
    rfit <- rbind(rfit, sweep(lfit, 2, apply(lfit, 2, mean) - mean(lyy[li])) )
    lpress <- apply(lpres,2,f.scl)
    rpress <- rbind(rpress,lpress)
  }
  rcompmn <- apply(rcomp,2,f.loc)
  dim(rcompmn) <- ldim <- dim(model[[comp]])
  rcompse <- sqrt(lnbl-1)*apply(rcomp,2,f.scl)
  dim(rcompse) <- ldim 
  lsigma <- exp(apply(log(rpress),2,f.loc))
  lsigsd <- lsigma*apply(log(rpress),2,f.scl)
  rr <- list(lsigma, lsigsd, exp(rpress), rcompmn, rcompse, rcomp, fit=rfit)
  names(rr)[1:6] <- c(t(outer(c('sigma',comp),c('cvmean','cvse','cvblock'),
                              paste,sep='.')))
  attr(rr,"robust") <- robust
  rr
}


lassogrp.formula <-
  function(formula, data, subset, weights = rep(1, length(y)), na.action,
           model='gaussian', nonpen = ~ 1, 
           lambda = NULL, lfac = 2^seq(-1,-10),
##         coef.init=NULL, 
##         penscale = sqrt, center=NA, standardize = TRUE,
           contrasts = NULL, 
##         control = lassoControl(),
           ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 27 Jun 2006, 14:52

#  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  ## Remove not-needed stuff to create the model-frame
  m$nonpen <- m$lambda <- m$coef.init <- m$penscale <- m$model <-
    m$standardize <- m$contrasts <- m$control <- m$... <- NULL
  l <- create.design(m, formula, nonpen, data, weights, subset, na.action,
                     contrasts, parent.frame())
  
  fit <- lassogrp.default(x = l$x, y = l$y, index = l$index, weights = l$w,
                          model = model, offset = l$off, lambda = lambda,
##-                           coef.init = coef.init,
##-                           penscale = penscale, 
##-                           standardize = standardize, center=center,
                          grpnames = attr(l$index,'grpnames'), 
##-                           control = control,
                          ...)
  ## subsetting has been done in create.design, do not use it in call again
  
  fit$terms <- l$Terms
  fit$contrasts <- attr(l$x, "contrasts")
  fit$xlevels <- .getXlevels(l$Terms, l$mf)
  fit$na.action <- attr(l$mf, "na.action")
  fit$call <- match.call() ## Overwrite lassogrp.default 
  structure(fit, class = "lassogrp")
}

lm.fit
function (x, y, offset = NULL, method = "qr", tol = 1e-07, singular.ok = TRUE, 
    ...) 
{
    if (is.null(n <- nrow(x))) 
        stop("'x' must be a matrix")
    if (n == 0L) 
        stop("0 (non-NA) cases")
    p <- ncol(x)
    if (p == 0L) {
        return(list(coefficients = numeric(0L), residuals = y, 
            fitted.values = 0 * y, rank = 0, df.residual = length(y)))
    }
    ny <- NCOL(y)
    if (is.matrix(y) && ny == 1) 
        y <- drop(y)
    if (!is.null(offset)) 
        y <- y - offset
    if (NROW(y) != n) 
        stop("incompatible dimensions")
    if (method != "qr") 
        warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            method), domain = NA)
    if (length(list(...))) 
        warning("extra arguments ", paste(names(list(...)), sep = ", "), 
            " are just disregarded.")
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    z <- .Fortran("dqrls", qr = x, n = n, p = p, y = y, ny = ny, 
        tol = as.double(tol), coefficients = mat.or.vec(p, ny), 
        residuals = y, effects = y, rank = integer(1L), pivot = 1L:p, 
        qraux = double(p), work = double(2 * p), PACKAGE = "base")
    if (!singular.ok && z$rank < p) 
        stop("singular fit encountered")
    coef <- z$coefficients
    pivot <- z$pivot
    r1 <- seq_len(z$rank)
    dn <- colnames(x)
    if (is.null(dn)) 
        dn <- paste("x", 1L:p, sep = "")
    nmeffects <- c(dn[pivot[r1]], rep.int("", n - z$rank))
    r2 <- if (z$rank < p) 
        (z$rank + 1L):p
    else integer(0L)
    if (is.matrix(y)) {
        coef[r2, ] <- NA
        coef[pivot, ] <- coef
        dimnames(coef) <- list(dn, colnames(y))
        dimnames(z$effects) <- list(nmeffects, colnames(y))
    }
    else {
        coef[r2] <- NA
        coef[pivot] <- coef
        names(coef) <- dn
        names(z$effects) <- nmeffects
    }
    z$coefficients <- coef
    r1 <- y - z$residuals
    if (!is.null(offset)) 
        r1 <- r1 + offset
    qr <- z[c("qr", "qraux", "pivot", "tol", "rank")]
    colnames(qr$qr) <- colnames(x)[qr$pivot]
    c(z[c("coefficients", "residuals", "effects", "rank")], list(fitted.values = r1, 
        assign = attr(x, "assign"), qr = structure(qr, class = "qr"), 
        df.residual = n - z$rank))
}
