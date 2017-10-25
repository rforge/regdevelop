lasso <- function(x,...)  UseMethod("lasso")

lasso.default <-
  function(x, y, index, subset, weights = rep(1, length(y)), model='gaussian',
           lambda=NULL, lstep = 21,
           adaptive = FALSE, cv.function = cv.lasso, cv=NULL,
           adaptcoef = NULL, adaptlambda = NULL, penscale = sqrt,
           center=NA, standardize = TRUE,
           save.x = TRUE, control = lassoControl(), ...)
{
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 29 Nov 2007, 08:22
  fit <- lassoGrpFit(x, y, index, subset=subset, weights=weights, model=model,
                     lambda=lambda, lstep=lstep,
                     penscale=penscale,
                     center=center, standardize = standardize,
                     save.x = if(adaptive) TRUE else save.x,
                     control = control, ...)
  fit$innercall <- fit$call
  fit$call <- match.call() ## Overwrite lassoGrpFit

  if (adaptive) {
    if (control@trace > 0)
      cat('\n*** calling lasso.lassogrp to adapt to results of first call\n\n')

    fit <- lasso.lassogrp(fit, lambda=lambda, lstep=lstep,
                          cv=cv, cv.function = cv.function,
                          adaptcoef = NULL, adaptlambda = NULL,
                          save.x = FALSE, control = control, ...)
    ## bookkeeping
    if (save.x) fit$x <- x
  }
  fit
}
## ==================================================================
lasso.formula <-
  function(formula, data, subset, weights, na.action, offset, nonpen = ~ 1,
           model='gaussian', lambda=NULL, lstep = 21, adaptive = FALSE,
           cv.function = cv.lasso, penscale = sqrt,
           cv=NULL, adaptcoef = NULL, adaptlambda = NULL,
           contrasts = NULL, save.x = TRUE,
           center=NA, standardize = TRUE,
           control = lassoControl(), ...)
{
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 29 Nov 2007, 08:22
  ## same as lassogrp.formula
  m <- match.call(expand.dots = FALSE)
  ## Remove not-needed stuff to create the model-frame
  m$nonpen <- m$model <- m$lambda <- m$lstep <- m$adaptive <- m$cv.function <-
    m$coef.init <- m$penscale <- m$cv <- m$adaptcoef <- m$adaptlambda <-
      m$contrasts <- m$save.x <-
        m$center <- m$standardize <- m$control <- m$... <- NULL

  ## Workaround subtle bug (which was *not* in grplasso; see
  ## "offset(" in ../tests/test_grplasso.R :
  if( ("offset" %in% all.names(formula)) &&
     !("offset" %in% all.vars (formula)))
    stop("Cannot (yet) use  offset(.) in lasso() formula;\n",
         " use 'offset' argument instead")
  
  environment(nonpen) <- environment(formula)
  
  l <- lassogrpModelmatrix(m, formula, nonpen, data, weights, subset, na.action,
                     contrasts, env = parent.frame())
##-   if(is.null(coef.init))  coef.init <- rep(0, ncol(l$x))

  ## 'subset': has been done in lassogrpModelmatrix, do not use it again below:
  fit <- if(is.null(adaptcoef))
    lasso.default(x = l$x, y = l$y, index = l$index,
                       weights = l$w,
                       model = model, lambda = lambda,
                       adaptive = adaptive,
                       cv.function = cv.function,
                       offset = l$off,
                       penscale=penscale,
                       center=center, standardize=standardize,
                       save.x = save.x, control = control, ...)
  else {
    if (is.list(adaptcoef)&&"coefficients"%in%names(adaptcoef))
      adaptcoef <- adaptcoef$coefficients
    l$model <- model
    lasso.lassogrp(l, lambda=lambda, lstep=lstep,
                   cv = cv, cv.function = cv.function,
                   adaptcoef = adaptcoef, adaptlambda = NULL,
                   penscale = penscale, ## weights = NULL,
                   center = center, standardize = standardize,
                   save.x = save.x, control = control, ...)
  }
  fit$terms <- l$Terms
  ## !!! Terms and x do not match index and coefs if adaptive==T
  fit$contrasts <- attr(l$x, "contrasts")
  fit$xlevels <- .getXlevels(l$Terms, l$mf)
  fit$na.action <- attr(l$mf, "na.action")
  fit$formula <- formula
  fit$innercall <- fit$call
  fit$call <- match.call() ## Overwrite lassoGrpFit
  structure(fit, class = "lassogrp")
}

## ==================================================================
lasso.lassogrp <- function(x, lambda=NULL, lstep=21,
                           cv = NULL, cv.function = cv.lasso,
			   adaptcoef = NULL, adaptlambda = NULL,
                           penscale = sqrt, weights = NULL,
                           center = NA, standardize = TRUE,
                           save.x = TRUE, control = lassoControl(), ...)
{ ## adaptive lasso
  if (length(x$x)==0 || length(x$y)==0)
      stop("must have an 'x' and 'y' component")
  lx <- x$x
  ly <- x$y
  lindex <- x$index
  lcall <- match.call(expand.dots=TRUE)
## coefficients to adapt to
  lcf <- adaptcoef
  if (length(lcf)==0) {
    if(length(adaptlambda) > 1)
      stop("length(adaptlambda) is larger than one")
    if (length(adaptlambda) == 0) { ## use CV
      lcv <- if (is.null(cv)) cv.function(x, se=FALSE) else cv
      adaptlambda <- x$lambda[which.min(lcv$mse)]
    }
    else if (adaptlambda < 0) {
      adaptlambda <- x$lambda[-adaptlambda]
    }
    lil <- which(abs(adaptlambda-x$lambda) <= 0.01*x$lambda)
    if (length(lil)==0)
      stop('lambda not suitable. not yet programmed') # call lassogrp
    lcf <- x$coef[,lil[1]]
  } else # length(lcf) >= 1 :
    if (length(lcf)!=ncol(lx)) stop('wrong number of coefficients in adaptcoef')
  if (is.matrix(lcf)) lcf <- structure(c(lcf), names=dimnames(lcf)[[1]])
  if (is.null(names(lcf))) {
    if (length(lcf)>1)
      stop('coefficients must have names')
    else { ## index of coef
      lcf <- x$coefficients[,lcf]
      if (is.null(lcf)) stop("'adaptcoef' not suitable")
      names(lcf) <- dimnames(x$coefficients)[[1]]
    }
  }
  ## drop variables with coef==0
  ltdrop <- aggregate(lcf==0, list(lindex), all)
  ltdr <- ltdrop[ltdrop[,2],1]
  lv <- which(!lindex%in%ltdr)
  ltrm <- x$lasso.terms
  if (length(ltdr)) {
    ltrm[ltdr+1] <- 0
    lindex <- structure(lindex[lv], term.labels=attr(lindex, 'term.labels'))
    # do not drop term.labels!
  }
  if (all(lindex<=0))
    stop('no nonzero coefficients available to adapt to')
  lxx <- lx <- lx[,lv,drop=FALSE]
  lsc <- lcf <- lcf[lv]
  ## same scale within groups
  lgrp <- table(lindex)>1
  if (any(lgrp)) {
    lgrp <- as.numeric(names(lgrp)[lgrp])
    for (li in lgrp) {
      lig <- which(lindex==li)
      lsc[lig] <- sqrt(sum(lcf[lig]^2))
    }
  }
  lsc <- lsc[lindex>0]
  lj <- match(names(lsc),dimnames(lx)[[2]])
  if (any(is.na(lj)))
    stop(paste("coefficients", names(lsc)[is.na(lj)]," not in model matrix"))
## --- change scales
  lx[,lj] <- sweep(lx[,lj,drop=FALSE],2,lsc,"*")
##-   attr(lindex, 'grpnames') <-
##-     grpnames  ## contains names even for eliminated groups
  if (is.null(weights)) weights <- x$weights
  fit <- lassoGrpFit(lx,ly, index=lindex, weights=weights,
                     model=x$model, lambda=lambda, lstep=lstep,
                     penscale=penscale, standardize=FALSE,
                     save.x = save.x, control = control, ...)
## adjust scales
  fit$coefficients[lj,] <- sweep(fit$coefficients[lj, ,drop=FALSE],
                                 1, lsc, "*")
  if (save.x) fit$x <- lxx
  fit$adaptcoef <- lcf

## terms
  lt <- fit$lasso.terms[!is.na(fit$lasso.terms)]
  ltrm[names(lt)] <- lt
  fit$lasso.terms <- ltrm
  if(!is.null(x$formula)) { ## i.e. *NOT* normally (called from lasso.default()):
    ## formula (without dropped terms)
    ltr <- setdiff(names(ltrm)[(!is.na(ltrm))&ltrm!=0],"(Intercept)")
    fit$formula <- update(x$formula, paste('~', paste(ltr, collapse="+")))
  }
  fit$call <- lcall
##  fit
  structure(fit, class = "lassogrp")
}

## =========================================================
## lassogrp <- function(x, ...)
##   UseMethod("lassogrp")

## lassogrp.formula <-
##   function(x, data, subset, weights = NULL, na.action,
##            model='gaussian', nonpen = ~ 1,
##            lambda = NULL, lfac = 2^seq(-1,-10),
## ##         coef.init=NULL,
## ##         penscale = sqrt, center=NA, standardize = TRUE,
##            contrasts = NULL,
## ##         control = lassoControl(),
##            ...)
## {
##   ## Purpose:
##   ## ----------------------------------------------------------------------
##   ## Arguments:
##   ## ----------------------------------------------------------------------
##   ## Author: Lukas Meier, Date: 27 Jun 2006, 14:52

##   m <- match.call(expand.dots = FALSE)
##   ## Remove not-needed stuff to create the model-frame
##   m$nonpen <- m$lambda <- m$coef.init <- m$penscale <- m$model <-
##     m$standardize <- m$contrasts <- m$control <- m$... <- NULL
##   l <- lassogrpModelmatrix(m, x, nonpen, data, weights, subset, na.action,
##                      contrasts, parent.frame())

##   fit <- lassogrp.default(x = l$x, y = l$y, index = l$index, weights = l$w,
##                           model = model, offset = l$off, lambda = lambda,
## ##-                           coef.init = coef.init,
## ##-                           penscale = penscale,
## ##-                           standardize = standardize, center=center,
##                           grpnames = attr(l$index,'grpnames'),
## ##-                           control = control,
##                           ...)
##   ## subsetting has been done in lassogrpModelmatrix, do not use it in call again

##   fit$terms <- l$Terms
##   fit$contrasts <- attr(l$x, "contrasts")
##   fit$xlevels <- .getXlevels(l$Terms, l$mf)
##   fit$na.action <- attr(l$mf, "na.action")
##   fit$call <- match.call() ## Overwrite lassogrp.default
##   structure(fit, class = "lassogrp")
## }
## ===================================================================
## lassogrp.default <-
lassoGrpFit <-
  function(x, y, index, subset, weights = rep(1, length(y)), model='gaussian',
           offset = rep(0, length(y)), lambda = NULL, lstep = 21,
           coef.init = rep(0, ncol(x)), penscale = sqrt,
           center = NA, standardize = TRUE,
           save.x = NULL, control = lassoControl(), ...)
{
  ## Purpose: Function to fit a solution (path) of a group lasso problem
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x: design matrix (including intercept), already rescaled and
  ##    possibly blockwise orthonormalized.
  ## y: response vector
  ## index: vector which defines the grouping of the variables. Components
  ##        sharing the same number build a group. Non-penalized
  ##        coefficients are marked with "NA".
  ## subset:  either logical or index vector that selects rows.
  ##          needed at this low level mainly for cross validation
  ## weights: vector of observation weights.
  ## offset: vector of offset values; needs to have the same length as the
  ##         response vector.
  ## lambda: vector of penalty parameters. Optimization starts with the
  ##         first component. See details below.
  ## coef.init: initial vector of parameter estimates corresponding to the
  ##            first component in the vector "lambda".
  ## penscale: rescaling function to adjust the value of the penalty
  ##           parameter to the degrees of freedom of the parameter group.
  ##           See the reference below.
  ## model: an object of class "lassoModel" implementing the negative
  ##        log-likelihood, gradient, hessian etc. See the documentation
  ##        of "lassoModel" for more details.
  ## control: options for the fitting algorithm, see "lassoControl".
  ## ...: additional arguments to be passed to the functions defined
  ##      in "model".
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 30 Aug 2005, 09:02

  ## Do some argument checking first

  if (is.character(model)) {
    lmn <- c(gaussian='LinReg', binomial='LogReg', poisson='PoissReg')[model]
    if (is.na(lmn)) {
      warning('unsuitable argument "model"')
      lmn <- 'LinReg'
    }
    model <- get(lmn)(...)
  } else if(length(list(...)))## catch typos, etc (!!)
      warning("arguments in ", deparse(list(...)), " are disregarded")

  ## Check the design matrix
  if(!is.matrix(x)) stop("x has to be a matrix")
  p <- ncol(x)
  if(any(is.na(x))) stop("Missing values in x not allowed!")
  ## Check the response
  if(!is.numeric(y)) stop("y has to be of type 'numeric'")
  nobs <- length(y)
  if(nrow(x) != nobs) stop("x and y have not correct dimensions")
  if(!model@check(y)) stop("y has wrong format")

  ## Check the other arguments
  if(length(weights)==0) weights <- rep(1,nobs) else
  if(length(weights)!=nobs)
    stop("length(weights) not equal length(y)")
  if(any(weights < 0))
    stop("negative weights not allowed")

  if(length(offset) != nobs)
    stop("length(offset) not equal length(y)")

  if(missing(coef.init) || is.null(coef.init))  coef.init <- rep(0, p)
  else if(length(coef.init) != p)
    stop("length(coef.init) not equal ncol(x)")

  if(!is.numeric(index))
    stop("argument  'index'  has to be of type 'numeric'!")
  if(length(index) != p)
    stop("length(index) not equal ncol(x)!")
  ## !!! WSt: the following statement might cause a bug
  if(any(iina <- is.na(index))) ## recode 'NA' to '0' , for back compatibility,
    ## as in grplasso,  NA's (only!) where used to indicate "non.pen."
    index[iina] <- 0L
  if(all(in0 <- index <= 0))
    stop("None of the predictors are penalized.")
  if(!is.function(penscale))
    stop("'penscale' is not a function")
  check <- validObject(control) ## will stop the program if error occurs

  ## subset
  if(!missing(subset)) {
    if(is.logical(subset)) {
      if(length(subset)!=nobs)
        stop("length of logical vector  subset  not equal length(y)")
    } else {
      if (!(all(range(c( subset,1,nobs))==c(1,nobs)) ||
            all(range(c(-subset,1,nobs))==c(1,nobs))))
        stop("argument  'subset'  not suitable")
    }
    x <- x[subset, ,drop=FALSE]
    nobs <- length(y <- y[subset])
    weights <- weights[subset]
    offset <- offset[subset]
  }

  ## lambda
  if (length(lambda)) lambda <- lambda[!is.na(lambda)]
  if (length(lambda)==0) {
    lambdamax <-
      lambdamax(x, y, index, weights=weights, offset = offset,
                coef.init = coef.init, penscale = penscale,
                model = model, standardize = standardize, ...)
    if (length(lstep)==1)
      lstep <- if(lstep<0.5)
        seq(1,0,-max(lstep,0.001)) else seq(1,0,length=max(3,lstep))^2 # !!!
    lambda <- lambdamax*lstep
  }
  lambda <- unique(lambda)
  if(any(diff(lambda)>0)) {
    warning("lambda values will be sorted in decreasing order")
    lambda <- sort(lambda, decreasing=TRUE)
  }

  ## Extract the control information
  update.hess  <- control@update.hess
  update.every <- control@update.every
  inner.loops  <- control@inner.loops
  line.search  <- control@line.search
  max.iter     <- control@max.iter
  lower        <- control@lower
  upper        <- control@upper
  if (is.null(save.x)) save.x <- control@save.x
  tol          <- control@tol
  trace        <- control@trace
  beta         <- control@beta
  sigma        <- control@sigma

  nrlambda <- length(lambda)
  if(nrlambda > 1 && update.hess == "always"){
    warning("More than one lambda value and update.hess = \"always\". You may want to use update.hess = \"lambda\"")
  }

  ## For the linear model, the Hessian is constant and has hence to be
  ## computed only *once*

  if(model@name == "Linear Regression Model"){
    if(update.hess != "lambda"){
      update.hess <- "lambda"
      if(trace >= 1)
        cat("Setting update.hess = 'lambda'\n")
    }
    if(update.every <= nrlambda){
      update.every <- nrlambda + 1
      if(trace >= 1)
        cat("Setting update.every = length(lambda) + 1\n")
    }
  }

  ## keep original x if wanted
  x.old <- if(save.x) x else NULL

  ## Which are the non-penalized parameters?
  any.notpen    <- any(in0, na.rm=TRUE)
  inotpen.which <- which(in0)
  nrnotpen      <- length(inotpen.which)
  interc.which  <-
    if(all(x[,1] == 1)) 1L ## = 99% of cases
    else which(apply(x==1, 2, all))
  has.interc <- length(interc.which) > 0
  notpenintonly <- nrnotpen==1 && has.interc
  if (is.na(center) || is.null(center)) center <- has.interc
  if(center) {
    if (!notpenintonly)
      warning('penalization not adjusted to non-penalized carriers')
    if (!has.interc)
      message("centering without intercept .. ")
  } else { ## not center
    if(notpenintonly)
      warning('Are you sure you want uncentered carriers with model with intercept?')
  }
  ## Index vector of the penalized parameter groups
  if(any.notpen) {
    ipen <- index[-inotpen.which]
    ipen.which <- split((1:p)[-inotpen.which], ipen)
  } else {
    if(has.interc)
      warning("All groups are penalized, including the intercept.")
    ipen <- index
    ipen.which <- split((1:p), ipen)
  }

  nrpen      <- length(ipen.which)
  dict.pen   <- sort(unique(ipen))

  ## Table of degrees of freedom
  ipen.tab   <- table(ipen)[as.character(dict.pen)]

  ## Center
  if (center) {
    if(has.interc) {
      ctr <- colMeans(x. <- x[,-interc.which, drop = FALSE])
      x[,-interc.which] <- sweep(x., 2, ctr)
    } else {
      x[] <- sweep(x, 2, colMeans(x))
    }
  }
  ## Standardize the design matrix -> blockwise orthonormalization
  if(standardize) {
    ##warning("...Using standardized design matrix.\n")
    stand        <- varBlockStand(x, ipen.which, inotpen.which)
    x            <- stand$x
    scale.pen    <- stand$scale.pen
    scale.notpen <- stand$scale.notpen
  }
  ## From now on x is the *normalized* design matrix!

  ## Extract the columns into lists, works faster for large matrices
  if(any.notpen){
    x.notpen <- list(); length(x.notpen) <- nrnotpen
    for(i in seq_along(inotpen.which))
      x.notpen[[i]] <- x[,inotpen.which[[i]], drop = FALSE]
  }

  x.pen <- list(); length(x.pen) <- length(nrpen)
  for(i in seq_along(ipen.which))
    x.pen[[i]] <- x[,ipen.which[[i]], drop = FALSE]

  ## Extract the needed functions
  check     <- validObject(model)
  invlink   <- model@invlink
  nloglik   <- model@nloglik
  ngradient <- model@ngradient
  nhessian  <- model@nhessian

  coef      <- coef.init
  coef.pen  <- coef.init
  if(any.notpen)
    coef.pen  <- coef[-inotpen.which]

##-   lambdanm <- format(lambda)
##-   if (any(duplicated(lambdanm))) lambdanm <- as.character(lambda)
  lambdanm <- paste('l',1:nrlambda,sep='.')
  names(lambda) <- lambdanm

  norms.pen    <- c(sqrt(rowsum(coef.pen^2, group = ipen)))

  norms.pen.m  <- matrix(0, nrow = nrpen, ncol = nrlambda,
                         dimnames = list(NULL, lambdanm))
  norms.npen.m <- matrix(0, nrow = nrnotpen, ncol = nrlambda,
                         dimnames = list(NULL, lambdanm))
  nloglik.v <- fn.val.v <- numeric(nrlambda)
  coef.m    <- grad.m   <-
               matrix(0, nrow = p, ncol = nrlambda,
                       dimnames = list(colnames(x), lambdanm))
  fitted    <- linear.predictors <-
               matrix(0, nrow = nobs, ncol = nrlambda,
                      dimnames = list(rownames(x), lambdanm))

  converged <- rep(TRUE, nrlambda)

  ## *Initial* vector of linear predictors (eta) and transformed to the
  ## scale of the response (mu)
  eta <- offset + c(x %*% coef)
  mu <- invlink(eta)

  ## Create vectors for the Hessian approximations
  if(any.notpen)
    nH.notpen <- numeric(nrnotpen)
  nH.pen <- numeric(nrpen)

  for(pos in 1:nrlambda){
    l <- lambda[pos]

    if(trace >= 2)
      cat("\nLambda:", l, "\n")

    ## Initial (or updated) Hessian Matrix of the *negative* log-likelihood
    ## function (uses parameter estimates based on the last penalty parameter
    ## value)

    if(update.hess == "lambda" && (pos %% update.every == 0 || pos == 1)) {
      ## Non-penalized groups
      if(any.notpen){
        for(j in 1:nrnotpen){ ## changed
          Xj <- x.notpen[[j]]
          nH.notpen[j] <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
        }
      }
      ## Penalized groups
      for(j in 1:nrpen){
        ind <- ipen.which[[j]]
        Xj  <- x.pen[[j]]
        diagH <- numeric(length(ind))
        for(i in seq_along(ind)){
          diagH[i] <- nhessian(Xj[, i, drop = FALSE], mu, weights, ...)
        }
        nH.pen[j] <- min(max(diagH, lower), upper)
      }
    }

    ## Start the optimization process
    fn.val <- nloglik(y, eta, weights, ...) +
      l * sum(penscale(ipen.tab) * norms.pen)

    ## These are needed to get into the while loop the first time
    do.all <- FALSE
    d.fn   <- d.par <- 1

    counter    <- 1 ## Count the sub-loops
    iter.count <- 0 ## Count the loops through *all* groups

    ## Stop the following while loop if the convergence criterion is fulfilled
    ## but only if we have gone through all the coordinates

    ##while(d.fn > tol | d.par > sqrt(tol) | !do.all){
    while(d.fn > tol || d.par > sqrt(tol) || !do.all){
      ## Escape loop if maximal iteration reached
      if(iter.count >= max.iter){
        converged[pos] <- FALSE
        warning("Maximal number of iterations reached for lambda[", pos, "]")
        break
      }

      ## Save the parameter vector and the function value of the previous step
      fn.val.old <- fn.val
      coef.old   <- coef

      ## Check whether we have some useful information from the previous step

      ## Go through all groups if counter == 0 or if we have exceeded the
      ## number of inner loops (inner.loops)
      if(counter == 0 || counter > inner.loops){
        do.all <- TRUE
        guessed.active <- 1:nrpen
        counter <- 1
        if(trace >= 2)
          cat("...Running through all groups\n")
      }else{## Go through the groups which were identified at the previous step
        guessed.active <- which(norms.pen != 0)
        if(length(guessed.active) == 0){
          guessed.active <- 1:nrpen
          do.all <- TRUE
          if(trace >= 2)
            cat("...Running through all groups\n")
        }else{
          do.all <- FALSE
          if(counter == 1 && trace >= 2)
            cat("...Starting inner loop\n")
          counter <- counter + 1
        }
      }
      if(do.all)
        iter.count <- iter.count + 1

      ## These are used for the line search, start at initial value 1
      ## They are currently here for security reasons
      start.notpen <- rep(1, nrnotpen)
      start.pen    <- rep(1, nrpen)

      if(any.notpen){
        ## Optimize the *non-penalized* parameters
        for(j in 1:nrnotpen){
          ind <- inotpen.which[j]
          Xj  <- x.notpen[[j]]

          ## Gradient of the negative log-likelihood function
          ngrad <- c(ngradient(Xj, y, mu, weights, ...))

          ## Update the Hessian if necessary
          if(update.hess == "always"){
            nH <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
          }else{
            nH <- nH.notpen[j]
          }

          ## Calculate the search direction
          d <- -(1 / nH) * ngrad
          ## Set to 0 if the value is very small compared to the current
          ## coefficient estimate

          d <- zapsmall(c(coef[ind], d))[2]

          ## If d != 0, we have to do a line search
          if(d != 0){
            scale <- min(start.notpen[j] / beta, 1) ##1
            coef.test      <- coef
            coef.test[ind] <- coef[ind] + scale * d

            Xjd <- Xj * d
            eta.test     <- eta + Xjd * scale

            if(line.search){
              qh    <- sum(ngrad * d)

              fn.val0      <- nloglik(y, eta, weights, ...)
              fn.val.test  <- nloglik(y, eta.test, weights, ...)

              qh <- zapsmall(c(qh, fn.val0))[1]

              ## Armijo line search. Stop if scale gets too small (10^-30).
              while(fn.val.test - fn.val0 > sigma * scale * qh & scale > 10^-30){
                ##cat("Doing line search (nonpen)\n")
                scale          <- scale * beta
                coef.test[ind] <- coef[ind] + scale * d
                eta.test       <- eta + Xjd * scale
                fn.val.test    <- nloglik(y, eta.test, weights, ...)
              }
            } ## end if(line.search)
            if(scale <= 10^-30){ ## Do nothing in that case
              ## coef.test <- coef
              ## eta.test  <- eta
              ## mu        <- mu
              start.notpen[j] <- 1
            }else{ ## Update the information
              coef <- coef.test
              eta  <- eta.test
              mu   <- invlink(eta)
              start.notpen[j] <- scale
            }

            ## Save the scaling factor for the next iteration (in order that
            ## we only have to do very few line searches)
            ## start.notpen[j] <- scale

            ## Update the remaining information
            ## coef <- coef.test
            ## eta  <- eta.test
            ## mu   <- invlink(eta)
          } ## end if(abs(d) > sqrt(.Machine$double.eps))
        } ## end for(j in 1:nrnotpen)
      } ## if(any.notpen)

      ## Optimize the *penalized* parameter groups
      for(j in guessed.active){
        ind  <- ipen.which[[j]]
        npar <- ipen.tab[j]

        coef.ind       <- coef[ind]
        cross.coef.ind <- crossprod(coef.ind)

        ## Design matrix of the current group
        Xj <- x.pen[[j]]

        ## Negative gradient of the current group
        ngrad <- c(ngradient(Xj, y, mu, weights, ...))

        ## Update the Hessian if necessary
        nH <-
          if(update.hess == "always") {
            dH <- numeric(length(ind))
            for(i in seq_along(ind)){ ## for loop seems to be faster than sapply
              dH[i] <- nhessian(Xj[,i,drop = FALSE], mu, weights, ...)
            }
            min(max(dH, lower), upper)
          }
          else nH.pen[j]

        cond       <- -ngrad + nH * coef.ind
        cond.norm2 <- crossprod(cond)

        ## Check the condition whether the minimum is at the non-differentiable
        ## position (-coef.ind) via the condition on the subgradient.

### __ FIXME  use 'border' instead of  'l * penscale(npar)'  below -__________

        border <- penscale(npar) * l
        if(cond.norm2 > border^2){
          d <- (1 / nH) *
            (-ngrad - l * penscale(npar) * (cond / sqrt(cond.norm2)))
          ##d <- zapsmall(c(coef.ind, d))[-(1:npar)]
        }else{
          d <- -coef.ind
        }

        ## If !all(d == 0), we have to do a line search
        if(!all(d == 0)){
          scale <- min(start.pen[j] / beta, 1)

          coef.test      <- coef
          coef.test[ind] <- coef.ind + scale * d
          Xjd            <- c(Xj %*% d)
          eta.test       <- eta + Xjd * scale

          if(line.search){
            qh <- sum(ngrad * d) +
              l * penscale(npar) * sqrt(crossprod(coef.ind + d)) -
                l * penscale(npar)* sqrt(cross.coef.ind)

            fn.val.test    <- nloglik(y, eta.test, weights, ...)
            fn.val0        <- nloglik(y, eta, weights, ...)

            left <- fn.val.test - fn.val0 +
              l  * penscale(npar) * sqrt(crossprod(coef.test[ind])) -
                l  * penscale(npar) * sqrt(cross.coef.ind)
            right <- sigma * scale * qh
            while(left > right & scale > 10^-30){
              ##cat("Doing line search (pen)\n")
              scale          <- scale * beta
              coef.test[ind] <- coef.ind + scale * d
              eta.test       <- eta + Xjd * scale
              fn.val.test    <- nloglik(y, eta.test, weights, ...)

              left <- fn.val.test - fn.val0 +
                l  * penscale(npar) * sqrt(crossprod(coef.test[ind])) -
                  l  * penscale(npar) * sqrt(cross.coef.ind)
              right <- sigma * scale * qh
            } ## end while(left > right & qh != 0)
          } ## end if(line.search)
          ## If we escaped the while loop because 'scale' is too small
          ## (= we add nothing), we just stay at the current solution to
          ## prevent tiny values
          if(scale <= 10^-30){ ## Do *nothing* in that case
            ##coef.test <- coef
            ##eta.test  <- eta
            ##mu        <- mu
            start.pen[j] <- 1
          }else{
            coef <- coef.test
            eta  <- eta.test
            mu   <- invlink(eta)
            start.pen[j] <- scale
          }
        } ## end if(!all(d == 0))
        norms.pen[j] <- sqrt(crossprod(coef[ind]))
      } ## end for(j in guessed.active)

      fn.val <- nloglik(y, eta, weights, ...) +
        l * sum(penscale(ipen.tab) * norms.pen)

      ## Relative difference with respect to parameter vector
      d.par <- sqrt(crossprod(coef - coef.old)) / (1 + sqrt(crossprod(coef)))

      ## Relative difference with respect to function value (penalized
      ## likelihood)
      d.fn <- (fn.val.old - fn.val) / (1 + abs(fn.val))

      ## Print out improvement if desired (trace >= 2)
      if(trace >= 2){
        cat("d.fn:", d.fn, " d.par:", d.par,
            " nr.var:", sum(coef != 0), "\n")
      }

      ## If we are working on a sub-set of predictors and have converged
      ## we stop the optimization and will do a loop through all
      ## predictors in the next run. Therefore we set counter = 0.

      ##if(d.fn <= tol && d.par <= sqrt(tol)){
      if(d.fn <= tol && d.par <= sqrt(tol)){
        counter <- 0 ## will force a run through all groups
        if(trace >= 2 && !do.all)
          cat("...Subproblem (active set) solved\n")
      }
    } ## end of while(d.fn > tol | d.par > sqrt(tol) | !do.all)

    if(trace == 1)
      cat("Lambda:", l, " nr.var:", sum(coef != 0), "\n")

    coef.m[,pos]            <- coef
    fn.val.v[pos]           <- fn.val
    norms.pen.m[,pos]       <- norms.pen
    nloglik.v[pos]          <- nloglik(y, eta, weights, ...)
    grad.m[,pos]            <- ngradient(x, y, mu, weights, ...)
    linear.predictors[,pos] <- eta
    fitted[,pos]            <- invlink(eta)
  } ## end for(pos in 1:nrlambda)
## ---
  nterms <- max(abs(index))
  terml <- c('(Intercept)', attr(index, 'term.labels'))
  nterms <- max(abs(index)+1,length(terml))
  if (length(terml)<nterms)
    terml <- c('(Intercept)', paste('G',1:(nterms-1),sep='.') )
  dimnames(norms.pen.m)[[1]] <- terml[dict.pen+1]
  termt <- rep(NA,nterms)
  termt[dict.pen+1] <- 1
  termt[unique(-index[inotpen.which])+1] <- -1
  names(termt) <- terml
  ## Transform the coefficients back to the original scale if the design
  ## matrix was standardized
  if(standardize) {
    if(any.notpen)
      coef.m[inotpen.which,] <- coef.m[inotpen.which, ,drop=FALSE] / scale.notpen
    ## For df > 1 we have to use a matrix inversion to go back to the
    ## original scale
    for(j in seq_along(ipen.which)){
      ind <- ipen.which[[j]]
      coef.m[ind,] <- solve(scale.pen[[j]], coef.m[ind,,drop = FALSE])
    }
  }
  ## correct intercept for centering
  if(center && has.interc)
    coef.m[interc.which,] <-
      coef.m[interc.which,] - ctr %*% coef.m[-interc.which, ,drop=FALSE]

  structure(list(coefficients = coef.m,
                 norms.pen    = norms.pen.m,
                 nloglik      = nloglik.v,
                 fn.val       = fn.val.v,
                 fitted       = fitted,
                 linear.predictors = linear.predictors,
                 call         = match.call(),
                 x            = x.old, ## use untransformed values
		 y	      = y,
                 index        = index,
                 lasso.terms  = termt,
                 weights      = weights,
                 model        = model,
                 offset       = offset,
                 lambda       = lambda,
                 penscale     = penscale,
                 center=center, standardize=standardize,
                 ngradient    = grad.m,
                 converged    = converged,
                 control      = control),
            class = "lassogrp")
} ## {lassoGrpFit}

## ==================================================================

print.lassogrp <- function(x, coefficients=TRUE, doc = options("doc")[[1]], ...)
{
  ## Purpose: Print an object of class "lassogrp"
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x:      Object of class "lassogrp"
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel
  if (length(doc)==0) doc <- 0
  if (doc >= 1) {
    if (length(tit(x))) cat("\nlasso: ",tit(x),"\n")
    if (doc >= 2 && length(descr(x))) cat(paste(descr(x),"\n"),"\n")
  }

  cat("Call:\n  ", deparse(x$call), sep = "")
  cat("\n* Number of observations:", length(x$y))
  cat("\n\n* Penalty parameter lambda:  ")
  lambda <- x$lambda
  if (length(lambda)>5)
    cat(length(lambda),' values between ',
        format(min(lambda)), ' and ',
        format(max(lambda)), '\n')  else{
          cat("\n")
          print(lambda)
        }
  if (!is.null(x$adaptcoef)) {
    cat("* Adaptation coefficients: \n")
    print(x$adaptcoef) ## cat('  ',format(x$adaptcoef), '\n')
  }
  cat("\n* Predictors: groups     :",
      length(unique(na.omit(x$index))), "\n")
  llt <- x$lasso.terms
  lltn <- names(llt)
  cat("* Penalized:\n")
  print(lltn[!is.na(llt) & llt>0], quote=FALSE)
  cat("  Not penalized:\n")
  print(lltn[!is.na(llt) & llt<0], quote=FALSE)
  if (any(ldr <- !is.na(llt) & llt==0)) {
    cat("  Not included:\n")
    print(lltn[ldr], quote=FALSE)
  }
  if (coefficients) {
    p <- ncol(lcf <- cbind(x$coefficients))
    lsel <- p > 5
    ljlam <- if (lsel) round(seq(1,p,length=5)) else 1:p
    cat("\n* Coefficients", if(lsel) "for selected lambdas",
        "(*N* = not penalized) :\n")
    lind <- x$index
    lord <- c(which(lind== 0),
              which(lind < 0)[order(lind[lind<0])],
              which(lind > 0)[order(lind[lind>0])])
    lcf <- lcf[lord, ljlam, drop=FALSE]
    lind <- lind[lord]
    dimnames(lcf)[[1]] <-
      ##- paste(c('not pen',lnm)[lind+1], dimnames(lcf)[[1]], sep=": ")
      paste(ifelse(lind <= 0, '*N* ', '    '), dimnames(lcf)[[1]], sep=" ")
    print(rbind(lambda=x$lambda[ljlam], lcf))
  }
  invisible(x)
}

## ==================================================================
predict.lassogrp <- function(object, newdata = NULL,
                             type = c("link", "response"),
                             na.action = na.pass, ...)
{
  ## Purpose: Obtains predictions from a "lassogrp" object.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## object:  a "lassogrp" object
  ## newdata: data.frame or design matrix of observations at which
  ##          predictions are to be made.
  ## type: the type of prediction. type = "link" is on the
  ##       scale of linear predictors, whereas type = "response" is on
  ##       the scale of the response variable, i.e. type = "response"
  ##       applies the inverse link function on the linear predictors.
  ## ...:  other options to be passed to the predict function.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  7 Apr 2006, 09:04

  type <- match.arg(type)

  na.act <- object$na.action

  ## If no new data is available, use the information in the fit object
  if(missing(newdata) || is.null(newdata)) {
    pred <- switch(type,
                   link = object$linear.predictors,
                   response = fitted(object))
    if(!is.null(na.act))
      pred <- napredict(na.act, pred)
  } else {
    tt <- object$terms
    if (is.null(tt)) tt <- object$formula
    if (!is.null(tt)) { ## if we have a terms object in the fit
      newdata <- as.data.frame(newdata)
      Terms <- delete.response(tt)
      m <- model.frame(Terms, newdata, na.action = na.action,
                       xlev = object$xlevels)
      offset <- attr(tt, "offset")

      if(!is.null((cl <- attr(Terms, "dataClasses"))))
        .checkMFClasses(cl, m)
      x <- model.matrix(Terms, m, contrasts = object$contrasts)
      if (ncol(x)!=nrow(coef(object))) {
        warning('ncol(x)!=nrow(coef(object)): check!')
        x <- x[,row.names(coef(object))]
      }
      pred <- x %*% coef(object)
      if(!is.null(offset)){
        offset <- eval(attr(tt, "variables")[[offset]], newdata)
        pred <- pred + offset
      }
    } else { ## if the object comes from lassoGrpFit
      ## Hmm, this does happen (FIXME?)
      x <- as.matrix(newdata)
      pred <- x %*% coef(object)
      if(any(object$offset != 0))
        warning("Possible offset not considered!")
    }

    pred <- switch(type,
                   link = pred,
                   response = object$model@invlink(pred)
                   ##apply(pred, 2, object$model@invlink))
                   )
    if(!is.null(na.action))
      pred <- napredict(na.action, pred)
  }
  if(dim(pred)[2] == 1)
    pred <- structure(c(pred), names=row.names(pred))

  attr(pred, "lambda") <- object$lambda
  pred
}

## ==================================================================
fitted.lassogrp <- function(object, ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 26 Jun 2006, 12:11

  out <- object$fitted
  attr(out, "lambda") <- object$lambda
  out
}

## ==================================================================
"[.lassogrp" <- function(x, i) {

  ## First get dimensions of the original object x

  nrlambda <- length(x$lambda)

  if(missing(i))
      i <- seq_len(nrlambda)

  ## Error checking
  ## ...

  ## Subset the object
  fit.red <- x

  fit.red$coefficients <- coef(x)[,i,drop = FALSE]
  if(length(fit.red$coefficients) == 0)
    stop("Not allowed to remove everything!")

  fit.red$call$lambda  <- fit.red$lambda <- x$lambda[i]
  fit.red$ngradient    <- x$ngradient[,i,drop = FALSE]
  fit.red$nloglik      <- x$nloglik[i]
  fit.red$fitted       <- fitted(x)[,i]
  fit.red$linear.predictors <- x$linear.predictors[,i,drop = FALSE]
  fit.red$fn.val       <- x$fn.val[i]
  fit.red$fn.val       <- x$fn.val[i]
  fit.red$norms.pen    <- x$norms.pen[,i,drop = FALSE]
  fit.red
}


## =========================================================
plot.lassogrp <-
  function(x, type = c("norms","coefficients","criteria"),
           cv = NULL, se = TRUE,
           col = NULL, lty=NULL, lwd = 1.5, mar = NULL, axes = TRUE, 
           ylim = NULL,
           legend = TRUE, main = NULL, xlab = "Lambda", ylab = NULL, ...)
{
  ## Purpose: Plots the solution path of a "lassogrp" object. The x-axis
  ##          is the penalty parameter lambda, the y-axis can be
  ##          coefficients or the l2-norm of the coefficient groups.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x:    a lassogrp object
  ## type: what should be on the y-axis? Coefficients (dummy
  ##       variables)?
  ## col:  a vector indicating the color of the different solution
  ##       paths. The length should equal the number of coefficients.
  ## ...:  other parameters to be passed to the plotting functions.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  7 Apr 2006 / Werner Stahel Aug 2009

  lambda <- x$lambda
  if(length(lambda) == 1)
    stop("Plot function not available for a single lambda")

  type <- match.arg(type)

  xlim <- rev(range(sqrt(lambda)))
  ind <- unique(x$index)

  nr.npen <- sum(x$index<=0)
  dict.pen <- na.omit(ind)
  dict.pen.ord <- nr.npen + 1 : length(dict.pen)
## ----------------------------
  if (type == "norms" || type == "coefficients") {
    if (type == "norms") {
      yy <- x$norms.pen
    ## colors and ltypes
      col <- if(is.null(col))  1:nrow(yy)  else rep(col,length=nrow(yy))
      lty <- if(is.null(lty))  1:nrow(yy)  else rep(lty,length=nrow(yy))
      if (is.null(main)) main <- "Paths of norms of penalized coefficients"
    } else { # if(type == "coefficients")
      if (is.null(main)) main <- "Coefficient paths"
      col <- if(is.null(col))  1:length(ind)  else rep(col,length=length(ind))
      lty <- if(is.null(lty))  1:length(ind)  else rep(lty,length=length(ind))
      ## possibly too complicated
      index.ord <- numeric()
      index.ord[x$index<=0] <- 1:nr.npen

      dict.pen.sort <- sort(dict.pen)
      for(j in seq_along(dict.pen))
        index.ord[x$index == dict.pen.sort[j]] <- dict.pen.ord[j]

      col <- col[index.ord]
      lty <- lty[index.ord]
      yy <- coef(x)
    }
    if (is.null(ylab)) ylab <- type
    if (is.null(ylim)) ylim <- range(yy[,lambda>0])
    matplot(sqrt(lambda), t(yy), type = "l", axes = FALSE,
            xlab = xlab, ylab = ylab, col = col, lty=lty,
            main = main, xlim = xlim, ylim = ylim, lwd=lwd, ...)
    if (axes[length(axes)]) axis(4)
#    axis(4, mgp=c(3,2,0), at = yy[, ncol(yy)], labels = rownames(yy))
    if (legend) legend(x="topleft", legend=rownames(yy), col=col, lty=lty)
## ----------------------------
  } else if(type=='criteria') {

    lmain <- if (is.null(main)) "log-likelihood and penalty" else main
    col <- if(is.null(col))  c(1,2,2,3)  else rep(col,length=4)
    lty <- if(is.null(lty))  c(1,5,6,2)  else rep(lty,length=4)
    mar <- if(is.null(mar))  c(4,4,4,4)  else rep(mar,length=4)
    l1 <- colSums(x$norms.pen)
    yy <- cbind(x$nloglik/length(x$y))
    lylab <- if(is.null(ylab)) "-loglik/n" else ylab
    if (is.null(cv)) cv <- x$cv
    if (!is.null(cv)) {
      if (is.logical(cv)&&cv) cv <- cv.lasso(x)
      if (!is.list(cv)) {
        warning('argument  "cv"  not suitable')
      } else {
        yy <- cbind(yy, cv$mse)
        if (se) yy <- cbind(yy, cv$mse+outer(cv$mse.se,c(-1,1)))
        if (is.null(ylab)) lylab <- paste(lylab, " || MSE (cv)")
        if (is.null(main)) lmain <- "log-likelihood, MSE (cv), and penalty"
      }
    }
    matplot(sqrt(lambda), yy, type = "l", axes=FALSE,
            xlab = xlab, ylab = lylab,
            col = col[c(1,2,3,3)], lty=lty[c(1,2,3,3)],
            main = lmain, xlim = xlim, mar=mar, lwd=lwd, ...)
    par(usr=c(par('usr')[1:2], 0,1.05*max(l1)))
    lines(sqrt(lambda), l1, col=col[4], lty=lty[4], lwd=lwd)
    if (axes[length(axes)]) {
      axis(4, col=col[4])
      mtext("l1 norm",4,1.8)
    }
  } else warning(':plot.lassogrp: invalid argument "type". No plot')
  if (axes) {
    box()
    axis(2)
    xlb <- pretty(lambda)
    axis(1, at=sqrt(xlb), labels=format(xlb))
    axis(3, at=sqrt(lambda), labels=rep('',length(lambda)), tcl=0.5,
         mgp=c(1,0.5,0), xpd=TRUE)
    ilam <- c(1,length(lambda))
    if(ilam[2]>8)
      ilam <- c(seq(1,length(lambda)-3,by=5),length(lambda))
    la <- lambda[ilam]
    axis(3, at=sqrt(la), labels=rep("",length(la)), tcl=-0.3, xpd=TRUE)
    mtext(names(la),3,0.3,at=sqrt(la))
  }
  stamp(sure=FALSE)
  invisible(yy)
}

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
  lform <- formula(object)
  lcall <- object$call
  if (is.null(lform)) lform <- lcall$x
  if (is.null(lform))
    stop('extract.lassogrp needs an object with formula')
  if (is.null(i)) {
    if (is.null(lambda))
      stop('!extract.lassogrp! Either arg.  i  or  lambda  must be given')
    if (length(lambda)>1 || lambda<0 || lambda>max(object$lambda))
      stop('!extract.lassogrp! argument "lambda" not suitable')
#    lamlam <- outer(lambda,object$lambda,'/')
    i <- which.min(abs(sqrt(lambda)-sqrt(object$lambda)))
##-       unique(apply( lamlam>0.99&lamlam<1.01, 1, which))
##-     if (is.null(i))
##-       stop('!extract.lassogrp!  lambda  not equal to any lambdas in lassogrp object')
  }
  ldata <- data
  if (is.null(ldata)) ldata <- object$data
  if (is.null(ldata)) ldata <- eval(object$call$data)
  if (is.null(ldata)) ldata <- eval(eval(object$call$x)$call$data)
  if (is.null(ldata)) stop("no 'data' found")
##  environment(lform) <- ldata
  ni <- length(i)
  if (ni!=1) {
    warning('extract.lassogrp() not programmed for extracting more than 1 model')
    i <- i[1]
  }
  lpen <- object$norms.pen[,i]
  ltrms <- names(lpen[lpen>0])
  lform <- update(eval(lform),
                  if (is.null(lcall$nonpen)) ~1 else lcall$nonpen )
  lform <- update(lform, paste('~.+',paste(ltrms, collapse=' + ')) )
  result <- NULL
  lmod <- object$model
  if (!is.null(lmod)) {
    if (is.character(lmod))
      lmod <- pmatch(lmod,c('gaussian','binomial','poisson'))
    else
      lmod <- pmatch(substring(lmod@name,1,3),c('Linear','Logistic','Poisson'))
  }
  lmeth <- c("lm","binomial","glm")[lmod]
  if (!is.null(data)) lcall$data <- substitute(data)
  if (is.null(lcall$data))
    stop ("!extract.lassogrp! I do not find the data. Specify the 'data' argument")
  if (is.null(lcall$na.action)) lcall$na.action <- as.name('nainf.exclude')
  fcall <- if (fitfun=='regr')
    call(fitfun, formula=lform, data=lcall$data, method = lmeth,
               family=c('gaussian','binomial','poisson')[lmod],
               subset=lcall$subset, weights=lcall$weights,
               na.action=lcall$na.action) # , '...'=...
  else {
    if (fitfun=='lm')
    call(fitfun, formula=lform, data=lcall$data,
               subset=lcall$subset, weights=lcall$weights)
             #  ,na.action=lcall$na.action) # , '...'=...
    else
    call(fitfun, formula=lform, data=lcall$data,
               family=c('gaussian','binomial','poisson')[lmod],
               subset=lcall$subset, weights=lcall$weights,
               na.action=lcall$na.action) # , '...'=...
  }
  lreg <- eval(fcall, parent.frame())
  for (li in seq_along(i)) { ## for loop not in use!
    lr <- lreg
    lk <- i[li]
    lcoef <- lr$coefficients <- object$coefficients[,lk]
    lr$fitted.values <- object$fitted[,lk]
    lrsd <- lr$residuals <- object$y - lr$fitted.values
    lr$df.residual <- lreg$df.residual - sum(lcoef==0)
    ## only good for raw residuals
    lcl <- lreg$call
    lcl[1] <- paste('lasso',fitfun,sep='.')
    attr(lcl,'comment') <- 'call not R usable'
    lr$call <- lcl
    lrss <- sum(lrsd^2, na.rm=TRUE)
    lr$sigma <- sqrt(lrss/lr$df.residual)
    lr$r.squared <- 1-lrss/sum(object$y^2,na.rm=TRUE)
    lr$stres <- lr$testcoef <- lr$adj.r.squared <- lr$fstatistic <-
      lr$covariance <- lr$correlation <- NULL
    lr$fitfun <- fitfun
    lr$faceff <- ## lfaceff <- if (exists("factoreffects"))
      factoreffects(lr) ## else    try(dummy.coef(lr), silent=TRUE)
    ## if (is.list(lfaceff)) lr$faceff <- lfaceff
    lr$fit.unpen <- lreg

    result <- c(result,lr)
  }
  result <- if (length(result)==1) result[[1]] else result
  class(result) <- c('lassofit', class(lreg))
  result
}
## ===================================================================
print.lassofit <-
  function(x, # dummycoef = NULL, #digits = max(3, options("digits")[[1]] - 3),
    residuals=FALSE, doc = options("doc")[[1]], ...)
{
  ## Author: Werner Stahel, Date: 24 Aug 2009, 12:27
  if(length(doc)==0) doc <- 0
  if(doc >= 1) {
    if (length(tit(x))) cat("\n",tit(x),"\n")
    if (doc >= 2 && length(descr(x))) cat(paste(descr(x),"\n"),"\n")
  }
##-   if (length(dummycoef)==0)
##-     dummycoef <- c(options("show.dummy.coef")[[1]],TRUE)[1]
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),"\n", sep = "")
  cat("Fitting function ",x$fitfun,"\n\n")
  cf <- coef(x)
  print(list('nonzero coefficients'=cf[cf!=0],
             'zero coefficents'=names(cf[cf==0])))
  cat("Root Mean Square of Error: ", format(x$sigma),'\n')
  if (residuals) {
    cat("\nResiduals:\n")
    print(summary(x$residuals))
  }
  invisible(x)
}

## Needed, e.g., for  update(., )  to work:
formula.lassogrp <- ## identical to stats:::formula.glm
function (x, ...) 
{
    form <- x$formula
    if (!is.null(form)) {
        form <- formula(x$terms)
        environment(form) <- environment(x$formula)
        form
    }
    else formula(x$terms)
}
## ===================================================================
cv.lasso <-
  function (object, blocks = 10,
            trace=FALSE, control=lassoControl(trace=0), env=globalenv(),
            plot.it = NULL, se = TRUE, ...)
  ## modified cv.lars, WSt
  ## allows for blocks to be specified rather than random
{
  x <- object$x
  y <- object$y
  n <- length(y)
  blocksinmodel <- FALSE
  ## argument blocks -> generate blocks
  if (is.character(blocks)) {
    bl <- get.blocks (object, blocks, env=env)
    blocklist <- bl$blocklist
    mf <- bl$model.frame
    blocksinmodel <- blocks %in% all.vars(as.formula(object$call$formula))
  } else if (length(blocks)==1) {
    stopifnot(blocks >= 1)
    blocklist <- split(sample(1:n), rep(1:blocks, length = n))
  } else if (length(blocks)!=n) {
    stop('!cv.lasso! argument "blocks" not suitable')
  } else {
    blocklist <- split(sample(1:n), blocks)
  }
  ## fetch x
  lj <- match(names(object$index), colnames(x))
  if (any(is.na(lj)))
    stop (if(length(x)==0) "!cv.lasso! no x component" else
          "!cv.lasso! length(index) not equal ncol(x).")
  lx <- x[,lj]
  ## crossvalidate
  K <- length(blocklist)
  blockmse <- matrix(0, K, length(object$lambda))
  fitted <- matrix(0, n, length(object$lambda))
  for (i in seq(K)) {
    omit <- blocklist[[i]]
    fit <-
      lassoGrpFit(lx, y, object$index, subset=-omit, weights=object$weights,
                  model=object$model, offset=object$offset, lambda=object$lambda,
                  penscale=object$penscale, center=object$center,
                  standardize=object$standardize,
                  save.x=FALSE, control=control, ...)
##-     fit <- update(object, subset=-omit, model=FALSE,
##-                   lambda=object$lambda, adaptive=FALSE, cv.function=NULL,
##-                   save.x=FALSE, weights=object$eights, offset=object$offset,
##-                   control=control)
    if (blocksinmodel)  {
      xnew <- mf[omit, , drop = FALSE]
      xnew[,blocks] <- factor(mf[-omit,blocks][1])
    } else {
      xnew <- lx[omit, , drop = FALSE]
      fit$terms <- NULL
    }
    pred <- rbind(predict(fit, xnew))
       ## out of sample predictions for all lambda values
    blockmse[i,] <-
      if (blocksinmodel) {
        pred <- sweep(pred, 2, colMeans(pred) - mean(y[omit]))
        apply(y[omit] - pred, 2, var)
      }
      else
        colMeans((y[omit] - pred)^2)
    fitted[omit,] <- pred
    if (trace)  cat("\n CV Fold", i, "\n\n")
  }
  ## summarize
  mse <- colMeans(blockmse)
  mse.se <- sqrt(apply(blockmse, 2, var)/K)
  result <- list(mse = mse, mse.se = mse.se, mse.blocks = blockmse,
                 fitted = fitted, lambda = object$lambda,
                 blocksinmodel = blocksinmodel)
  ## plot
  if (is.null(plot.it)) plot.it <- identical(parent.frame(), globalenv())
  if (plot.it)   plot.lassogrp(object, type='crit', cv = result, se = se)
  ##
  invisible(structure(result, class="lassoCV"))
}
## ===================================================================
lassogrpModelmatrix <-
  function(m, formula, nonpen = ~1, data, weights, subset, na.action,
           contrasts, env)
{
  ## Purpose:  Generates the model matrix and adds information about
  ## variables to be penalized or non-penalized by the L1 term in the
  ## lasso fitting.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 30 Jun 2006: create.design() in helpers.R

  if(!inherits(formula, "formula") || length(formula) != 3)
    stop("Argument 'formula' is of wrong type or length")
  any.nonpen <- !is.null(nonpen)
  ## Case where some variables won't be penalized -> merge formulas,
  ## also check the environments (is this the correct way ???)

  m$formula <- ## <-> lasso.formula() has first arg 'formula' !
    if(any.nonpen) {
      if(!inherits(nonpen, "formula"))
        stop("Argument 'nonpen' must be a formula")

      is.gEnv <- function(e) identical(e, .GlobalEnv)

      ## Paste the two formulas together ## !!! changed
      nonp <- ~ . + dummy
      nonp[[2]][[3]] <- nonpen[[length(nonpen)]]
      fterms <- terms(formula, data=eval(m$data))
      f <- update(fterms, nonp)
      ## Get the environment of the formulas
      env.formula <- environment(formula)
      env.nonpen <- environment(nonpen)

      ## If env. of 'formula' is not global, check if env. of 'nonpen' differs.
      ## If yes give warning and use env. of 'formula'
      if(!is.gEnv(env.formula)){
        ## environment(f) <- env.formula
        if(!is.gEnv(env.nonpen) && !identical(env.formula, env.nonpen))
          warning("'formula' and 'nonpen' have different environments. I use environment of 'formula'")

        environment(f) <- environment(formula)
      } else if (!is.gEnv(env.nonpen)){ ## if env. of 'nonpen' is not global
        ## (but env. of 'formula' is), use env. of 'nonpen'
        environment(f) <- env.nonpen
      }
      f
    }
    else formula

  m$drop.unused.levels <- TRUE

  ## Create model-frame
  m[[1]] <- as.name("model.frame")

  mf <- eval(m, env)

  ## Create design matrix, na.action handles the missing values, therefore
  ## weights and offset may be of shorter length (mf does this for us)
  Terms <- attr(mf, "terms") ##terms(m$formula, data = data)
  x <- model.matrix(Terms, data = mf, contrasts = contrasts)
  y <- model.response(mf)
  w <- model.weights(mf)
  off <- model.offset(mf)

  if (!is.null(w) && any(w < 0))
    stop("Negative weights not allowed")

  if(!is.numeric(off))
    off <- rep(0, length(y))
  if(!length(w))
    w <- rep(1, length(y))

  ## Handle the non-penalized coefficients
  if(any.nonpen){
    tmp <- terms(nonpen, data = data)
    co <- contrasts[attr(tmp, "term.labels")]
    used.co <- if(length(co)) co[!unlist(lapply(co, is.null))] # else NULL
    ## also uses the response...to be changed
    x.nonpen <- model.matrix(tmp, data = mf, contrasts = used.co)
    matches <- match(colnames(x.nonpen), colnames(x))
  }
  index <- attr(x, "assign")
  names(index) <- dimnames(x)[[2]]
  if(any.nonpen)
    index[matches] <- - index[matches] # was 'NA'
  attr(index,'term.labels') <- attr(Terms,'term.labels')
  list(x = x,
       y = y,
       w = w,
       off = off,
       mf = mf,
       index = index,
       Terms = Terms)
}

## ===================================================================
get.blocks <- function(object, blocks, env=globalenv())
{
  ## Author: Werner Stahel
  m <- object$call
  if (is.null(m$formula))
    stop('!get.blocks! call of "object" does not contain a formula')
  m$formula <- update(as.formula(m$formula), paste('~.+',blocks))
  m$nonpen <- m$lambda <- m$model <- m$adaptive <- m$cv.function <-
    m$coef.init <- m$penscale <- m$standardize <- m$contrasts <- m$control <-
    m$trace <- m$... <- NULL
  m[[1]] <- as.name("model.frame")
  mf <- eval(m, env)
  if (nrow(mf)!=nrow(object$x))
    stop('!get.blocks! blocks and data incompatible')
##-   Terms <- attr(mf, "terms") ##terms(m$formula, data = data)
##-   mm <- model.matrix(Terms, data = mf, contrasts = contrasts)
##-   mmb <- model.matrix(as.formula(paste('~',blocks)),
##-                       data = mf, contrasts = contrasts)
  list(blocklist = split(1:nrow(mf), factor(mf[,blocks])), model.frame = mf)
}

## ==================================================================
varBlockStand <- function(x, ipen.which, inotpen.which)
{
  ## Purpose:   standardize matrix, respecting blocks of variables
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  4 Aug 2006, 08:50

  n <- nrow(x)
  x.ort <- x
  scale.pen <- vector(mode='list', length = length(ipen.which))
  if(length(inotpen.which) > 0){
    x. <- x[,inotpen.which, drop=FALSE]
    scale.notpen <- sqrt(colMeans(x.^2))
    x.ort[,inotpen.which] <- scale(x., FALSE, scale.notpen)
  } else scale.notpen <- NULL
  cnms <- colnames(x)
  if(is.null(cnms)) cnms <- as.character(seq_len(ncol(x)))
  rt.n <- sqrt(n)
  for(j in seq_along(ipen.which)){
    p. <- length(ind <- ipen.which[[j]])
    decomp <- qr(x[,ind])
    if(decomp$rank < min(n,p.)) ## block does not have full rank
        stop("Block belonging to columns ",
             paste(cnms[ind], collapse = ", "),
             sprintf(" has qr()-rank = %d  <  min(n=%d, p'=%d)",
                     decomp$rank, n, p.))
    scale.pen[[j]] <- qr.R(decomp) * 1 / rt.n
    x.ort[,ind] <- qr.Q(decomp) * rt.n
  }
  list(x = x.ort, scale.pen = scale.pen, scale.notpen = scale.notpen)
}

## ==================================================================
lambdamax <- function(x, ...)
  UseMethod("lambdamax")

lambdamax.formula <-
  function(formula, nonpen  = ~ 1, data, weights, subset, na.action, offset,
           coef.init, penscale = sqrt, model = LogReg(),
           center = NA, standardize = TRUE, contrasts = NULL, nlminb.opt = list(), ...)
{
  ## Purpose: Function to find the maximal value of the penalty parameter
  ##          lambda
  ## ----------------------------------------------------------------------
  ##
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 20 Apr 2006, 11:24

  call <- match.call()
  m <- match.call(expand.dots = FALSE)

  ## Remove not-needed stuff to create the model-frame
  m$nonpen <- m$coef.init <- m$penscale <- m$model <-
    m$center <- m$standardize <- m$contrasts <- m$nlminb.opt <- m$... <- NULL

  l <- lassogrpModelmatrix(m, formula, nonpen, data, weights, subset, na.action,
                     contrasts, parent.frame())

  if(missing(coef.init))
    coef.init <- rep(0, ncol(l$x))

  lambdamax.default(l$x, y = l$y, index = l$index, weights = l$w,
                    offset = l$off, coef.init = coef.init,
                    penscale = penscale, model = model,
                    center = center, standardize = standardize,
                    nlminb.opt = nlminb.opt, ...)
}

## ==================================================================
lambdamax.default <-
  function(x, y, index, weights = NULL, offset = rep(0, length(y)),
           coef.init = rep(0, ncol(x)), penscale = sqrt, model = LogReg(),
           center = NA, standardize = TRUE, nlminb.opt = list(), ...)
{
  ## Purpose: Function to find the maximal value of the penalty parameter
  ##          lambda
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## X: design matrix (including intercept), already rescaled and
  ##    possibly blockwise orthonormalized.
  ## y: response vector
  ## index: vector which defines the grouping of the variables. Components
  ##        sharing the same number build a group. Non-penalized
  ##        coefficients are marked with "NA".
  ## weights: vector of observation weights.
  ## offset: vector of offset values.
  ## coef.init: initial parameter vector. Penalized groups are discarded.
  ## penscale: rescaling function to adjust the value of the penalty
  ##           parameter to the degrees of freedom of the parameter group.
  ## model: an object of class "lassoModel" implementing
  ##        the negative log-likelihood, gradient, hessian etc. See
  ##        "lassoModel" for more details.
  ## nlminb.opt: arguments to be supplied to "nlminb".
  ## ... : additional arguments to be passed to the functions defined in
  ##       model.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 20 Apr 2006, 11:24

  if(any(iina <- is.na(index))) ## recode 'NA' to '0' , for back compatibility,
    ## as in grplasso,  NA's (only!) where used to indicate "non.pen."
    index[iina] <- 0L
  any.notpen <- any(inp <- index <= 0)
  coef.npen <- coef.init[inp] ## unpenalized parameters
  inotpen.which <- which(inp)

  ## Index vector of the penalized parameter groups
  ipen <- index[!inp]

  ## Table of degrees of freedom
  dict.pen <- sort(unique(ipen))
  ipen.tab <- table(ipen)[as.character(dict.pen)]

  ## Indices of parameter groups
  ipen.which <- split((1:ncol(x))[index>0], ipen)

  ## The 'center' part was missing in earlier versions of pkg 'lassogrp'
  ## This is "cut & paste" from lassoGrpFit() above
  interc.which  <-
    if(all(x[,1] == 1)) 1L ## = 99% of cases
    else which(apply(x==1, 2, all))
  n.int <- length(interc.which)
  if(n.int > 1)
    stop("Multiple intercepts (columns of 1) in 'x'")
  has.interc <- n.int == 1
  if (is.na(center) || is.null(center)) center <- has.interc
  else if(center && !has.interc) {
    message("Couldn't find intercept column. Setting center = FALSE.")
    center <- FALSE
  }
  if(center) { ## also  has.interc (!)
    ctr <- colMeans(x. <- x[,-interc.which, drop = FALSE])
    x[,-interc.which] <- sweep(x., 2, ctr)
  }

  if(standardize) {
    stand <- varBlockStand(x, ipen.which, inotpen.which)
    x     <- stand$x
  }

  x.npen <- x[,inotpen.which, drop = FALSE]
  if (length(weights)==0) weights <- rep(1, length(y))

  nlogLikFUN <- function(par)
    model@nloglik(y, offset + x.npen %*% par, weights, ...)

  mu0 <-
    if(any.notpen){
      par0 <- do.call(nlminb, c(list(start = coef.npen,
                                     objective = nlogLikFUN), nlminb.opt))$par
      model@invlink(offset + x.npen %*% par0)
    }
    else model@invlink(offset)

  ngrad0 <- model@ngradient(x, y, mu0, weights, ...)[index>0]

  ##gradnorms <- numeric(length(dict.pen))

  gradnorms <- c(sqrt(rowsum(ngrad0^2, group = ipen))) / penscale(ipen.tab)

  ##for(j in seq_along(dict.pen)){
  ##  gradnorms[j] <- sqrt(crossprod(ngrad0[which(index == dict.pen[j])])) /
  ##    penscale(sum(index == dict.pen[j], na.rm = TRUE))
  ##}
  max(gradnorms)
}
## ==================================================================
## ==================================================================
## control
setClass("lassoControl",
         representation = representation(
           save.x       = "logical",
           update.hess  = "character",
           update.every = "numeric",
           inner.loops  = "numeric",
           line.search  = "logical",
           max.iter     = "numeric",
           tol          = "numeric",
           lower        = "numeric",
           upper        = "numeric",
           beta         = "numeric",
           sigma        = "numeric",
           trace        = "numeric"),

         prototype = list(
           save.x       = FALSE,
           update.hess  = "lambda",
           update.every = 3,
           inner.loops  = 10,
           line.search  = TRUE,
           max.iter     = 500,
           tol          = 5 * 10^-8,
           lower        = 10^-2,
           upper        = 10^9,
           beta         = 0.5,
           sigma        = 0.1,
           trace        = 0),

         validity = function(object){
           if(ceiling(object@update.every) != floor(object@update.every) ||
              object@update.every <= 0)
             return("update.every has to be a natural number")

           if(ceiling(object@inner.loops) != floor(object@inner.loops) ||
              object@inner.loops < 0)
             return("inner.loops has to be a natural number or 0")

           if(ceiling(object@max.iter) != floor(object@max.iter) ||
              object@max.iter <= 0)
             return("inner.loops has to be a natural number or greater than 0")

           if(object@beta <= 0 || object@beta >= 1)
             return("beta has to be in (0, 1)")

           if(object@sigma <= 0 || object@sigma >= 1)
             return("sigma has to be in (0, 1)")

           if(object@tol <= 0)
             return("tol has to be positive")

           if(object@lower > object@upper)
             return("lower <= upper has to hold")

           return(TRUE)
         }
)

## ==================================================================
lassoControl <- function(save.x = FALSE,
                         update.hess = c("lambda", "always"),
                         update.every = 3, inner.loops = 10,
                         line.search = TRUE, max.iter = 500,
                         tol = 5 * 10^-8, lower = 10^-2, upper = Inf,
                         beta = 0.5, sigma = 0.1, trace = 0){

  ## Purpose: Options for the Group Lasso Algorithm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## save.x: a logical indicating whether the design matrix should be saved.
  ## update.hess: should the hessian be updated in each
  ##              iteration ("always")? update.hess = "lambda" will update
  ##              the Hessian once for each component of the penalty
  ##              parameter "lambda" based on the parameter estimates
  ##              corresponding to the previous value of the penalty
  ##              parameter.
  ## inner.loops: how many loops should be done (at maximum) when solving
  ##              only the active set (without considering the remaining
  ##              predictors)
  ## tol: convergence tolerance; the smaller the more precise, see
  ##      details below.

  ## lower: lower bound for the diagonal approximation of the
  ##        corresponding block submatrix of the Hessian of the negative
  ##        log-likelihood function.
  ## upper: upper bound for the diagonal approximation of the
  ##        corresponding block submatrix of the Hessian of the negative
  ##        log-likelihood function.
  ## beta: scaling factor beta < 1 of the Armijo line search.
  ## sigma: 0 < \sigma < 1 used in the Armijo line search.
  ## trace: integer. "0" omits any output,
  ##        "1" prints the current lambda value,
  ##        "2" prints the improvement in the objective function after each
  ##        sweep through all the parameter groups and additional
  ##        information.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  1 Jun 2006, 10:02


  update.hess <- match.arg(update.hess)

  RET <- new("lassoControl",
             save.x       = save.x,
             update.hess  = update.hess,
             update.every = update.every,
             inner.loops  = inner.loops,
             line.search  = line.search,
             max.iter     = max.iter,
             tol          = tol,
             lower        = lower,
             upper        = upper,
             beta         = beta,
             sigma        = sigma,
             trace        = trace)
  RET
}
## ===========================================================================
setClass("lassoModel", representation = representation(
    invlink          = "function",
    link             = "function",
    nloglik          = "function",
    ngradient        = "function",
    nhessian         = "function",
    check            = "function",
    name             = "character",
    comment          = "character"
))

setMethod("show", "lassoModel", function(object) {
    cat("Model:", object@name, "\n")
    cat("Comment:", object@comment, "\n\n")
})

lassoModel <- function(invlink, link, nloglik, ngradient, nhessian,
                       check, name = "user-specified",
                       comment = "user-specified"){
  ## Purpose: Generates models to be used for the Group Lasso algorithm.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## invlink: a function with arguments "eta" implementing the inverse link
  ##          function.
  ## link:    a function with arguments "mu" implementing the link
  ##          function.
  ## nloglik: a function with arguments "y", "mu" and "weights"
  ##          implementing the negative log-likelihood function.
  ## ngradient: a function with arguments "x", "y", "mu" and "weights"
  ##            implementing the negative gradient of the log-likelihood
  ##            function.
  ## nhessian: a function with arguments "x", "mu" and "weights"
  ##           implementing the negative} hessian of the log-likelihood
  ##           function.
  ## name: a character name
  ## comment: a character comment
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  1 Jun 2006, 10:12

    RET <- new("lassoModel",
               invlink   = invlink,
               link      = link,
               nloglik   = nloglik,
               ngradient = ngradient,
               nhessian  = nhessian,
               check     = check,
               name      = name,
               comment   = comment)
    RET
}
## ===========================================================================
## family.R
## ===========================================================================
## Logistic Regression
LogReg <- function(){
  lassoModel(invlink   = function(eta) 1 / (1 + exp(-eta)),
             link      = function(mu) log(mu / (1 - mu)),
             nloglik   = function(y, eta, weights, ...)
               -sum(weights * (y * eta - log(1 + exp(eta)))),
             ngradient = function(x, y, mu, weights, ...)
               -crossprod(x, weights * (y - mu)),
             nhessian  = function(x, mu, weights, ...)
               crossprod(x, weights * mu * (1 - mu) * x),
             check     = function(y) all(y %in% c(0, 1)),
             name      = "Logistic Regression Model",
             comment   = "Binary response y has to be encoded as 0 and 1")
}
## ===========================================================================
## Linear Regression
LinReg <- function(){
  lassoModel(invlink  = function(eta) eta,
             link  = function(mu) mu,
             nloglik   = function(y, eta, weights, ...)
               sum(weights * (y - eta)^2),
             ngradient = function(x, y, mu, weights, ...)
               -2 * crossprod(x, weights * (y - mu)),
             nhessian  = function(x, mu, weights, ...)
               2 * crossprod(x, weights * x),
             check     = function(y) TRUE,
             name      = "Linear Regression Model",
             comment   = "Use update.hess=\"lambda\" in lassoControl because the Hessian is constant")
}
## ===========================================================================
## Poisson Regression
PoissReg <- function(){
  lassoModel(invlink    = function(eta) exp(eta),
             link       = function(mu) log(mu),
             nloglik    = function(y, eta, weights, ...)
               sum(weights * (exp(eta) - y * eta)),
             ngradient  = function(x, y, mu, weights, ...)
               -crossprod(x, weights * (y - mu)),
             nhessian   = function(x, mu, weights, ...)
               crossprod(x, weights * mu * x),
             check      = function(y) all(y >= 0) & all(y == ceiling(y)),
             name       = "Poisson Regression Model",
             comment    = "")
}
