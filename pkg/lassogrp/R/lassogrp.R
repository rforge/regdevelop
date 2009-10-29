lasso <- function(x,...)  UseMethod("lasso")
  
lasso.default <-
  function(x, y, index, subset, model='gaussian', lambda=NULL, lstep = 21,
           adaptive = FALSE, cv.function = cv.lasso,
           save.x = TRUE, ...)
{
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 29 Nov 2007, 08:22
  fit <- lassogrp.default(x, y, index, subset=subset, model=model, 
                          lambda=lambda, save.x = save.x, ...)
  fit$innercall <- fit$call
  fit$call <- match.call() ## Overwrite lassogrp.default
  fit
}
## ==================================================================
lasso.formula <-
  function(x, data, subset, weights, na.action, model='gaussian',
           offset, nonpen = ~ 1, lambda=NULL, lstep = 21, adaptive = FALSE, 
           cv.function = cv.lasso, contrasts = NULL, save.x = TRUE,
           control = lassoControl(), ...)
{
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 29 Nov 2007, 08:22
  ## same as lassogrp.formula
  m <- match.call(expand.dots = FALSE)
  ## Remove not-needed stuff to create the model-frame
  m$nonpen <- m$lambda <- m$model <- m$adaptive <- m$cv.function <- 
    m$coef.init <- m$penscale <- m$standardize <- m$contrasts <- m$control <-
    m$... <- NULL
  l <- lassogrpModelmatrix(m, x, nonpen, data, weights, subset, na.action,
                     contrasts, env = parent.frame())
##-   if(is.null(coef.init))  coef.init <- rep(0, ncol(l$x))

  fit <- lasso.default(x = l$x, y = l$y, index = l$index, weights = l$w,
                       model = model, lambda = lambda, adaptive = adaptive,
                       cv.function = cv.function, offset = l$off,
##-                           penscale = penscale, 
##-                           standardize = standardize, center=center,
                       save.x = save.x, control = control, ...)
  ## subsetting has been done in lassogrpModelmatrix, do not use it in call again
  if (adaptive) {
    if (control@trace>0)
      cat('\n*** calling lasso.lassogrp to adapt to results of first call\n\n')
    fit <- lasso.lassogrp(fit, cv.function = cv.function, save.x = FALSE,
                          control = control)
##- bookkeeping
    if (save.x) fit$x <- l$x
  }
  fit$terms <- l$Terms
  ## !!! Terms and x do not match index and coefs if adaptive==T
  fit$contrasts <- attr(l$x, "contrasts")
  fit$xlevels <- .getXlevels(l$Terms, l$mf)
  fit$na.action <- attr(l$mf, "na.action")
  fit$innercall <- fit$call
  fit$call <- match.call() ## Overwrite lassogrp.default
  structure(fit, class = "lassogrp")
}

## ==================================================================
lasso.lassogrp <- function(x, lambda=NULL, lstep=21, cv=NULL, 
                  adaptcoef = NULL, adaptlambda = NULL, ...)
{ # adaptive lasso
  if (length(x$x)==0|length(x$y)==0) {
    m <- match.call(expand.dots = FALSE)
    m$lambda <- m$coef.init <- m$penscale <- m$model <-
      m$standardize <- m$contrasts <- m$control <- m$... <- NULL
    lds <- lassogrpModelmatrix(m, formula=x$formula, env=parent.frame())
    lx <- cbind(lds$x)
    ly <- lds$y
    lindex <- lds$index
  } else {
    lx <- x$x
    ly <- x$y
    lindex <- x$index
  }
## coefficients to adapt to
  lcf <- adaptcoef
  if (length(lcf)==0) {
    if (length(adaptlambda)>1)
      warning('length(adaptlambda >1. I only use first element')
    if (length(adaptlambda)==0) {
      lcv <- if (is.null(cv)) cv.lasso(x, se=FALSE) else cv
      adaptlambda <- x$lambda[which.min(lcv$rmse)]
    }
    else if (adaptlambda<0) {
      adaptlambda <- x$lambda[-adaptlambda]
    }
    lil <- which(abs(adaptlambda-x$lambda)<=0.01*x$lambda)
    if (length(lil)==0) 
      stop('lambda not suitable. not yet programmed') # call lassogrp
    lcf <- x$coef[,lil]
  } else # if (length(lcf)==1)
    if (length(lcf)!=ncol(lx)) stop('wrong number of coefficients in  adaptcoef')
  if (length(names(lcf))==0) 
    if (length(lcf)>1)
      stop('coefficients must have names') else {
        lcf <- x$coefficients[,lcf]
        names(lcf) <- dimnames(x$coefficients)[[1]]
      }
  ## drop variables with coef==0
  ltdrop <- aggregate(lcf==0, list(lindex), all)
  lt <- ltdrop[!ltdrop[,2],1]
  lv <- which(lindex%in%lt)
  ltrm <- x$lasso.terms
  ltrm[ltdrop[ltdrop[,2],1]+1] <- 0
  lindex <- lindex[lv]

  if (all(lindex<=0))
    stop('no nonzero coefficients available to adapt to')
  lx <- lx[,lv,drop=FALSE]
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
  rslt <- lassogrp.default(lx,ly, index=lindex, model=x$model,
                           standardize=FALSE, ...)
## adjust scales
  rslt$coefficients[lj,] <- sweep(rslt$coefficients[lj,],1,lsc,"/")
  rslt$adaptcoef <- lcf
  lt <- rslt$lasso.terms[!is.na(rslt$lasso.terms)]
  ltrm[names(lt)] <- lt
  rslt$lasso.terms <- ltrm
  rslt
}

## =========================================================
lassogrp <- function(x, ...)
  UseMethod("lassogrp")

lassogrp.formula <-
  function(x, data, subset, weights = NULL, na.action,
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

  m <- match.call(expand.dots = FALSE)
  ## Remove not-needed stuff to create the model-frame
  m$nonpen <- m$lambda <- m$coef.init <- m$penscale <- m$model <-
    m$standardize <- m$contrasts <- m$control <- m$... <- NULL
  l <- lassogrpModelmatrix(m, x, nonpen, data, weights, subset, na.action,
                     contrasts, parent.frame())
  
  fit <- lassogrp.default(x = l$x, y = l$y, index = l$index, weights = l$w,
                          model = model, offset = l$off, lambda = lambda,
##-                           coef.init = coef.init,
##-                           penscale = penscale, 
##-                           standardize = standardize, center=center,
                          grpnames = attr(l$index,'grpnames'), 
##-                           control = control,
                          ...)
  ## subsetting has been done in lassogrpModelmatrix, do not use it in call again
  
  fit$terms <- l$Terms
  fit$contrasts <- attr(l$x, "contrasts")
  fit$xlevels <- .getXlevels(l$Terms, l$mf)
  fit$na.action <- attr(l$mf, "na.action")
  fit$call <- match.call() ## Overwrite lassogrp.default 
  structure(fit, class = "lassogrp")
}
## ===================================================================
lassogrp.default <-
  function(x, y, index, subset, weights = rep(1, length(y)), model='gaussian',
           offset = rep(0, length(y)), lambda = NULL, lstep = 21,
           coef.init = rep(0, ncol(x)), penscale = sqrt,
           center=NULL, standardize = TRUE, 
           save.x = NULL, save.y = NULL, control = lassoControl(), ...)
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
    model <- get(lmn)()
  }

  if(is.null(standardize)) standardize <- TRUE # the default
  
  ## Check the design matrix
  if(!is.matrix(x))
    stop("x has to be a matrix")

  if(any(is.na(x)))
    stop("Missing values in x not allowed!")

  ## Check the response
  if(!is.numeric(y))
    stop("y has to be of type 'numeric'")

  nobs <- length(y)
  
  if(NROW(x) != nobs)
    stop("x and y have not correct dimensions")
  
  if(!model@check(y))
    stop("y has wrong format")

  ## Check the other arguments
  if(length(weights) != nobs)
    stop("length(weights) not equal length(y)")

  if(any(weights < 0))
    stop("Negative weights not allowed")
  
  if(length(offset) != nobs)
    stop("length(offset) not equal length(y)")

  if(is.null(coef.init))  coef.init <- rep(0, NCOL(x))
  if(length(coef.init) != NCOL(x))
    stop("length(coef.init) not equal ncol(x)")

  if(!is.numeric(index))
    stop("argument  'index'  has to be of type 'numeric'!")
  if(length(index)!=NCOL(x))
    stop("length(index) not equal ncol(x)!")

  if(all(index<=0))
    stop("None of the predictors are penalized.")
  
  check <- validObject(control) ## will stop the program if error occurs

  ## subset
  if(!missing(subset)) {
    if(is.logical(subset)) {
      if(length(subset)!=nobs)
        stop("length of logical vector  subset  not equal length(y)")
    } else {
      if (!(all(range(c(subset,1,nobs))==c(1,nobs))|
            all(range(c(-subset,1,nobs))==c(1,nobs))))
        stop("argument  'subset'  not suitable")
    }
    x <- cbind(x)[subset,]
    y <- y[subset]
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
  if (is.null(save.x)) save.x       <- control@save.x 
  if (is.null(save.y)) save.y       <- control@save.y
  tol          <- control@tol
  trace        <- control@trace
  beta         <- control@beta
  sigma        <- control@sigma
  
  nrlambda <- length(lambda)
  ncolx    <- ncol(x)
  nrowx    <- nrow(x)

  if(nrlambda > 1 & update.hess == "always"){
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
    if(update.every <= length(lambda)){
      update.every <- length(lambda) + 1
      if(trace >= 1)
        cat("Setting update.every = length(lambda) + 1\n")
    }
  }
  
  ## keep original x if wanted
  x.old <- if(save.x) x else NULL

  ## Which are the non-penalized parameters?
  any.notpen    <- any(index<=0)
  inotpen.which <- which(index<=0)
  nrnotpen      <- length(inotpen.which)
  interc.which  <- which(apply(x==1,2,all))
  notpenintonly <- nrnotpen==1 && length(interc.which)>0
  if (is.null(center)) center <- length(interc.which)>0
  if (notpenintonly&!center)
    warning('Are you sure you want uncentered carriers with model with intercept?')
  if (center&!notpenintonly)
    warning('penalization not adjusted to non-penalized carriers')
  
  ## Index vector of the penalized parameter groups
  if(any.notpen) {
    ipen <- index[-inotpen.which]
    ipen.which <- split((1:ncolx)[-inotpen.which], ipen)
  } else {
    if (length(interc.which)) 
    warning("All groups are penalized, including the intercept.")
    ipen <- index
    ipen.which <- split((1:ncolx), ipen)
  }

  nrpen      <- length(ipen.which)
  dict.pen   <- sort(unique(ipen))
  
  ## Table of degrees of freedom
  ipen.tab   <- table(ipen)[as.character(dict.pen)]
  
  ## Center
  if (center) {
    ctr <- apply(x[,-interc.which],2,mean)
    x[,-interc.which] <- sweep(x[,-interc.which],2,ctr)
  }
  ## Standardize the design matrix -> blockwise orthonormalization
  if(standardize){
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
    for(i in 1:length(inotpen.which))
      x.notpen[[i]] <- x[,inotpen.which[[i]], drop = FALSE]
  }
  
  x.pen <- list(); length(x.pen) <- length(nrpen)
  for(i in 1:length(ipen.which))
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
  lambdanm <- paste('l',1:length(lambda),sep='.')
  names(lambda) <- lambdanm
  
  norms.pen    <- c(sqrt(rowsum(coef.pen^2, group = ipen)))

  norms.pen.m  <- matrix(0, nrow = nrpen, ncol = nrlambda,
                         dimnames = list(NULL, lambdanm))
  norms.npen.m <- matrix(0, nrow = nrnotpen, ncol = nrlambda,
                         dimnames = list(NULL, lambdanm))
  nloglik.v <- fn.val.v <- numeric(nrlambda)
  coef.m    <- grad.m   <-
               matrix(0, nrow = ncolx, ncol = nrlambda,
                       dimnames = list(colnames(x), lambdanm))
  fitted    <- linear.predictors <-
               matrix(0, nrow = nrowx, ncol = nrlambda,
                      dimnames = list(rownames(x), lambdanm))

  converged <- rep(TRUE, nrlambda)
  
  ## *Initial* vector of linear predictors (eta) and transformed to the
  ## scale of the response (mu)
  eta <- offset + c(x %*% coef)
  mu <- invlink(eta)

  ## Create vectors for the Hessian approximations
  if(any.notpen){
    nH.notpen <- numeric(nrnotpen)
  }
  nH.pen <- numeric(nrpen)

  for(pos in 1:nrlambda){
    l <- lambda[pos]

    if(trace >= 2)
      cat("\nLambda:", l, "\n")

    ## Initial (or updated) Hessian Matrix of the *negative* log-likelihood
    ## function (uses parameter estimates based on the last penalty parameter
    ## value)

    if(update.hess == "lambda" & pos %% update.every == 0 | pos == 1){
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
        for(i in 1:length(ind)){
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
    while(d.fn > tol | d.par > sqrt(tol) | !do.all){
      ## Escape loop if maximal iteration reached
      if(iter.count >= max.iter){
        converged[pos] <- FALSE
        warning(paste("Maximal number of iterations reached for lambda[", pos,
                      "]", sep = ""))
        break
      }
      
      ## Save the parameter vector and the function value of the previous step
      fn.val.old <- fn.val
      coef.old   <- coef

      ## Check whether we have some useful information from the previous step

      ## Go through all groups if counter == 0 or if we have exceeded the
      ## number of inner loops (inner.loops)
      if(counter == 0 | counter > inner.loops){
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
          if(counter == 1 & trace >= 2)
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
        if(update.hess == "always"){
          diagH <- numeric(length(ind))
          for(i in 1:length(ind)){ ## for loop seems to be faster than sapply
            diagH[i] <- nhessian(Xj[,i,drop = FALSE], mu, weights, ...)
          }
          nH <- min(max(diagH, lower), upper)
        }else{
          nH <- nH.pen[j]
        }
        
        cond       <- -ngrad + nH * coef.ind
        cond.norm2 <- crossprod(cond)
        
        ## Check the condition whether the minimum is at the non-differentiable
        ## position (-coef.ind) via the condition on the subgradient.
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
      
      ##if(d.fn <= tol & d.par <= sqrt(tol)){
      if(d.fn <= tol & d.par <= sqrt(tol)){
        counter <- 0 ## will force a run through all groups
        if(trace >= 2 & !do.all)
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
  if (length(terml)<nterms+1)
    terml <- c('(Intercept)', paste('G',1:nterms,sep='.') )
  dimnames(norms.pen.m)[[1]] <- terml[dict.pen+1]
  termt <- rep(NA,nterms+1)
  termt[dict.pen+1] <- 1
  termt[unique(-index[inotpen.which])+1] <- -1  
  names(termt) <- terml
  ## Transform the coefficients back to the original scale if the design
  ## matrix was standardized
  if(standardize) {
    if(any.notpen)
      coef.m[inotpen.which,] <- coef.m[inotpen.which,] / scale.notpen
    ## For df > 1 we have to use a matrix inversion to go back to the
    ## original scale
    for(j in 1:length(ipen.which)){
      ind <- ipen.which[[j]]
      coef.m[ind,] <- solve(scale.pen[[j]], coef.m[ind,,drop = FALSE])
    }
  }
  ## correct intercept for centering
  if(center) coef.m[interc.which,] <-
    coef.m[interc.which,]- ctr%*%coef.m[-interc.which,]

  if(!save.y)
    y <- NULL

  out <- list(coefficients = coef.m,
              norms.pen    = norms.pen.m,
              nloglik      = nloglik.v,
              fn.val       = fn.val.v,
              fitted       = fitted,
              linear.predictors = linear.predictors,
              call         = match.call(),
              x = x.old, ## use untransformed values
              y = y, 
              index        = index,
              lasso.terms  = termt, 
              weights      = weights,
              model        = model,
              offset       = offset,
              lambda       = lambda,
              penscale     = penscale,
              ngradient    = grad.m,
              converged    = converged,
              control      = control)
  structure(out, class = "lassogrp")
}
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
  if (doc>=1) if (length(tit(x)))
    cat("\nlasso: ",tit(x),"\n")
  if (doc>=2) if (length(descr(x)))
    cat(paste(descr(x),"\n"),"\n")

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
#    cat('  ',format(x$adaptcoef), '\n')
    print(x$adaptcoef)
  }
  cat("\n* Predictor groups      :",
      length(unique(na.omit(x$index))), "\n")
  llt <- x$lasso.terms
  lltn <- names(llt)
  cat("* Penalized predictors  :\n")
    print(lltn[(!is.na(llt))&llt>0], quote=FALSE)
  cat("  not penalized:\n")
    print(lltn[(!is.na(llt))&llt<0], quote=FALSE)
  ldr <- (!is.na(llt))&llt==0
  if (any(ldr)) {
    cat("  not included:\n")
    print(lltn[ldr], quote=FALSE)
  }
  if (coefficients) {
    lcf <- cbind(x$coefficients)
    lsel <- ncol(lcf)>5
    ljlam <- if (lsel) round(seq(1,ncol(lcf),length=5)) else 1:ncol(lcf)
    cat("\n* Coefficients", if(lsel) "for selected lambdas",
        "(*: p = penalized) :\n")
    lind <- x$index
    lord <- order(lind)
    lcf <- lcf[lord, ljlam, drop=FALSE]
    lind <- lind[lord]
    dimnames(lcf)[[1]] <-
##-       paste(c('not pen',lnm)[lind+1], dimnames(lcf)[[1]], sep=": ")
      paste(ifelse(lind<=0,'   ','  p'), dimnames(lcf)[[1]], sep=" ")
    print(rbind('  * lambda'=x$lambda[ljlam], lcf))
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
  }  else {
  if (!is.null(tt <- object$terms)) { ## if we have a terms object in the fit
    newdata <- as.data.frame(newdata)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action,
                     xlev = object$xlevels)
    offset <- attr(tt, "offset")

    if(!is.null((cl <- attr(Terms, "dataClasses"))))
      .checkMFClasses(cl, m)
    x <- model.matrix(Terms, m, contrasts = object$contrasts)
    pred <- x %*% coef(object)
    if(!is.null(offset)){
      offset <- eval(attr(tt, "variables")[[offset]], newdata)
      pred <- pred + offset
    }
  } else { ## if the object comes from lassogrp.default
    x <- as.matrix(newdata)
    pred <- x %*% coef(object)
    if(any(object$offset != 0))
      warning("Possible offset not considered!")
  }
    
  pred <- switch(type,
                 link = pred,
                 response = object$model@invlink(pred))
                 ##apply(pred, 2, object$model@invlink))

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
"[.lassogrp" <- function(x, i){
  
  ## First get dimensions of the original object x

  nrlambda <- length(x$lambda)

  if(missing(i))
    i <- 1:nrlambda

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
           col = NULL, lty=NULL, mar=NULL, main=NULL, cv=NULL, se=TRUE,
           ylim = NULL, legend = TRUE, ...)
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
  xlb <- pretty(lambda)

  ind <- unique(x$index)

  nr.npen <- sum(x$index<=0)
  dict.pen <- na.omit(ind)
  dict.pen.ord <- nr.npen + 1 : length(dict.pen)
## ----------------------------    
  if (type == "norms" | type == "coefficients") {
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
      for(j in 1:length(dict.pen))
        index.ord[x$index == dict.pen.sort[j]] <- dict.pen.ord[j]
        
      col <- col[index.ord]
      lty <- lty[index.ord]
      yy <- coef(x)
    }
    if (is.null(ylim)) ylim <- range(yy[,lambda>0])
    matplot(sqrt(lambda), t(yy), type = "l", axes = FALSE,
            xlab = "Lambda", ylab = type, col = col, lty=lty, 
            main = main, xlim = xlim, ylim = ylim, ...)
    box()
    axis(1, at=sqrt(xlb), labels=as.character(xlb))
    axis(3, at=sqrt(lambda), labels=rep('',length(lambda)), tcl=0.5, xpd=TRUE)
    la <- lambda[c(1,length(lambda))]
    mtext(names(la),3,0.5,at=sqrt(la))
    axis(2)
    axis(4)
#    axis(4, mgp=c(3,2,0), at = yy[, ncol(yy)], labels = rownames(yy))
    if (legend) legend(x="topleft", legend=rownames(yy), col=col, lty=lty)
## ----------------------------    
  } else if(type=='criteria') {
    
    if (is.null(main)) main <- "log-likelihood and penalty"
    col <- if(is.null(col))  c(1,2,2,3)  else rep(col,length=4)
    lty <- if(is.null(lty))  c(1,5,6,2)  else rep(lty,length=4)
    mar <- if(is.null(mar))  c(4,4,4,4)  else rep(mar,length=4)
    l1 <- apply(x$norms.pen,2,sum)
    yy <- cbind(x$nloglik/length(x$y))
    if (!is.null(cv)) {
      if (is.logical(cv)&&cv) cv <- cv.lasso(x)
      if (!is.list(cv)) {
        warning('argument  "cv"  not suitable')
      } else {
        yy <- cbind(yy, cv$rmse)
        if (se) yy <- cbind(yy, cv$rmse+outer(cv$rmse.se,c(-1,1)))
      }
    }
    matplot(sqrt(lambda), yy, type = "l", axes=FALSE,
            xlab = "Lambda", ylab = "-loglik",
            col = col[c(1,2,3,3)], lty=lty[c(1,2,3,3)],
            main = main, xlim = xlim, mar=mar, ...)
    box()
    axis(1, at=sqrt(xlb), labels=as.character(xlb))
    axis(3, at=sqrt(lambda), labels=rep('',length(lambda)))
    axis(2)
    par(usr=c(par('usr')[1:2], 0,1.05*max(l1)))
    lines(sqrt(lambda), l1, col=col[4], lty=lty[4])
    axis(4, col=col[4])
    mtext('penalty',4,3)
  } else warning(':plot.lassogrp: invalid argument "type". No plot')
}

## =========================================================
extract.lassogrp <-
  function(object, i=NULL, lambda=NULL, fitfun='lm', ...)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## does not make sense yet for more than one i
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 21 Aug 2009, 08:55
  lform <- object$formula
  lcall <- object$call
  if (is.null(lform)) lform <- lcall$x
  if (is.null(lform))
    stop('extract.lassogrp needs an object with formula')
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
  result <- NULL
  lmod <- lcall$model
  if (!is.null(lmod)) {
    if (is.character(lmod))
      lmod <- pmatch(lmod,c('gaussian','binomial','poisson'))
    else {
      lmod <- lmod@name
      lmod <- pmatch(substring(lmod,1,3),c('Linear','Logistic','Poisson'))
    }
  }
  lmeth <- c("lm","binomial","glm")[lmod]
  call <- if (fitfun=='regr')
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
  lreg <- eval(call, parent.frame())
  for (li in seq_along(i)) {
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
## ===================================================================
print.lassofit <-
  function(x, # dummycoef = NULL, #digits = max(3, options("digits")[[1]] - 3),
    residuals=FALSE, doc = options("doc")[[1]], ...)
{
  ## Author: Werner Stahel, Date: 24 Aug 2009, 12:27
  if (length(doc)==0) doc <- 0
  if (doc>=1) if (length(tit(x)))
    cat("\n",tit(x),"\n")
  if (doc>=2) if (length(descr(x)))
    cat(paste(descr(x),"\n"),"\n")
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
#  invisible(x)
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
    blocksinmodel <- blocks%in%all.vars(as.formula(object$call$formula))
  } else  {
    if (length(blocks)==1)
      blocklist <- split(sample(1:n), rep(1:blocks, length = n))
    else  {
      if (length(blocks)!=n)
        stop('!cv.lasso! argument "blocks" not suitable')
      blocklist <- split(sample(1:n), blocks)
    }
  }
  ## crossvalidate
  K <- length(blocklist)
  blockrmse <- matrix(0, K, length(object$lambda))
  fitted <- matrix(0, n, length(object$lambda))
  for (i in seq(K)) {
    omit <- blocklist[[i]]
    fit <-
      lassogrp.default(x, y, object$index, subset=-omit, weights=object$weights,
                       offset=object$offset, lambda=object$lambda,
                       penscale=object$penscale, center=object$innercall$center,
                       standardize=object$innercall$standardize,
                       save.x=FALSE, save.y=FALSE, control=control, 
                       model=object$model, adaptive=FALSE, cv.function=NULL,  
                       ...)
##-     fit <- update(object, subset=-omit, model=FALSE,
##-                   lambda=object$lambda, adaptive=FALSE, cv.function=NULL,
##-                   save.x=FALSE, weights=object$eights, offset=object$offset,
##-                   control=control)
    if (blocksinmodel)  {
      xnew <- mf[omit, , drop = FALSE]
      xnew[,blocks] <- factor(mf[-omit,blocks][1])
    } else {
      xnew <- x[omit, , drop = FALSE]
      fit$terms <- NULL
    }
    pred <- rbind(predict(fit, xnew))
    if (blocksinmodel) {
      pred <- sweep(pred, 2, apply(pred, 2, mean) - mean(y[omit]))
      blockrmse[i,] <- apply(y[omit] - pred, 2, sd)
    } else
      blockrmse[i,] <- sqrt(apply((y[omit] - pred)^2, 2, mean))
    fitted[omit,] <- pred
    if (trace)  cat("\n CV Fold", i, "\n\n")
  }
  ## summarize
  rmse <- apply(blockrmse, 2, mean)
  rmse.se <- sqrt(apply(blockrmse, 2, var)/K)
  result <- list(rmse = rmse, rmse.se = rmse.se, rmse.blocks = blockrmse,
                 fitted = fitted, lambda = object$lambda,
                 blocksinmodel = blocksinmodel)
  ## plot
  if (is.null(plot.it)) plot.it <- identical(parent.frame(), globalenv())
  if (plot.it)   plot.lassogrp(object, type='crit', cv = result, se = se)
  ##
  invisible(result)
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
  ## Author: Lukas Meier, Date: 30 Jun 2006, 09:46 (createDesign)

  if(!inherits(formula, "formula") || length(formula) != 3)
    stop("Argument 'formula' is of wrong type or length")
  any.nonpen <- !is.null(nonpen)
  ## Case where some variables won't be penalized -> merge formulas,
  ## also check the environments (is this the correct way ???)
  if(any.nonpen){
      if(!inherits(nonpen, "formula"))
        stop("Argument 'nonpen' of wrong type")

      is.gEnv <- function(e) identical(e, .GlobalEnv)

      ## Paste the two formulas together
      f <- as.formula(paste(c(deparse(formula[[2]]), "~",
                              deparse(formula[[3]]), "+",
                              deparse(nonpen[[length(nonpen)]])),
                            collapse = ""))

      ## Get the environment of the formulas
      env.formula <- environment(formula)
      env.nonpen <- environment(nonpen)

      ## If env. of 'formula' is not global, check if env. of 'nonpen' differs.
      ## If yes give warning and use env. of 'formula'
      if(!is.gEnv(env.formula)){
        environment(f) <- env.formula
        if(!is.gEnv(env.nonpen) && !identical(env.formula, env.nonpen))
          warning("'formula' and 'nonpen' have different environments. I use environment of 'formula'")
        
        environment(f) <- environment(formula)
      } else if (!is.gEnv(env.nonpen)){ ## if env. of 'nonpen' is not global
                 ## (but env. of 'formula' is), use env. of 'nonpen'
        environment(f) <- env.nonpen
      }
      m$formula <- f
  } else {
    m$formula <- formula
  }

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
    
    if(length(co))
      used.co <- co[!unlist(lapply(co, is.null))]
    else
      used.co <- NULL

    ## also uses the response...to be changed
    x.nonpen <- model.matrix(tmp, data = mf, contrasts = used.co)
    matches <- match(colnames(x.nonpen), colnames(x))
  }
  index <- attr(x, "assign")
  names(index) <- dimnames(x)[[2]]
  if(any.nonpen)
    index[matches] <- - index[matches]
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
varBlockScale <- function(x, ipen.which, inotpen.which, coef = NULL)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  4 Aug 2006, 08:50

  jcoef <- length(coef)>0
  if (jcoef&length(coef)!=ncol(x))
    stop('length(coef) not equal ncol(x)!')
  n <- nrow(x)
  x.ort <- x
  scale.pen <- vector(mode = 'list', length = length(ipen.which))
  scale.notpen <- NULL
  
  if(length(inotpen.which) > 0 & !jcoef){
    scale.notpen <- sqrt(apply(x[,inotpen.which]^2, 2, mean))
    x.ort[,inotpen.which] <- scale(x[,inotpen.which], FALSE, scale.notpen)
  }
##-     sc <- if (jcoef) coef[inotpen.which] else
##-       1/sqrt(apply(x[,inotpen.which]^2, 2, mean))
##-     x.ort[,inotpen.which] <- sweep(x[,inotpen.which], 2, sc, '*')
##-     scale.notpen <- sc

  ##!!! penalty for groups
  for(j in 1:length(ipen.which)){
    ind <- ipen.which[[j]]
    decomp <- qr(x[,ind])
    if(decomp$rank < length(ind)) ## Warn if block has not full rank
      stop("Block belonging to columns ", paste(ind, collapse = ", "),
              " has not full rank! \n")
    scale.pen[[j]] <- qr.R(decomp) * 1 / sqrt(n)
    x.ort[,ind] <- qr.Q(decomp) * sqrt(n)
  }
  list(x = x.ort, scale.pen = scale.pen, scale.notpen = scale.notpen)
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
  scale.notpen <- NULL
  
  if(length(inotpen.which) > 0){
    one <- rep(1, n)
    scale.notpen <- sqrt(drop(one %*% (x[,inotpen.which]^2)) / n)
    x.ort[,inotpen.which] <- scale(x[,inotpen.which], FALSE, scale.notpen)
  }
    
  for(j in 1:length(ipen.which)){
    ind <- ipen.which[[j]]
    decomp <- qr(x[,ind])
    if(decomp$rank < length(ind)) ## Warn if block has not full rank
      stop("Block belonging to columns ", paste(ind, collapse = ", "),
              " has not full rank! \n")
    scale.pen[[j]] <- qr.R(decomp) * 1 / sqrt(n)
    x.ort[,ind] <- qr.Q(decomp) * sqrt(n)
  }
  list(x = x.ort, scale.pen = scale.pen, scale.notpen = scale.notpen)
}

## ==================================================================
lambdamax <- function(x, ...)
  UseMethod("lambdamax")

lambdamax.formula <-
  function(formula, nonpen  = ~ 1, data, weights, subset, na.action,
           coef.init, penscale = sqrt, model = LogReg(),
           standardize = TRUE, contrasts = NULL, nlminb.opt = list(), ...)
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
  m$nonpen <- m$lambda <- m$coef.init <- m$penscale <- m$model <-
    m$standardize <- m$contrasts <- m$... <- NULL

  l <- lassogrpModelmatrix(m, formula, nonpen, data, weights, subset, na.action,
                     contrasts, parent.frame())

  if(missing(coef.init))
    coef.init <- rep(0, ncol(l$x))
  
  lambdamax.default(l$x, y = l$y, index = l$index, weights = l$w,
                    offset = l$off, coef.init = coef.init,
                    penscale = penscale,
                    model = model, standardize = standardize,
                    nlminb.opt = nlminb.opt, ...)
}

## ==================================================================
lambdamax.default <-
  function(x, y, index, weights = NULL, offset = rep(0, length(y)),
           coef.init = rep(0, ncol(x)), penscale = sqrt, model = LogReg(),
           standardize = TRUE, nlminb.opt = list(), ...)
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

  any.notpen <- any(index<=0)
  coef.npen <- coef.init[index<=0] ## unpenalized parameters
  
  inotpen.which <- which(index<=0)
  
  ## Index vector of the penalized parameter groups
  ipen <- index[index>0]
  
  ## Table of degrees of freedom
  dict.pen <- sort(unique(ipen))
  ipen.tab <- table(ipen)[as.character(dict.pen)] 

  ## Indices of parameter groups
  ipen.which <- split((1:ncol(x))[index>0], ipen)

  if(standardize){
    stand        <- varBlockStand(x, ipen.which, inotpen.which)
    x            <- stand$x
  }

  x.npen <- x[,inotpen.which, drop = FALSE]
  if (length(weights)==0) weights <- rep(1, length(y))

  helper <- function(par)
    model@nloglik(y, offset + x.npen %*% par, weights, ...)

  if(any.notpen){
    par0 <- do.call(nlminb, args = c(list(start = coef.npen,
                              objective = helper), nlminb.opt))$par
    mu0  <- model@invlink(offset + x.npen %*% par0)
  }else{
    mu0 <- model@invlink(offset)
  }
  
  ngrad0 <- model@ngradient(x, y, mu0, weights, ...)[index>0]
  
  ##gradnorms <- numeric(length(dict.pen))

  gradnorms <- c(sqrt(rowsum(ngrad0^2, group = ipen))) / penscale(ipen.tab)
                                                       
  ##for(j in 1:length(dict.pen)){
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
           save.y       = "logical",
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
           save.y       = TRUE,
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
           if(ceiling(object@update.every) != floor(object@update.every) |
              object@update.every <= 0)
             return("update.every has to be a natural number")

           if(ceiling(object@inner.loops) != floor(object@inner.loops) |
              object@inner.loops < 0)
             return("inner.loops has to be a natural number or 0")

           if(ceiling(object@max.iter) != floor(object@max.iter) |
              object@max.iter <= 0)
             return("inner.loops has to be a natural number or greater than 0")

           if(object@beta <= 0 | object@beta >= 1)
             return("beta has to be in (0, 1)")
           
           if(object@sigma <= 0 | object@sigma >= 1)
             return("sigma has to be in (0, 1)")

           if(object@tol <= 0)
             return("tol has to be positive")

           if(object@lower > object@upper)
             return("lower <= upper has to hold")

           return(TRUE)
         }
)

## ==================================================================
lassoControl <- function(save.x = FALSE, save.y = TRUE,
                         update.hess = c("lambda", "always"),
                         update.every = 3, inner.loops = 10,
                         line.search = TRUE, max.iter = 500,
                         tol = 5 * 10^-8, lower = 10^-2, upper = Inf,
                         beta = 0.5, sigma = 0.1, trace = 0){
  
  ## Purpose: Options for the Group Lasso Algorithm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## save.x: a logical indicating whether the design matrix should be saved.
  ## save.y: a logical indicating whether the response should be saved.
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
             save.y       = save.y,
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
## ===========================================================================
## ===========================================================================
descr <- function(x) attr(x,"descr")
## ---
"descr<-" <- function(x, value)
{
  ##-- Create descr attribute or  PREpend  new descr to existing one.
  value <- as.character(value)
  attr(x, "descr") <- if (length(value)==0) NULL else
  if(value[1]=="^") value[-1] else c(value, attr(x, "descr"))
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
