lassoselect <- function(object, data, s=NA, cv=is.na(s), cv.k=10,
                        cv.s=seq(0,1,length=100),
                        adaptive=TRUE, coefficients=object$coefficients, ...)
{
  ## Purpose:   model selection by lasso
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   object     either a formula or a result of an earlier call to
  ##              lassoselect
  ##   data       data frame in which the formula is evaluated
  ##   s          importance of l1 penalty given as in predict.lars
  ##              (mode="fraction")
  ##              if NA, s is determined by cross validation
  ##   cv         calculate cross validated ...
  ##              (see  cv.lars )
  ##   adaptive   if TRUE, the l1 penalty is sum(abs(beta/coefficients))
  ##              were  coefficients  is given as the next argument or
  ##              as  object$coefficients
  ##   coefficients  named vector of initial coeffficients used as
  ##              inverse weights in the l1 penalty.
  ##              Only needed if adaptive is TRUE.
  ##              Defaults to  object$coefficients
  ##   plot.it    passed to  cv.lars
  ## ----------------------------------------------------------------------
  ##   requires   library(lars)
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 29 Nov 2007, 08:22
  require(lars)
  ljlars <- TRUE
  if (is.formula(object)) {
    lfo <- object
    ldt <- data
    lmf <- model.frame(lfo, ldt)
  } else  {
    if (!inherits(object,"lassoselect")) stop("unsuitable first argument")
    lfo <- object$fullformula
    ldt <- object$data
    lcf0 <- object$coefficients 
    if (!is.na(s)) {
      lla <- object$lars
      lcv <- object$lcv
      lfo <- object$fullformula
      ljlars <- FALSE
    }
  }
  lmm <- model.matrix(lfo, ldt)[,-1]
  ly <- ldt[dimnames(lmm)[[1]],1]
  if (ljlars) {
    if (adaptive) {
      lcf0 <- coefficients
      if (length(lcf0)==0) { # stop("no initial coefficients found")
        lr <- lassoselect(object, data, adaptive=FALSE)
        lcf0 <- lr$coefficients
      }
      if (length(names(lcf0))==0) stop("initial coefficients must have names")
      if (any(lj <- !names(lcf0)%in%dimnames(lmm)[[2]]))
        stop(paste("term(s)", names(lcf0)[lj]," not in model matrix"))
      lmm <- sweep(lmm[,names(lcf0)],2,lcf0,"*")
    }
    lla <- lars(lmm,ly, normalize=!adaptive)  ## drop intercept column
    lcv <- NULL
    if (!is.na(s)&&(s<=0|s>1)) {
      warning("s value not suitable. I use minimizer of cross validation")
      s <- NA
    }
    if (is.na(s)|cv) {
      lcv <- data.frame( cv.lars(lmm,ly, K=cv.k, fraction=cv.s, plot.it=FALSE) )
      if (is.na(s)) s <- lcv$fraction[which.min(lcv$cv)]
    }
  }
  lcf <- predict.lars(lla, type="coefficients", mode="fraction",
                      s=s)$coefficients
  lcf <- lcf[lcf!=0]  ##  select variables
  lfit <- c(lmm[,names(lcf)]%*%lcf)
  if (adaptive) lcf <- lcf*lcf0[names(lcf)]
  lfonew <- update(lfo, if (length(lcf)>0)
    paste("~",paste(names(lcf), collapse=" + ")) else ~1) 
  lfit <- lfit+mean(ly)-mean(lfit)
  names(lfit) <- dimnames(lmm)[[1]]
  result <- list(lars=lla, data=ldt, formula=lfonew, fullformula=lfo,
                 s=s, coefficients=lcf,
                 fitted=lfit, adaptive=adaptive, cv=lcv)
  class(result) <- "lassoselect"
  result
}
## ===========================================================================
plot.lassoselect <-
  function(x, plot.s=TRUE, plot.cv=TRUE, plot.cvse=TRUE, main="", ...)
{
  ## Purpose:   plot coefficients and possibly cross validation trace
  ##            for a lassoselect object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 30 Nov 2007, 14:06
  lcl <- !(x$adaptive|plot.cv)
  lmain <- if (length(main))  
  i.plotlars(x$lars, plot.coeflabels=lcl, main=main, ...) # modification 
  if (plot.s>0) {
    lss <- x$s
    abline(v=lss,lty=6,col=plot.s)
    mtext(format(lss),3,0.2,at=lss,col=plot.s,cex=par("cex"))
  }
  if (plot.cv) {
    lcv <- x$cv
    if (!is.null(lcv)) {
      par(usr=c(par("usr")[1:2],0,max(lcv$cv+lcv$cv.error)))
      lines(lcv$fraction,lcv$cv,lty=3,col=plot.cv)
      if (plot.cvse) segments(lcv$fraction,lcv$cv+lcv$cv.error,
                             lcv$fraction,lcv$cv-lcv$cv.error)
      axis(4)
      mtext("MSPE",4,par("mgp")[2],adj=1.1)
    }}
}
## ===========================================================================
i.plotlars <- 
function (x, xvar = c("norm", "df", "arc.length", "step"), breaks = TRUE, 
    plottype = c("coefficients", "Cp"), omit.zeros = TRUE, eps = 1e-10, 
    plot.coeflabels=TRUE, ...) 
{
    object <- x
    plottype <- match.arg(plottype)
    xvar <- match.arg(xvar)
    coef1 <- object$beta
    if (x$type != "LASSO" && xvar == "norm") 
        coef1 = betabreaker(x)
    stepid = trunc(as.numeric(dimnames(coef1)[[1]]))
    coef1 <- scale(coef1, FALSE, 1/object$normx)
    if (omit.zeros) {
        c1 <- drop(rep(1, nrow(coef1)) %*% abs(coef1))
        nonzeros <- c1 > eps
        cnums <- seq(nonzeros)[nonzeros]
        coef1 <- coef1[, nonzeros, drop = FALSE]
    }
    else cnums <- seq(ncol(coef1))
    s1 <- switch(xvar, norm = {
        s1 <- apply(abs(coef1), 1, sum)
        s1/max(s1)
    }, df = object$df, arc.length = cumsum(c(0, object$arc.length)), 
        step = seq(nrow(coef1)) - 1)
    xname <- switch(xvar, norm = "|beta|/max|beta|", df = "Df", 
        arc.length = "Arc Length", step = "Step")
    if (plottype == "Cp") {
        Cp <- object$Cp
        plot(s1, Cp, type = "b", xlab = xname, main = object$type, 
            ...)
    }
    else {
        matplot(s1, coef1, xlab = xname, ..., type = "b", pch = "*", 
            ylab = "Standardized Coefficients")
        title(object$type, line = 2.5)
        abline(h = 0, lty = 3)
        if (plot.coeflabels)
        axis(4, at = coef1[nrow(coef1), ], label = paste(cnums), 
            cex = 0.8, adj = 0)
        if (breaks) {
            axis(3, at = s1, labels = paste(stepid), cex = 0.8)
            abline(v = s1)
        }
    }
    invisible()
}
## ===========================================================================
