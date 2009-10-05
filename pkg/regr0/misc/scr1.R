plot.regr <-
function(x, data=NULL, markprop=NULL, lab=NULL, cex.lab=0.7, 
         mf = NULL, mfcol=FALSE, mar=c(3,3,2,1), mgp=c(2,0.7,0),
         oma = 2*(prod(mf)>1), cex=par("cex"), ask = NULL,
         multnrows = 0, multncols=0, 
         lty = c(1,2,5,3,4,6), lwd = c(1,1,2,1,1.5,1),
         colors = options("colors.ra")[[1]], pch=NULL, 
         main = NULL, cex.title = NULL, wsymbols=NULL, symbol.size=NULL, 
         smooth = TRUE, smooth.par=NA, smooth.iter=NA, 
         smooth.sim=19, nxsmooth=51, 
         plotselect = NULL, 
         weights = NULL, hat.cooklim = 1:2, res.lim = TRUE, y.lim=TRUE,
         glm.restype = "deviance", sequence=NA, 
	 xplot = TRUE, x.se=FALSE, x.smooth = smooth, addcomp=FALSE,
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
##  res.lim    plotting limits for residuals
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
  lpolr <- inherits(x,"polr")
  lfcount <- lfam == "binomial" | lfam == "poisson" | lfam == "multinomial" |
             lpolr
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
  if(is.na(pmatch(rtype,"cond.quant")))
    lres <- residuals(x, type=rtype) else {
      lres <- residuals.polr(x)
      lpolr <- NCOL(lres)>1
    }
  if (NCOL(lres)==1) lres <- cbind(lres)  ## do not destroy class attr
  if (var(lres[,1],na.rm=TRUE)==0)
    stop("!plot.regr! all residuals are equal -> no residual plots")
  lmres <- if (lpolr) 1 else ncol(lres)
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
  if (is.logical(res.lim)&&!is.na(res.lim)) res.lim <- 
    if (res.lim) {
      if(lpolr) robrange(c(lresa)) else apply(lresa,2,robrange)
    } else  NULL
  if (!is.null(res.lim))
    if (any(dim(cbind(res.lim))!=c(2,lmres))) {
      warning("!plot.regr! unsuitable argument  res.lim ")
      res.lim <- NULL
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
  if(length(lhat)==0) lhat <- hat(x$qr) # [lnna]
##-   if (lwgt) lhat <- lhat/lwgts use almost-unweighted h for plotting
## standardized residuals
  lstres <- x$stres
  lstres <- if (length(lstres)>0)
    cbind(lstres)[lnna,,drop=FALSE] else {
      if (all(lsigma>0)) {
        if (any(is.na(lhat))) sweep(lresa,2,lsigma,"/") else
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
## smooth
  if (is.logical(smooth)) smooth <- if (smooth) 
    function(x,y,weights,par,iter)
          loess(y~x, weights=weights, span=par, iter=iter)$fitted
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
      order(abs(lstres[,1]))[1:(lna*(1-markprop))]  # [,1] for lpolr
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
    c("black","gray","blue","cyan","red","magenta","darkgreen","gray30")  else
    rep(colors,length=8)
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
        warning("!plot.regr! I cannot plot categorical Y on fit")
      else {
      ly <- x[["y"]]
      if (length(ly)==0) ly <- lf + lresa
      lsimr <- c(lf)+lsimres
      i.plotlws(lf,ly, lfname,lyname,main,outer.margin,cex.title,
                colors,lty,lwd,lpch, 
                lpty, llabna,cex.lab,ltxt, lwsymbols,lwgt,lwgts,liwgt,lsyinches,
                lpllevel>1,smooth,lsmpar,smooth.iter, ylim=y.lim,
                reflinex=0,refliney=1,lnsims=lnsims, simres=lsimr)
  } } }
## ---
  if(lpls=="ta") {
    if (is.na(lpllevel)) lpllevel <- 3*(lplsel["yfit"]==0)
    if (is.na(lpllevel)) lpllevel <- 3
    if (lpllevel>0) 
      i.plotlws(lf,lresa, lfname,lrname,main,outer.margin,cex.title,
                colors,lty,lwd,lpch, lpty, llabna,cex.lab,ltxt, 
                lwsymbols,lwgt,lwgts,liwgt,lsyinches,
                lpllevel>1,smooth,lsmpar,smooth.iter, ylim=res.lim,
                reflinex=mean(lf),refliney=-1,
                lnsims=lnsims, simres=lsimres)
  }
## --- 
  if(lpls=="tascale")  {
    if (is.na(lpllevel)) lpllevel <- 3*(!lfcount)
    if (lpllevel>0) 
      i.plotlws(lf,lrabs, lfname, paste("|",lstrname,"|"),
        main, outer.margin, cex.title,
        colors,lty,lwd,lpch, lpty, llabna,cex.lab,ltxt,
        FALSE,FALSE, # lplwgt=F,wgt=F,
        lwgts,liwgt,lsyinches, lpllevel>1, smooth,lsmpar,smooth.iter, 
        smooth.power=0.5,
        ylim=c(0,1.05*max(lrabs,na.rm=TRUE)),yaxs="i")
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
      warning("!plot.regr! no weights foound. cannot plot absres on weights")
    else {
      i.plotlws(lweights,lrabs, "weight (log)",
        paste("|",lstrname,"|"),main,outer.margin, cex.title,
        colors,lty,lwd,lpch, lpty, llab,cex.lab,ltxt, FALSE,FALSE, # lplwgt=F,wgt=F,
        lwgts,liwgt,lsyinches,
        lpllevel>1,smooth,lsmpar,smooth.iter, smooth.power=0.5,
        log="x", ylim=c(0,1.05*max(lrabs,na.rm=TRUE)),yaxs="i")
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
        li <- order(lhat, decreasing=TRUE)[1:(lna*(1-markprop/2))]
        llabh[li] <- llabels[lnna][li]
      }
      lhattit <- paste("hat diagonal", if(lwgt) "(unweighted)")
      for (lj in 1:lmres) {
        i.plotlws(lhat, lstres[,lj], lhattit, lstrname[lj],
              main,outer.margin,cex.title,
              colors,lty,lwd,lpch, lpty, llabh,cex.lab,ltxt,
              lwsymbols,lwgt, lwgts,liwgt,lsyinches,
              FALSE, ylim=res.lim) 
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
##-               colors,lty,lpch, 
##-             lpty, llabna,cex.lab,ltxt, lwsymbols,lwgt,lwgts,liwgt,lsyinches,
##-             lpllevel>1,smooth,lsmpar,smooth.iter, ylim=res.lim)
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
  llabs <- rep(NA,ln)
#  llabs[lnna] <- llabna
  lylim <- if (addcomp) y.lim else res.lim
  if (lxpl | is.na(lseq)|lseq) {
    plresx(x, partial.resid=TRUE, glm.restype = glm.restype, lab=llabs,
           cex.lab=cex.lab,
           vars=xvars, sequence=lseq, se=x.se, addcomp = addcomp,
           smooth=x.smooth, smooth.par=lsmpar, 
           smooth.sim=lsimres, lty=lty, lwd=lwd, colors = colors, 
           main=main, cex.title=cex.title, ylim = lylim, 
           wsymbols=lwsymbols, symbol.size=symbol.size, pch=lpch, ...)
  }
## --- end
  par(loldpar)
  "plot.regr done"
}
