## relevance.R
## two groups and one sample
twosamples <- #f
  function(x, ...) UseMethod("twosamples")
twosamples.default <- #f
  function(x, y=NULL, paired=FALSE, var.equal=TRUE, testlevel=0.05, rlvThres=NULL, ...)
{ ## effect (group difference) and relevance
  lpvs <- rlvoptions("pSymbols")
  lrls <- rlvoptions("pSymbols")
  lonegroup <- FALSE
  if (length(y)==0) {
    if (NCOL(x)==2) {
      y <- x[,2]
      x <- x[,1]
    } else {
      if (!is.atomic(x))
        stop("!twosamples! argument 'x' not suitable")
##      else lonegroup <- TRUE
    }
  }
  if (paired) {
    if (length(x)!=length(y))
      stop("!twosamples! 'x' and 'y' must have equal length if 'paired' is TRUE")
    x <- y-x
    lonegroup <- TRUE
  }
  if (length(y)==0) lonegroup <- TRUE
  lpq <- 1-testlevel/2
  lrlvth <- i.def(i.def(rlvThres, i.getopt("rlvThres")["stand"]), 0.1)
  x <- x[is.finite(x)]
  ln <- lnx <- length(x)
  ldf <- lnx-1
  ld <- lmnx <- mean(x)
  lssx <- sum((x-lmnx)^2)
  lv <- lssx/ldf
  lny <- NULL
  method <- "One Sample t inference"
  if (!lonegroup) {
    y <- y[is.finite(y)]
    lny <- length(y)
    ln <- lnx+lny
    lmny <- mean(y)
    ld <- lmny-lmnx
    lssy <- sum((y-lmny)^2)
    if (var.equal) {
      lv <- ln*(1/lnx+1/lny)*(lssx+lssy)/(ln-2)
      ldf <- ln-2
      method <- "Two Sample t inference"
    } else { ## welch
      lv <- ln*(lssx/(lnx-1)/lnx+lssy/(lny-1)/lny)
      lsex2 <- lssx/(lnx-1)/lnx
      lsey2 <- lssy/(lny-1)/lny
      ldf <- (lsex2+lsey2)^2/(lsex2^2/(lnx - 1) + lsey2^2/(lny - 1))
      method <- "Two Sample t inference, unequal variances (Welch)"
    }
  }
  lse <- sqrt(lv/ln)
  lq <- qt(lpq, ldf)
  lciwid <- lq*lse
  leffci <- ld + c(0,-1,1)*lciwid
  lpv <- 2*pt(-abs(leffci[1]/lse), ldf)
  lrlvci <- leffci/sqrt(lv)/lrlvth
  structure(
    as.data.frame(
      list(effect=leffci[1], ciLow=leffci[2], ciUp=leffci[3],
           Rle=lrlvci[1], Rls=lrlvci[2], Rlp=lrlvci[3],
           Rls.symbol=getsymbol(lrlvci[2], lrls),
           Sig0=c(leffci[1]/lciwid), Sigth=c((lrlvci[1]-1)*2/diff(lrlvci[2:3])),
           p.value=lpv, p.symbol=getsymbol(lpv, lpvs))),
    class=c("effecttable", "data.frame"),
    rlvThres=structure(lrlvth, type="standardized"),
    n=c(lnx, lny), df=ldf, means=if (!lonegroup) c(lmnx,lmny), V=lv,
    pSymbols=rlvoptions("pSymbols"), rlvSymbols=rlvoptions("rlvSymbols"))
}
## ============================================================================
twosamples.formula <- 
  function (formula, data, subset, na.action, ...) 
{ ## adapted from t.test.formula
  if (missing(formula) || (length(formula) != 3L)) 
    stop("!twosamples! 'formula' must have left and right term")
  oneSampleOrPaired <- FALSE
  if (length(attr(terms(formula[-2L]), "term.labels")) != 1L) 
    if (formula[[3]] == 1L) 
      oneSampleOrPaired <- TRUE
    else stop("!twosamples! 'formula' incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) 
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(paste("`",names(mf),"'",sep=""), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  if (!oneSampleOrPaired) { ## two indep. samples
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2L) 
      stop("grouping factor must have exactly 2 levels")
    DATA <- setNames(split(mf[[response]], g), c("x", "y"))
    rr <- do.call("twosamples", c(DATA, paired=FALSE, list(...)))
    names(attr(rr, "n")) <- levels(g)
  } else { ## one sample
    respVar <- mf[[response]]
    if (inherits(respVar, "Pair")) {
      DATA <- list(x = respVar[, 1], y = respVar[, 2]) 
      rr <- do.call("twosamples", c(DATA, paired = TRUE, list(...)))
    } else {
      DATA <- list(x = respVar)
      rr <- do.call("twosamples", c(DATA, paired = FALSE, list(...)))
    }
  }
  rr$data.name <- DNAME
  rr
}
## ===========================================================================
## regression
ciSgRl <-
  function (coef=NULL, se=NULL, df=Inf, testlevel=0.05, stcoef=TRUE, rlv=TRUE,
            rlvThres=0.1, object=NULL)
{ ## for a coefficients table,
  ## calculate confidence interval, significance and relevance
  if (length(coef)==0) coef <- object
  if (inherits(coef, regrModelClasses)) {
    object <- coef
    coef <- getcoeftable(object)
    df <- df.residual(object)
  }
  if (!(is.atomic(coef)||length(dim(coef))==2))
    stop("!ciSgRl! first argument not suitable")
  lpvs <- rlvoptions("pSymbols")
  lrls <- rlvoptions("rlvSymbols")
  if (is.null(se))
    if (NCOL(coef)>1) {
      se <- coef[,2]
      coef <- coef[,1]
    }
  if (is.null(se)) se <- attr(coef, "se")
  if (is.null(se))
    stop("!ciSgRl! no standard errors found")
  if (df==0||!any(is.finite(se))) {
    warning("!ciSgRl! no finite standard errors")
    return(cbind(coef))
  }
  ltq <- qt(1-testlevel/2, df)
  lci <- coef+outer(ltq*se, c(ciLow=-1,ciUp=1))
  ltst <- coef/se
  lsgf <- ltst/ltq
  lpv <- 2*pt(-abs(ltst), df)
  lpsymb <- getsymbol(lpv, lpvs)
  rr <- data.frame(coef=coef, se=se, lci, testst=ltst,
                   Sig0=lsgf, p.value=lpv, p.symbol=lpsymb)
  attr(rr, "pSymbols") <- lpvs
  ## --- standardized coefficients
  if (!u.notfalse(stcoef)) return (rr)
  ##    warning(":ciSigRlv: no standardized coefficients given. No relevances")
  if (length(stcoef)==0||u.true(stcoef)) {
    if (length(object)==0) {
      warning(":ciSgRl: argument 'object' needed. No relevances")
      return (rr)
    }
    if (inherits(object, "nls")) {
      warning(":ciSgRl: cannot calculate standardized coefficients for 'nls' objects.")
      return(rr)
    }
    lfac <- getcoeffactor(object)
    ##   lcls <- attr(lfac, "fitclass")
    stcoef <- coef*lfac
  }
  if (length(stcoef)!=length(coef)) {
    warning(":ciSgRl: argument 'stcoef' not suitable. No relevances")
    return (rr)
  }
  lfac <- stcoef/coef
  lstci <- cbind(stcoef=stcoef, st.Low=lfac*lci[,1], st.Up=lfac*lci[,2])
  rr <- cbind(rr, lstci)
  ## --- relevance
  if (!rlv) return(rr)
  if (length(rlvThres)>1) rlvThres <- rlvThres["coef"]
  if (is.na(rlvThres)) {
    warning(":ciSgRl: argument 'rlvThres' not suitable. No relevances")
    return (rr)
  }
  lrlv <- lstci/rlvThres
  li <- which(lstci[,1]<0)
  if (length(li)) lrlv[li,] <- - lrlv[li,c(1,3,2)]
  colnames(lrlv) <- c("Rle","Rls","Rlp")
  lrlvsy <-
    if (any(!is.na(ll <- lrlv[,"Rls"]))) getsymbol(ll, lrls) else NA
  cbind(rr, lrlv, Rls.symbol=lrlvsy)
}
## -------------------------------------------------------------------------
getcoeftable <-
  function (object)
{ ## get coefficient table from model fit
  if (inherits(object, "regr")) ltb <- object$coeftable
  else {
    if (inherits(object, "survreg")) {
      ltb <- summary(object)$table
      ltb <- ltb[-nrow(ltb),]
    } else {
      if (inherits(object, "coxph")) {
        lcoef <- object$coefficients
        se <- sqrt(diag(object$var))
        ltb <- cbind(lcoef, se, lcoef/se,  ## exp(lcoef),
                     pchisq((lcoef/se)^2, 1, lower.tail = FALSE))
        dimnames(ltb) <-
          list(names(lcoef), c("coef", "se", "z", "p")) ## "exp(coef)",
      } else  ltb <- summary(object)$coefficients
    }
    if (inherits(object, "polr")) ltb <- ltb[names(object$coefficients),]
  }
  lnm <- names(coefficients(object))
  if (any(is.na(match(lnm,row.names(ltb))))) {
    rr <- structure(matrix(NA, length(lnm), ncol(ltb)), dimnames=list(lnm,colnames(ltb)))
    rr[row.names(ltb),] <- ltb
    rr
  } else ltb
}
## --------------------------------------------------------------------------
termtable <-
  function (object, summary=NULL, testtype=NULL, r2x = TRUE, rlv = TRUE,
            rlvThres=c(rel=0.1, coef=0.1, drop=0.1, pred=0.05), testlevel=0.05)
{
  ## Purpose:  generate term table for various models
  ## --------------------------------------------------------
  ## prep
  if (length(summary)==0) summary <- summary(object)
  ## thresholds
  ltlev1 <- 1-testlevel
  lrlthrl <- rlvThres["rel"]
  lrlthcf <- rlvThres["coef"]
  lrlthdr <- rlvThres["drop"]
  lrlthpr <- rlvThres["pred"]
  ## --- sigma and threshold for standardized coefficient
  lsigma <- c(object$sigma, summary$sigma, 1)[1]
  lfamily <- if (is.character(lfm <- object$family)) lfm else lfm$family
  ldist <- object$dist
  if ((inherits(object, "glm") &&
       lfamily %in% c("poisson","quasipoisson")) ||
      (inherits(object, "survreg")&&
       ldist %in% c("weibull", "exponential", "lognormal"))
      )  lrlthcf <- lrlthrl
  ## --- testtype
  if (length(testtype)==0) {
    testtype <- "LRT"
    if (inherits(object, c("lm","lmrob"))) testtype <- "F"
    if (inherits(object, "glm")) {
      testtype <- 
        if (lfamily %in% c("quasibinomial","quasipoisson")) "F" else "LRT"
    }
    if (inherits(object, c("survreg"))) testtype <- "Chisq"
  }
  ## ---
  lterms <- terms(object)
  if(length(attr(lterms,"term.labels"))==0)
    return(data.frame(
      coef=c(object$coefficients,NA)[1], df = NA, se=NA, ciLow=NA, ciUp=NA,
      Sig0=NA, stcoef=NA, stciLow = NA, stciUp = NA,
      testst=NA, p.value=NA, p.symbol="", R2.x=NA,
      stringsAsFactors=FALSE)
      )
  ## degrees of freedom
  ldfres <- df.residual(object)
  if (ldfres<1) {
    warning(":termtable: no degrees of freedom left.")
    return(data.frame(
      coef=c(object$coefficients,NA)[1], df = NA, se=NA, ciLow=NA, ciUp=NA,
      Sig0=NA, stcoef=NA, stciLow = NA, stciUp = NA,
      testst=NA, p.value=NA, p.symbol="", R2.x=NA,
      stringsAsFactors=FALSE)
      )
  }
  ## --- coefficients
  lcoef <- object$coefficients
##-   lcoeftab <- object$coeftable
##-   if (length(lcoeftab)==0||ncol(lcoeftab)<10)
  lcoeftab <- ciSgRl(object=object)
  names(lcoeftab)[1] <- "coef"
  ## --- drop1
  ldr1 <-
    if (inherits(object, c("lm","lmrob"))&&!inherits(object, "glm")) {
      lcov <- summary$cov.unscaled
##-       if (u.debug())
##-         drop1Wald(object, test=testtype, scope=lterms) ## scale for AIC!!!
      ##-       else
      try(drop1Wald(object, test=testtype, scope=lterms), silent=TRUE)
    } else {
##-       if (u.debug())
##-         drop1(object, test=testtype, scope=lterms)
      ##-       else
      try(drop1(object, test=testtype, scope=lterms), silent=TRUE)
    }
  if (inherits(ldr1, "try-error")) {
    warning(":termtable: drop1 did not work. I return the codfficient table")
##                  produced by ", object$fitfun
    return(list(termtable=lcoeftab))
  }
  ldr1 <- ldr1[-1,]
  ldr1$RSS <- NULL # same ncol for lm and glm
  ldf <- ldr1[,1]
  if (inherits(object,"rlm"))  ldr1[,4] <- ldr1[,2]/ldf
  ## -- critical value for test
  ltstq <- if (testtype=="F") qf(ltlev1,c(1,pmax(1,ldf)),ldfres) else {
    if (testtype%in%c("Chisq", "LRT")) qchisq(ltlev1,c(1,pmax(1,ldf))) else NA }
  ltstq[which(ldf==0)+1] <- NA
  ltstq1 <- sqrt(ltstq[1]) ## 1 degree of freedom
  ltstq <- ltstq[-1]
##-   if (inherits(object,"mlm")||inherits(object,"manova"))
  ##-     return(list(termtable=ldr1))  ## !!! needs much more
  ## ---------------------
  ## model.matrix
  lmmt <- object[["x"]]
  if (length(lmmt)==0)  lmmt <- model.matrix(object)
  lasg <- attr(lmmt,"assign")[!is.na(lcoef)]
##  if (class(object)[1] %in% c("polr")) lasg <- lasg[-1] ## ,"coxph"
  ## terms without factor involvement
  lfactors <- attr(lterms,"factors")
  lvcont <- !attr(lterms,"dataClasses")[row.names(lfactors)] %in%
    c("numeric","logical") ## [...] excludes .weights. and possibly others
  ## terms only containing continuous variables
  lcont <- which( lvcont %*% lfactors ==0 )
  ## -----------------------------------
  ## --- r2x
  lr2 <- NA
  if (r2x) {
    lvift <-     ## lterms: n of levels for each term
##        if (u.debug()) vif.regr(object, mmat=lmmt) else
        try(vif.regr(object, mmat=lmmt), silent=TRUE)
    if (inherits(lvift, "try-error") || length(lvift)==0) {
      warning(":termtable: error in the calculation of R2.x")
      lvif <- NA
    } else lvif <- lvift[,3]^2
    lr2 <- 1-1/lvif
  }
  ## -----------------------------------
  ## --- prepare table
  lpvcol <- pmatch("Pr(",names(ldr1), nomatch=ncol(ldr1))
  lpv <- ldr1[,lpvcol]
  ldf <- ldr1[,1]
  ## !!!
  lcont <- which(ldf==1)
  lnobs1 <- ldfres+sum(ldf)-1
  ## drop effect relevance
  ltst <- ldr1[, ifelse(inherits(object, "polr"), 3, 4)]
  ldrncci <- rbind(confintF(ltst, ldf, ldfres, testlevel))
  ldreff2 <- cbind(ltst*ldf, ldrncci)/lnobs1
  ldrrl <- cbind(NA,NA,NA, sqrt(ldreff2)/lrlthdr,
                 pmax(0.5*log((ldfres+ldreff2*lnobs1)/(ldfres+ldf))/lrlthpr, 0))
  dimnames(ldrrl) <-
    list(NULL, c(t(outer(c("coef","drop","pred"),c("Rle","Rls","Rlp"),paste, sep=""))))
  ltst <- ldr1[,lpvcol-1]
  ## table, filled partially
  ltb <- data.frame(coef=NA, ldr1[,1,drop=FALSE], ## keep name if only 1 coef
                    se=NA, ciLow=NA, ciUp=NA,
                    Sig0=sqrt(pmax(0,ltst)/ltstq),
                    stcoef=NA, st.Low=NA, st.Up=NA,
                    testst=ltst, p.value=lpv, p.symbol="", R2.x=lr2, ldrrl,
                    stringsAsFactors=FALSE)
  names(ltb)[2] <- "df"
  ## intercept
  ljint <- "(Intercept)"==names(lcoef)[1]
  if (ljint) {
    ##    ltstint <- # if(class(object)[1] %in% c("lm","nls","rlm"))
    ##      lcoeftab[1,3]^2 # else lcoeftab[1,3]
    ltb <- rbind("(Intercept)"=ltb[1,],ltb)
    ltb[1,] <- NA
    ltb[1,"df"] <- 1
    ltb[1,"p.symbol"] <- ""
    ltstq <- c(ltstq1, ltstq)
    lcont <- c(0, lcont)
  } else if (!inherits(object, c("coxph", "polr")))
    warning(":termtable: No intercept. Statistics are difficult to interpret.")
  lcont1 <- lcont+ljint  # row number in dr1
  ## --- coefficients and statistics for terms with 1 df
  if (length(lcont)) { ## lcont refers to assign
    ltlb <- dimnames(ltb)[[1]]
    lclb <- ltlb[lcont1] ## lcont1 is the row in the coef table of summary(object)
    ljc <- match(lcont,lasg) # index of coefs for cont variables
    lj <- c("coef","se","ciLow","ciUp","stcoef", "st.Low","st.Up",
            "coefRle","coefRls","coefRlp")
    ljj <- c("coef","se","ciLow","ciUp","stcoef", "st.Low","st.Up",
            "Rle","Rls","Rlp")
    ltb[lcont1,lj] <- lcoeftab[ljc,ljj]
    ltb[lcont1,"Sig0"] <- sign(ltb[lcont1,"coef"])*ltb[lcont1,"Sig0"]
  }
  if (row.names(lcoeftab)[nrow(lcoeftab)]=="Log(scale)") { # survreg
    ltsc <- lcoeftab[nrow(lcoeftab),]
    lcont1 <- c(lcont1, nrow(lcoeftab))
    if (!u.true(object$dist=="weibull")) ltsc[2:4] <- NA
    lls <- ltb[1,]
    lls[1,] <- NA
    lq <- qnorm(1-testlevel/2)
    lls[1,1:7] <-
      c(ltsc[1],ltsc[2],ltsc[1]+c(-1,1)*lq*ltsc[2], 1, ltsc[3], ltsc[3]/lq)
    lls[,"p.value"] <- ltsc[4]
    ltb <- rbind(ltb,"log(scale)"= lls)
  }
  ## p-symbol
  lpvs <- rlvoptions("pSymbols")
  ltb[,"p.symbol"] <-
    if (any(!is.na(ll <- ltb$p.value))) getsymbol(ll, lpvs) else NA
  ## Rls-symbol
  lrls <- rlvoptions("rlvSymbols")
##-   ltb[,"Rls.symbol"] <-
##-     if (any(!is.na(ll <- ltb$Rls)))
  ##-       lrls$symbol[as.numeric(cut(ll, lrls$cutpoint))] else NA
  for (lr in c("coefRls","dropRls","predRls"))
    ltb[,paste(lr,"symbol",sep=".")] <-
      if (any(!is.na(ll <- ltb[,lr]))) getsymbol(ll, lrls) else NA
  ## ---
  structure(ltb, class=c("termtable", "data.frame"),
            testtype=testtype,
            fitclass=class(object), family=lfamily, dist=ldist,
            rlvThres=rlvThres,
            pSymbols=lpvs, rlvSymbols=lrls)
}
## end termtable
## ===========================================================
confintF <- function(f, df1, df2=Inf, testlevel=0.05) {
  ## confidence interval for non-centrality of F distribution
  p <- testlevel/2
  lf.fq <- function(x, fvalue, df1, df2, p) qf(p,df1,df2,x)-fvalue
  lf.ciup <- function(fvalue, df, p) { ## upper bound for upper limit
    lq <- 1.5*qnorm(p)
    lu <- lq^2*2/df
    df*(fvalue-1+lu+sqrt(lu*(lu+2*fvalue-1)))
  }
  ln <- max(length(f), length(df1), length(df2), length(p))
  f <- rep(f, length=ln)
  df1 <- rep(df1, length=ln)
  df2 <- rep(df2, length=ln)
  p <- rep(p, length=ln)
  p <- pmin(p,1-p)
  ## ---------------------------
  rr <- matrix(NA, ln, 2)
  for (li in 1:ln) {
    lx <- f[li]
    if (!is.finite(lx)) next
    if (lx>100)
      rr[li,] <- df1[li]*(sqrt(lx)+c(-1,1)*abs(qt(p[li],df2[li]))/sqrt(df1[li]))^2
    else {
      rr[li,1] <-  ## lower limit
        if (lf.fq(0, f[li], df1[li], df2[li], 1-p[li])>=0) 0
        else
          uniroot(lf.fq, c(0,df1[li]*f[li]),
                  fvalue=f[li], df1=df1[li], df2=df2[li], p=1-p[li])$root
      rr[li,2] <- ## upper limit
        if (pf(f[li], df1[li], df2[li])<=p[li]) 0  ## tiny F value
        else
        uniroot(lf.fq, interval=c(df1[li]*f[li], lf.ciup(f[li], df1[li], 1-p[li])),
                fvalue=f[li], df1=df1[li], df2=df2[li], p=p[li], extendInt="upX")$root
    }
  }
  if (ln==1) c(rr) else rr
}
## --------------------------------------------------------------------
print.termtable <-
  function(x, columns=NULL, printstyle = NULL, legend = NULL,
           digits = NULL, na.print = "  ", ...)
{
  digits <- i.getopt(digits)
  lcoldef <- rlvoptions(
    if (i.getopt(printstyle)=="relevance")
      "termcolumns.r" else "termcolumns.o")
  columns <- i.def(i.def(columns, rlvoptions("termcolumns")), lcoldef) 
  columns[columns=="estimate"] <- "coef"
  names(x)[names(x)=="estimate"] <- "coef"
  if (length(columns)==1 && columns=="all")
    columns <- names(x)
  if (length(ll <- setdiff(columns, names(x)))) {
    warning(":print.termtable: columns not found:  ", paste(ll, collapse=", "))
    if (length(lcol <- setdiff(columns, ll))) columns <- lcol else return()
  }
  lattr <- attributes(x)
## ---
  x <- data.frame(x)[,columns, drop=FALSE]
  lcnames <- colnames(x)
  ## --- round some columns to 3 digits
  ljrp <- lcnames[pmatch(c("R2","Sig0","p.v"), colnames(x), nomatch=0)]
  if (length(ljrp))
    x[,ljrp] <- round(as.matrix(x[,ljrp]),max(3,digits))
  if ("Sig0" %in% ljrp) x$Sig0 <- round(x$Sig0,i.last(digits)-1)
  ## --- paste symbols to numbers
  ljsy <- grep(".symbol", lcnames)
  if (any(ll <- lcnames[ljsy]=="p.symbol"))
    lcnames[ljsy[ll]] <-
      if ("p.value" %in% lcnames) "p.value.symbol"
      else {
        if ("Sig0" %in% lcnames) "Sig0.symbol" else "coef.symbol"
      }
  if (length(ljsy)) {
    ljc <- match(sub(".symbol", "", lcnames[ljsy]),lcnames)
    if (any(ll <- is.na(ljc))) {
      if (sum(ll)==1) ljc[ll] <- 1
      else {
        ljc <- ljc[!ll]
        ljsy <- ljsy[!ll]
      }
    }
    if (lnc <- length(ljc)) {
      x[,ljc] <-
        matrix(paste(sub("NA",na.print, format(as.matrix(x[,ljc]))),
                     sub("NA","  ",format(as.matrix(x[,ljsy], justify="left")))),ncol=lnc)
      x <- x[,-ljsy, drop=FALSE]
    }
  }
  ## -----------------------------
  xf <- format(x, na.encode=FALSE)
  xp <- data.frame(lapply(xf, function(x) sub("NA",na.print,x)),
                   row.names=row.names(x))
  print(if(ncol(xp)==1) t(xp) else xp, quote=FALSE, na.print=na.print, ...)
  ## --- legend(s)
  lprleg <- i.def(legend, rlvoptions("show.symbolLegend"))
  if (lprleg) {
    if (length(grep(".symbol", columns))) {
      ## cat("---")
      if ("p.symbol" %in% columns)
        cat("\nSignificance codes for p.value:  ", lattr[["pSymbols"]]$legend,"\n", sep = "")
      if (length(grep("rlvSymbol", columns)))
        cat("\nRelevance codes:    ", lattr[["rlvSymbols"]]$legend,"\n", sep = "")
      cat("\n")
    }
  }
  invisible(xp)
}
## ===========================================================================
termeffects <-
  function (object, se = 2, df = df.residual(object), rlv = TRUE)
    ## --------------------------------------------------------------
{
  if (is.atomic(object)||is.null(terms(object)))
      stop("!termeffects! inadequate first argument")
 ##  xl <- object$xlevels
  Terms <- delete.response(terms(object))
  tl <- attr(Terms, "term.labels")
  dcl <- attr(Terms,"dataClasses")[-1]
  if (all(dcl=="numeric")) ## ??? need coeftable. is it always available?
    return(as.list(coef(object))) ## !!! coeftable!
  ## result already available?
  allc <- object$termeffects
  if ((!is.null(allc))&&length(allc)==length(tl)&&
      (is.matrix(allc[[length(allc)]])|!se)) return(allc) ## !!! check!
  ## ---
  if (rlv) lsigma <- getscalepar(object)
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
    for (li in imat) lctr <- c(lctr, list(diag(length(xtlv[[li]]))))
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
  names(asgn) <- colnames(mm)
##-   lnna <- is.finite(coef)
##-   if (any(!lnna)){
##-     coef <- coef[lnna]
##-     mm <- mm[,lnna]
##-     asgn <- asgn[lnna]
##-   }
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
##-   licf <- pmatch(colnames(mm), names(coef))
  licf <- pmatch(names(coef), colnames(mm))
  lasgn <- asgn[licf]
  mm <- mm[,licf]
##  asgn <- asgn[names(coef)] ## !!!
  res <- setNames(vector("list", length(tl)), tl)
  ljfail <- NULL
  for (j in seq_along(tl)) {
    mmr <- rn == tl[j]  ## rows corresponding to the term
    mmc <- lasgn==j ## & !is.na(coef)
##-     lcf <- coef[licf[mmc]]
    lcf <- coef[mmc]
##    mmc <- names(asgn)[asgn == j & !is.na(coef)]  ## columns (logical fails for polr, vcov() too large) !!! was  which
    ##-     mmpart <- mm[mmr, mmc, drop=FALSE]
    if (all(is.finite(lcf))) {
      mmpart <- mm[mmr,mmc, drop=FALSE]
      rrj <- setNames(drop(mmpart %*% lcf), rnn[mmr]) ## coef[mmc]
      if (se) {
        sej <- sqrt(diag(mmpart %*% cov[mmc,mmc] %*% t(mmpart)))
        if (any(is.na(rrj))|any(!is.finite(sej))) {
          ##-         warning(":termeffects: missing coef or non-finite standard error for term '",
##-                 tl[j], "'. no standard errors etc")
          ljfail <- c(ljfail, tl[j])
        } else {
          rrj <- ciSgRl(rrj, sej, df, stcoef=rrj*0.5/lsigma, rlv=rlv)
        }
      }
      res[[j]] <- rrj
    }
  }
  if (length(ljfail))
    warning(":termeffects: error calculating se for terms  ",
            paste(ljfail, collapse=", "))
  if (int > 0) {
    res <- c(list(`(Intercept)` = coef[int]), res)
  }
##-   if (lIpolr) {
##-     lcfi <- object$intercepts[,1]
##-     res <- c(res,
##-              list("(Intercepts)"=ciSgRl(lcfi, object$intercepts[,2], df=df, stcoef=lcfi) ))
##-   }
  ##  class(res) <- "termeffects" ## don't do that:
  ##                                 want to be able to print the whole table
  res
}
## ---------------------------------------
print.effecttable <- function (x, columns=NULL, ...)
{
  lcoldef <-
    if (rlvoptions("printstyle")=="relevance")
      rlvoptions("effectcolumns.r") else rlvoptions("effectcolumns.o")
  columns <- i.def(i.def(columns, rlvoptions("effectcolumns")), lcoldef) 
  print.termtable(x, columns=columns, ...)
}
## ---------------------------------------
print.termeffects <- function (x, columns=NULL, printstyle=NULL, transpose=FALSE,
                               single=FALSE, ...)
{
  if (is.null(columns))
    columns <-
      if (i.getopt(printstyle)=="relevance")
        rlvoptions("termeffcolumns.r") else rlvoptions("termeffcolumns.o")
  lnam <- names(x)
  for (li in seq_along(x)) {
    xi <- x[[li]]
    if (is.null(dim(xi))) next
    if (nrow(xi)==1) {
      if (!single) next
      cat("\n")
      }
    else cat("\n",lnam[li],":\n")
    print.termtable(xi, columns=columns, legend=FALSE, ...)
  }
}
## ----------------------------------------------------------
getscalepar <-
  function(object)
{ ## get scale parameter of a fit
  lsry <- summary(object)
  sigma <- c(lsry$sigma, lsry$scale)[1]
  if (length(sigma)==0) sigma <- sqrt(c(lsry$dispersion,1)[1])
  sigma
}
## -----------------------------------------------------------
getcoeffactor <-
  function(object)
{
  ## get factor for converting coef to coef effect
  ## model matrix
  lcls <- class(object)
  lmmt <- object[["x"]]
  if (length(lmmt)==0)  object$x <- lmmt <- model.matrix(object)
  lfamily <- object$family$family 
  ldist   <- object$dist
  lsigma <- getscalepar(object)
  lsigma <-
    if (any(lcls=="glm")&&lfamily%in%c("binomial", "quasibinomial"))
      1.6683*lsigma   ## qlogis(pnorm(1))
    else {
      lif <- any(lcls%in%c("lm","lmrob","rlm"))||
        (any(lcls=="survreg")&&ldist=="gaussian")
      if (length(lif)!=1) stop(traceback())
      if (lif) lsigma
      else 1
    }
  lfac <- apply(lmmt, 2, sd)/lsigma
  lfac[lfac==0] <- NA
##  lnm <- names(object$coefficients)
##-   if (any(lna <- is.na(match(lnm,names(lfac))))) {
##-     warning(":getcoeffactor: error, possibly singular case", paste(lnm[lna], collapse=", "))
##-     lfac <- 1
##-   } else lfac <- lfac[lnm] ## needed for singular designs
  structure(lfac, sigma=lsigma, fitclass=lcls, family=lfamily, dist=ldist)
}
## =================================================================================
rlvoptions <- #f
  function (x=NULL, list=NULL, default=NULL, options = NULL, assign=TRUE, ...)
{ ##
  loptdef <- rlv.optionsDefault 
  lnewo <- loldo <- i.def(options, get("rlv.options", envir=rlv.envir))
  ##
  if (length(list)&&length(lold <- attr(list, "old"))) list <- lold
  largs <- c(list, list(...))
  lop <- c(names(largs))
  ##
  if (length(x)) {  ## asking for options
    if(!is.character(x)) {
      warning(":rlvoptions: First argument 'x' must be of mode character")
      return(NULL)
    }
    return(if(length(x)==1) loldo[[x]] else loldo[x])
  }
  ## --- get default values
  if (length(default) && u.notfalse(default)) {
    if (u.true(default)) default <- "all"
    if (!is.character(default))
      stop("!rlvoptions! Argument 'default' must be of mode character or logical")
    if (default[1]=="all") largs = c(loptdef, largs)
    ## resets all available components
    else {
      if (default[1]=="unset") {
        largs <- c(loptdef[names(loptdef)%nin%names(loldo)], largs)
        default <- default[-1]
      }
        if (any(default!=""))
          largs <- c(largs, loptdef[default[default%in%names(loptdef)]])
    }
  }
  ## --- set options
  ## check
  if (length(largs)) largs <- check.option(list=largs)
  if (length(largs)) lnewo[names(largs)] <- largs
  lo <- intersect(names(largs),names(loldo))
  if (length(lo)) attr(lnewo, "old") <- loldo[lo]
  ## set margin pars, whether changed or not
  if (assign) assign("rlv.options", lnewo, pos=rlv.envir)
  ## assignInMyNamespace does not work
  invisible(lnewo)
}
## end or rlvoptions
## ---------------------------------------------------------------------------------
##- getopt <- function(x, opt = NULL, optdefault = rlv.optionsDefault) {
##-   ## x is character, opt list or NULL
##-   if (is.null(opt))  opt <- options()
##-   llx <- length(x)
##-   if ((!is.atomic(x))||!is.character(x))
##-     stop("!getopt! Argument 'x' not suitable")
##-   if (llx>1)
##-     lx <- lapply(x, getopt, opt=opt)
##-   else {
##-     lx <- opt[[x]]
##-     lx <- if (length(lx)==0)  getOption(x) else lx <- check.option(x, lx)
##-     if (length(lx)==0)  lx <- optdefault[[x]]
##-   }
##-   lx
##- }
## ====================================================================
getsymbol <- function(value, symbols)
  symbols$symbol[as.numeric(cut(value, symbols$cutpoint))]

i.getopt <- function(x, options = NULL) {
  if (is.null(options))
    options <- get("rlv.options", envir=rlv.envir) ## list in calling fn
  lnam <- as.character(substitute(x))
  lopt <- x
##  if (is.function(options)) options <- NULL
  if (is.null(lopt)||(is.atomic(lopt)&&all(is.na(lopt))))
    lopt <- options[[lnam]]
  if (is.null(lopt)||(is.atomic(lopt)&&all(is.na(lopt))))
    lopt <- rlvoptions(lnam)
  else unlist(check.option(lnam, lopt))   ## check
  if (is.null(lopt)) lopt <- rlv.optionsDefault[[lnam]]
  lopt
}

## ===========================================================================
regrModelClasses <- c("regr","lm","lmrob","rlm","glm","survreg","coxph","rq","polr")
pSymbols <- list(symbol=c("***", "**", "*", ".", " ", ""),
                  cutpoint=c(0, 0.001, 0.01, 0.05, 0.1, 1) )
pSymbols$legend <- paste(rbind(pSymbols$cutpoint,pSymbols$symbol), collapse="  ")
rlvSymbols <- list(symbol=c("-", " ", ".", "+", "++", "+++", ""),
                  cutpoint=c(-Inf,-1,0,1,2,5,Inf) )
rlvSymbols$legend <-
  paste(i.last(rbind(rlvSymbols$cutpoint,rlvSymbols$symbol)[-1],-2), collapse="  ")
## -----------------------------------------------------------
rlv.optionsDefault <- list(
  digits = 4,
  testlevel = 0.05,
  rlvThres = c(stand=0.1, rel=0.1, coef=0.1, drop=0.1, pred=0.05),
  termtable = TRUE, vif = TRUE,
  show.termeffects = TRUE, show.coefcorr = FALSE,
  termcolumns = NA, termeffcolumns = NA, effectcolumns = NA,
  printstyle = "relevance",
  termcolumns.r = c("coef",  "df", "coefRlp", "coefRls", ## "dropRle", 
                    "dropRls", "dropRls.symbol", "predRle"),
  ##, "predRlp", "predRls", "predRls.symbol"
  termeffcolumns.r = c("coef","coefRls.symbol"),
  effectcolumns.r = c("effect", "Rle", "Rlp", "Rls", "Rls.symbol"),
  termcolumns.o = c("coef",  "df", "ciLow","ciUp","R2.x", "Sig0", "p.value",
                  "p.symbol"),
  termeffcolumns.o = c("coef","p.symbol"),
  effectcolumns.o = c("effect", "ciLow","ciUp", "Sig0", "p.value", "p.symbol"),
  show.symbolLegend = TRUE,
  na.print = ".",
  pSymbols = pSymbols,
  rlvSymbols = rlvSymbols
)
## -----------------------
rlvoptionsCheck <- list(
  digits = cnr(c(3,10)),
  testlevel = cnr(c(0.0001,0.5)),
  rlvThres = cnr(c(0.0001,1)),
  termtable = clg(), vif = clg(),
  show.termeffects = clg(), show.coefcorr = clg(),
  termcolumns = clg(), termeffcolumns = clg(), effectcolumns = clg(),
  printstyle = cch(),
  termcolumns.r = cch(),
  termeffcolumns.r = cch(),
  effectcolumns.r = cch(),
  termcolumns.o = cch(),
  termeffcolumns.o = cch(),
  effectcolumns.o = cch(),
  show.symbolLegend = clg(),
  na.print = cch(),
  pSymbols = cls(c("symbol","cutpoint","legend")),
  rlvSymbols = cls(c("symbol","cutpoint","legend"))
)
## -----------------------------------------------------------------------
rlv.envir <- new.env()
rlv.envir$rlv.options <- rlv.envir$rlv.optionsDefault <- rlv.optionsDefault
rlv.pSymbols <- pSymbols
rlv.rlvSymbols <- rlvSymbols
## ========================================================================
