## relevance.R
## two groups and one sample
twosamples <- #f
  function(x, ...) UseMethod("twosamples")

twosamples.default <- #f
  function(x, y=NULL, paired=FALSE, hypothesis=0, var.equal=TRUE,
           testlevel=getOption("testlevel"), rlvThres=getOption("rlvThres"), ...)
{ ## effect (group difference) and relevance
  ltlev2 <- testlevel/2
  ## --------------
  lonegroup <- FALSE
  if (length(y)==0) {
    if (NCOL(x)==2) {
      y <- x[,2]
      x <- x[,1]
    } else {
      if (!is.atomic(x))
        stop("!twosamples! argument 'x' not suitable")
    }
  }
  if (paired) {
    if (length(x)!=length(y))
      stop("!twosamples! 'x' and 'y' must have equal length if 'paired' is TRUE")
    x <- y-x
    lonegroup <- TRUE
  }
  else if (length(y)==0) lonegroup <- TRUE
  ## ---
  ltb <- table(c(x,y))
  lbin <- length(ltb)<=2 ## all(c(x,y)%in%c(0,1,NA))
  x <- x[is.finite(x)]
  lmnx <- mean(x)
  ln <- lnx <- length(x)
  lny <- NULL
  if (lbin) { ## binary data
    lrlvth <-
      i.def(rlvThres["prop"], c(rlvThres, rlv.optionsDefault[["rlvThres"]]["prop"])[1])
    lx <- sum(x)
    if (lonegroup) { ## one proportion
      if (paired) x <- (x[x!=0]+1)/2 ## McNemar
      lt <- binom.test(lx,lnx, conf.level=1-testlevel)
      leff <- unname(qlogis(lt$estimate))
      leffci <- qlogis(lt$conf.int) ## c("ci.low","ci.up") ))
      lar <-
        qlogis( (qintpol(c(ltlev2, 1-ltlev2), lnx, plogis(hypothesis))+0:1) / lnx )
##                      c("accreg.low","accreg.up") ))
      lsg <- (leff-hypothesis)/abs(lar-hypothesis)
      lsig <- ifelse(leff<0, lsg[1], lsg[2])
      effname <- paste("log odds", if(paired) "of changes")
      method <- "One proportion, binomial inference"
    } else { ## tow proportions
      y <- y[is.finite(y)]
      lmny <- mean(y)
      lny <- length(y)
      ly <- sum(y)
      lt <- fisher.test(table(x,y), conf.int=TRUE, conf.level=1-testlevel)
      leff <- unname(log(lt$estimate))
      leffci <- log(lt$conf.int) ## c("ci.low","ci.up") ))
      lsig <- NA
      effname <- "log odds ratio"
      method <- "Two proportions, Fisher's exact inference"
    }
    lrlvci <- c(leff, leffci)/lrlvth
    ltst <- lt$statistic
    lpv <- lt$p.value
    lsigth <- NA
    ldf <- lv <- NULL
  } else { ## quantitative data
    lrlvth <-
      i.def(rlvThres["stand"], c(rlvThres, rlv.optionsDefault[["rlvThres"]]["stand"])[1])
    lssx <- sum((x-lmnx)^2)
    lpq <- 1-ltlev2
    if (lonegroup) { ## one quantitative sample
      leff <- lmnx
      ldf <- lnx-1
      lv <- lssx/ldf
      lny <- NULL
      effname <- paste("mean", if(paired) "of differences")
      method <- "One Sample t inference"
    } else { # two quantitative samples
      y <- y[is.finite(y)]
      lny <- length(y)
      ln <- lnx+lny
      lmny <- mean(y)
      leff <- lmny-lmnx
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
      effname <- "difference of means"
    }
    lse <- sqrt(lv/ln)
    ltst <- leff/lse
    lq <- qt(lpq, ldf)
    lciwid <- lq*lse
    leffci <- leff + c(-1,1)*lciwid
    lpv <- 2*pt(-abs(ltst), ldf)
    lrlvci <- c(leff,leffci)/sqrt(lv)/lrlvth
    lsig <- leff/lciwid
    lsigth <- c((lrlvci[1]-1)*2/diff(lrlvci[2:3]))
  }
  if (leff<0) lrlvci <- -lrlvci[c(1,3,2)]
  structure(
      c(effect=leff, ciLow=leffci[1], ciUp=leffci[2],
        Rle=lrlvci[1], Rls=lrlvci[2], Rlp=lrlvci[3],
        Sig0=lsig, Sigth=lsigth, p.value=lpv),
    class = "inference",
    method = method, effectname=effname, hypothesis = hypothesis,
    n = c(lnx, lny), means = if(!lonegroup) c(lmnx,lmny),
    statistic = ltst, V = lv, df = ldf,
    rlvThres =
      structure(lrlvth, type=ifelse(lbin, "proportion","standardized"))
  )
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
  ## DNAME <- paste(paste("`",names(mf),"'",sep=""), collapse = " by ")
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
##  attr(rr, "data.name") <- DNAME
  structure(rr, formula=formula, data.name=as.character(substitute(data)))
}
## ========================================================================
qintpol <- function(p, par1, par2, dist="binom")
{
  ## interpolated (theoretical) quantile for discrete distributions
  lqfn <- get(paste("q",dist,sep=""))
  lpfn <- get(paste("p",dist,sep=""))
  lq <- lqfn(p, par1, par2)
  lp <- lpfn(lq, par1, par2)
  lmaxx <- if (dist=="binom") par1 else Inf
  lp0 <- lpfn(pmin(pmax(0,lq-1),lmaxx), par1, par2)
  lq - (lp-p)/(lp-lp0)
}
## ===========================================================================
## regression
ciSgRl <-
  function (estimate=NULL, se=NULL, df=Inf, testlevel=0.05, stcoef=TRUE, rlv=TRUE,
            rlvThres=0.1, object=NULL)
{ ## for a coefficients table,
  ## calculate confidence interval, significance and relevance
  if (length(estimate)==0) estimate <- object
  if (inherits(estimate, regrModelClasses)) {
    object <- estimate
    estimate <- getcoeftable(object)
    df <- df.residual(object)
  }
  if (!(is.atomic(estimate)||length(dim(estimate))==2))
    stop("!ciSgRl! first argument not suitable")
  if (is.null(se))
    if (NCOL(estimate)>1) {
      se <- estimate[,2]
      estimate <- estimate[,1]
    }
  if (is.null(se)) se <- attr(estimate, "se")
  if (is.null(se))
    stop("!ciSgRl! no standard errors found")
  if (df==0||!any(is.finite(se))) {
    warning("!ciSgRl! no finite standard errors")
    return(cbind(estimate))
  }
  ltq <- qt(1-testlevel/2, df)
  lci <- estimate+outer(ltq*se, c(ciLow=-1,ciUp=1))
  ltst <- estimate/se
  lsgf <- ltst/ltq
  lpv <- 2*pt(-abs(ltst), df)
  rr <- data.frame(estimate=estimate, se=se, lci, testst=ltst,
                   Sig0=lsgf, p.value=lpv)
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
    stcoef <- estimate*lfac
  }
  if (length(stcoef)!=length(estimate)) {
    warning(":ciSgRl: argument 'stcoef' not suitable. No relevances")
    return (rr)
  }
  lfac <- stcoef/estimate
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
  cbind(rr, lrlv)
}
## -------------------------------------------------------------------------
print.inference <- #f
  function (x, show = getOption("show.ifc"),
            digits = getOption("digits.reduced"), ...)
{
  lshow <- i.getshow(show, type="inference")
  if (any(c("statistic","p.value","Sig0")%in%lshow)) lshow <- c(lshow,"test")
  cat("\n")
  cat(attr(x, "method"),"\n")
  out <- c(
    if (length(ldn <- attr(x, "data.name"))) paste("data: ", ldn),
    if (length(lfo <- attr(x, "formula"))) paste("formula: ", format(lfo))
  )
  if (length(out)) cat(paste(out, collapse=";  "), "\n")
  ## ---
  if (!is.na(leff <- x["effect"]))
    cat(if(length(leffn <- attr(x, "effectname")))
          paste(leffn, ": ") else "effect:   ", format(leff))
  if (getOption("show.confint"))
    if (length(lci <- x[c("ciLow","ciUp")]))
      cat(";  confidence int.: [",
          paste(format(lci), collapse=", "),"]\n")
  ## ---
  if ("test"%in%lshow) {
    cat("Test:     hypothesis: effect = ", attr(x,"hypothesis"), "\n  ")
    lpv <- x["p.value"]
    lps <-
      if (length(lpv)& ("p.symbol"%in%lshow | getOption("show.signif.stars")) )
             symnum(lpv, pSymbols$cutpoint, pSymbols$symbol)
    out <- c(if(length(ltst <- attr(x, "statistic")))
               paste("statistic: ", round(ltst, digits)),
             ## !!! df
             if("p.value"%in%lshow && length(lpv <- x["p.value"]))
               paste("p value: ", round(lpv, digits+1), lps),
             if("Sig0"%in%lshow&&length(lsig <- x["Sig0"]))
               paste("Sig0: ", round(lsig, digits),
                     if (!"p.value"%in%lshow) lps)
             )
    if (length(out)) cat(paste(out, collapse=" ;  "), "\n")
  }
  ## ---
  if ("relevance"%in%lshow) {
    cat("Relevances (threshold: ", attr(x, "rlvThres"), "):\n  ")
    lrs <- if ("Rls.symbol"%in%lshow)
             symnum(x["Rls"], rlvSymbols$cutpoint, rlvSymbols$symbol)
    out <- c(paste("Rle: ", round(x["Rle"], digits)),
             paste("Rlp: ", round(x["Rlp"], digits)),
             paste("Rls: ", round(x["Rls"], digits), lrs)
             )
    if (length(out)) cat(paste(out, collapse=" ;  "), "\n")
  }
}
## -----------------------------------------------------------
i.getshow <- #f
  function(show, type=c("inference", "terms", "termeffects"))
{ ## collect items to be shown by print.inference
  ##-   lstyles <- c("test", "relevance", "classical")
##-   li <- pmatch(show, lstyles, nomatch=0)
##-   lcoll <- lstyles[li] ## careful: explicit columns might pmatch!
  lcoll <- intersect(show, c("test", "relevance", "classical"))
  if(length(lcoll)) {
    ltype <- pmatch(type, c("inference", "terms", "termeffects"), nomatch=0)
    ## if (length(ltype)==0)
    lc <- paste("show", c("ifc","term","termeff")[ltype], lcoll, sep=".")
    for (l in lc)
      show <- c(show, getOption(l))
  }
  setdiff(show,lcoll)
}
## =============================================================================
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
termtable <- #f
  function (object, summary=NULL, testtype=NULL, r2x = TRUE, rlv = TRUE,
            rlvThres=c(rel=0.1, coef=0.1, drop=0.1, pred=0.05), testlevel=0.05)
{
  ## Purpose:  generate term table for various models
  ## --------------------------------------------------------
  ## prep
  if (length(summary)==0) summary <- summary(object)
  ## thresholds
  ltlev1 <- 1-testlevel
  ## lrlthrl <- rlvThres[["rel"]]
  ## lrlthcf <- rlvThres[["coef"]]
  browser()
  lrlthdr <- rlvThres[["drop"]]
  lrlthpr <- rlvThres[["pred"]]
  ## --- sigma and threshold for standardized coefficient
  ## lsigma <- c(object$sigma, summary$sigma, 1)[1]
  lfamily <- if (is.character(lfm <- object$family)) lfm else lfm$family
  ldist <- object$dist
  ## if ((inherits(object, "glm") &&
  ##      lfamily %in% c("poisson","quasipoisson")) ||
  ##     (inherits(object, "survreg")&&
  ##      ldist %in% c("weibull", "exponential", "lognormal"))
  ##     )  lrlthcf <- lrlthrl
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
      testst=NA, p.value=NA, R2.x=NA,
      stringsAsFactors=FALSE)
      )
  ## degrees of freedom
  ldfres <- df.residual(object)
  if (ldfres<1) {
    warning(":termtable: no degrees of freedom left.")
    return(data.frame(
      coef=c(object$coefficients,NA)[1], df = NA, se=NA, ciLow=NA, ciUp=NA,
      Sig0=NA, stcoef=NA, stciLow = NA, stciUp = NA,
      testst=NA, p.value=NA, R2.x=NA,
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
      ## lcov <- summary$cov.unscaled
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
                    testst=ltst, p.value=lpv, R2.x=lr2, ldrrl,
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
    ltstq <- c(ltstq1, ltstq)
    lcont <- c(0, lcont)
  } else if (!inherits(object, c("coxph", "polr")))
    warning(":termtable: No intercept. Statistics are difficult to interpret.")
  lcont1 <- lcont+ljint  # row number in dr1
  ## --- coefficients and statistics for terms with 1 df
  if (length(lcont)) { ## lcont refers to assign
    ## ltlb <- dimnames(ltb)[[1]]
    ## lclb <- ltlb[lcont1] ## lcont1 is the row in the coef table of summary(object)
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
  structure(ltb, class=c("termtable", "data.frame"),
            testtype=testtype,
            fitclass=class(object), family=lfamily, dist=ldist,
            rlvThres=rlvThres)
}
## end termtable
## ===========================================================
confintF <- #f
  function(f, df1, df2=Inf, testlevel=0.05) {
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
print.termtable <- #f
  function(x, show = getOption("show.ifc"), transpose.ok = FALSE, legend = NULL,
           digits = getOption("digits.reduced"), na.print = "  ", ...)
{
  ## -------
  ff.mergesy <- function(lx, x, var, place)
  {
    lxx <- x[,var]
    ltypep <- substring(var,1,3)=="p.v"
    lxx[is.na(lxx)] <- ltypep
    lsymb <- if(ltypep) pSymbols else rlvSymbols
    ls <- symnum(lxx, lsymb$cutpoint, lsymb$symbol)
    lsy <- format(c("    ",ls))[-1]
    ## trick needed to ensure flush left priniting of the symbols
    lcl <- match(c(var, place, "coef"), colnames(lx), nomatch=0)[1]
    rr <- if(lcl<ncol(lx)) cbind(lx[,1:lcl], lsy, lx[,-(1:lcl),drop=FALSE])
          else  cbind(lx[,1:lcl], rsy=lsy)
    names(rr)[lcl+1] <- paste(substring(var,1,1),if (!ltypep) "R", "sy",sep="")
    attr(rr, if(ltypep) "pLegend" else "rlvLegend") <- attr(ls, "legend")
    rr
  }
    ## -------------------------------------
  digits <- pmax(digits, 3)
  show <- if (length(show)==1 && show=="all")
            c(colnames(x),
              "p.symbol", "coefRls.symbol", "dropRls.symbol", "predRls.symbol")
  else unique(c("coef", i.getshow(show, "terms")))
  show[show=="estimate"] <- "coef"
  colnames(x)[colnames(x)=="estimate"] <- "coef"
  lcols <- intersect(show, colnames(x))
  if (length(lcols)==0) {
    warning(":print.termtable: no columns selected")
    return(invisible(x))
  }
  ## ---
  lx <- as.data.frame(x)[,lcols, drop=FALSE]
  ## --- round some columns to 3 digits
  ljrp <- lcols[pmatch(c("R2","p.v"), lcols, nomatch=0)]
  if (length(ljrp)) lx[,ljrp] <- round(as.matrix(lx[,ljrp]),digits)
  ljrp <- lcols[c(grep("Rl", lcols),grep("Sig", lcols))]
  if (length(ljrp)) lx[,ljrp] <- round(as.matrix(lx[,ljrp]),i.last(digits)-1)
  ## --- paste symbols to numbers
  lpleg <- lrleg <- NULL
  if ("p.symbol"%in%show & "p.value"%in%colnames(x)) {
    lx <- ff.mergesy(lx, x, "p.value", "Sig0")
##-     lpv <- x[,"p.value"]
##-     lpv[is.na(lpv)] <- 1
##-     lpsy <- symnum(lpv, pSymbols$cutpoint, pSymbols$symbol)
##-     lcl <- match(c("p.value","Sig0","coef"), lcols, nomatch=0)[1]
##-     lx <- if(lcl<ncol(lx)) cbind(lx[,1:lcl], psy=lpsy, lx[,-(1:lcl),drop=FALSE])
##-           else  cbind(lx[,1:lcl], psy=lpsy)
  }
  li <- grep("Rls.symbol", show)
  if (length(li)) {
    lrs <- show[li]
    lrc <- sub(".symbol","",lrs)
    for (lii in seq_along(li)) {
      lx <- ff.mergesy(lx, x, lrc[lii], c("dropRls", "coefRls"))
    }
  }
  ## ---------------
  ## print
  if (transpose.ok &&
      ((lnc1 <- ncol(lx)==1)|| ncol(lx)==2 && length(grep(".symbol", show))) ) {
    if (lnc1) print(lx[[1]])
    else print(setNames(paste(format(lx[[1]]),lx[[2]]),
                        paste(row.names(lx),"    ")), quote=FALSE)
  } else
    print(lx, quote=FALSE, na.print=na.print)
  ## --- legend(s)
  if(u.notfalse(legend)) {
    if((u.true(legend)|(is.null(legend) & getOption("show.signif.stars"))) &&
       length(ll <- attr(lx, "pLegend")))
      cat("\nSignificance codes for p.value:  ", ll, "\n")
    if((u.true(legend)|(is.null(legend) & getOption("show.symbollegend"))) &&
       length(ll <- attr(lx, "rlvLegend")))
      cat("\nRelevance codes:    ", ll, "\n")
  }
  invisible(lx)
}
## ===========================================================================
termeffects <- #f
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
  if (inherits(object, "polr")) {
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
    } else if (inherits(object, "polr"))
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
##    mmc <- names(asgn)[asgn == j & !is.na(coef)]  ## lcols (logical fails for polr, vcov() too large) !!! was  which
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
          li <- match(c("Rle","Rlp","Rls"), names(rrj))
          names(rrj)[li] <- c("coefRle","coefRlp","coefRls")
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
##-   if (inherits(object, "polr")) {
##-     lcfi <- object$intercepts[,1]
##-     res <- c(res,
##-              list("(Intercepts)"=ciSgRl(lcfi, object$intercepts[,2], df=df, stcoef=lcfi) ))
##-   }
  ##  class(res) <- "termeffects" ## don't do that:
  ##                                 want to be able to print the whole table
  res
}
## ---------------------------------------
print.termeffects <- #f
  function (x, show = getOption("show.ifc"), transpose.ok=TRUE, single=FALSE, ...)
{
  show <- if (length(show)==1 && show=="all")
            c(colnames(x), "p.symbol", "coefRls.symbol")
  else unique(c("coef", i.getshow(show, "termeffects")))
  lnam <- names(x)
  for (li in seq_along(x)) {
    xi <- x[[li]]
    if (is.null(dim(xi))) next
    if (nrow(xi)==1) {
      if (!single) next
      cat("\n")
      }
    else cat("\n",lnam[li],":\n")
    print.termtable(xi, show=show, legend=FALSE, transpose.ok=transpose.ok, ...)
  }
}
## ----------------------------------------------------------
getscalepar <- #f
  function(object)
{ ## get scale parameter of a fit
  lsry <- summary(object)
  sigma <- c(lsry$sigma, lsry$scale)[1]
  if (length(sigma)==0) sigma <- sqrt(c(lsry$dispersion,1)[1])
  sigma
}
## -----------------------------------------------------------
getcoeffactor <- #f
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
## ====================================================================
pSymbols <- list(symbol=c("***", "**", "*", ".", " "),
                  cutpoint=c(0, 0.001, 0.01, 0.05, 0.1, 1) )
rlvSymbols <- list(symbol=c(" ", ".", "+", "++", "+++"),
                  cutpoint=c(-Inf,0,1,2,5,Inf) )
## -----------------------------------------------------------
rlv.optionsDefault <- list(
  digits.reduced = 3,
  testlevel = 0.05,
  rlvThres = c(stand=0.1, rel=0.1, prop=0.1, coef=0.1, drop=0.1, pred=0.05),
  show.confint = TRUE,
  termtable = TRUE, vif = TRUE,
##-   show.termeffects = TRUE, show.coefcorr = FALSE,
  show.ifc = "relevance",
  show.ifc.relevance = c("Rle", "Rlp", "Rls", "Rls.symbol"),
  show.ifc.test = c("Sig0", "p.symbol"),
  show.ifc.classical = c("statistic", "p.value", "p.symbol"),
  show.term.relevance = c("df", "R2.x", "coefRlp", "coefRls", ## "dropRle",
                         "dropRls", "dropRls.symbol", "predRle"),
  show.term.test = c("df", "ciLow","ciUp", "R2.x", "Sig0", "p.value",
                         "p.symbol"),
  show.term.classical = c("df", "se", "statistic", "p.value", "p.symbol"),
  show.termeff.relevance = c("coef","coefRls.symbol"),
  show.termeff.test = c("coef","p.symbol"),
  show.termeff.classical = c("coef","p.symbol"),
  show.symbollegend = TRUE,
  show.doc = TRUE,
  na.print = ".",
  pSymbols = pSymbols,
  rlvSymbols = rlvSymbols
)
.onLoad <- function(lib, pkg) options(rlv.optionsDefault)
## -----------------------
