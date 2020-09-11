## relevance.R

rlstats <-
  function (x=NULL, se=NULL, df=Inf, testlevel=NULL, steff=NULL, rlv=TRUE,
            rlvthres = NULL, object=NULL)
{ ## for a coefficients table,
  ## calculate confidence interval, significance and relevance
  testlevel <- i.getopt(testlevel)
  if (length(x)==0) x <- object
  if (inherits(x, regrModelClasses)) {
    object <- x
    x <- getcoeftable(object)
  }
  if (!(is.atomic(x)||length(dim(x))==2))
    stop("!rlstats! first argument not suitable")
  lpvs <- getopt(c("p.cutpoints","p.symbols"))
  lrls <- getopt(c("rlv.cutpoints","rlv.symbols"))
  if (u.isnull(se))
    if (NCOL(x)>1) {
      se <- x[,2]
      x <- x[,1]
    }
  if (u.isnull(se)) se <- attr(x, "se")
  if (u.isnull(se))
    stop("!rlstats! no standard errors found")
  if (df==0||!any(is.finite(se))) {
    warning("!rlstats! no finite standard errors")
    return(cbind(x))
  }
  ltq <- qt(1-testlevel/2, df)
  lci <- x+outer(ltq*se, c(ciLow=-1,ciUp=1))
  ltst <- x/se
  lsgf <- ltst/ltq
  lpv <- 2*pt(-abs(ltst), df)
  lpsymb <- lpvs[[2]][as.numeric(cut(lpv, lpvs[[1]]))]
  rr <- data.frame(effect=x, se=se, lci, testst=ltst,
                   signif0=lsgf, p.value=lpv, p.symbol=lpsymb)
  attr(rr, "p.legend") <- getopt("p.legend")
  ## --- standardized coefficients
  if (!u.notfalse(steff)) return (rr)
  ##    warning(":ciSigRlv: no standardized coefficients given. No relevances")
  ljst <- TRUE
  if (u.true(steff)&&length(object)==0) {
    warning(":rlstats: argument 'object' needed. No standardized effects")
    ljst <- FALSE
  } else {
    if ((length(steff)==0 & length(object)>0)||u.true(steff)) {
      if (inherits(object, "nls")) {
        warning(":rlstats: cannot calculate standardized coefficients for 'nls' objects.")
        ljst <- FALSE
      } else {
        lfac <- getcoeffactor(object)
        ##   lcls <- attr(lfac, "fitclass")
        steff <- x*lfac
      }
    }
  }
  if (ljst&&length(steff)!=length(x)) {
    warning(":rlstats: argument 'steff' not suitable. No standardized effects")
    ljst <- FALSE
  }
  if (ljst) {
    lfac <- steff/x
    lstci <- cbind(steffect=steff, st.Low=lfac*lci[,1], st.Up=lfac*lci[,2])
    rr <- cbind(rr, lstci)
  } else lstci <- rr[,c("effect","ciLow","ciUp")]
  ## --- relevance
  if (!rlv) return(rr)
  rlvthres <- i.getopt(rlvthres)
  if (length(rlvthres)>1) rlvthres <- i.def(rlvthres["stand"], getopt("rlvthres")$stand)
  if (is.na(rlvthres)) {
    warning(":rlstats: argument 'rlvthres' not suitable. No relevances")
    return (rr)
  }
  lrlv <- lstci/rlvthres
  li <- which(lstci[,1]<0)
  if (length(li)) lrlv[li,] <- - lrlv[li,c(1,3,2)]
  colnames(lrlv) <- c("Rle","Rls","Rlp")
  lrlvsy <-
    if (any(!is.na(ll <- lrlv[,"Rls"])))
      lrls[[2]][as.numeric(cut(ll, lrls[[1]]))] else NA
  rr <- cbind(rr, lrlv, Rls.symbol=lrlvsy)
  attr(rr, "p.legend") <- getopt("p.legend")
  attr(rr, "rlv.legend") <- getopt("rlv.legend")
  rr
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
  ltb
}
## --------------------------------------------------------------------------
termtable <-
  function (object, summary=NULL, testtype=NULL, r2x = TRUE, rlv = TRUE,
            rlvthres = NULL, testlevel = NULL)
{
  ## Purpose:  generate term table for various models
  ## --------------------------------------------------------
  ## prep
  if (length(summary)==0) summary <- summary(object)
  ## thresholds
  testlevel <- i.getopt(testlevel)
  ltlev1 <- 1-testlevel
  rlvthres <- i.getopt(rlvthres)
  lrlthrl <- rlvthres["rel"]
  lrlthcf <- rlvthres["stand"]
  lrlthdr <- rlvthres["drop"]
  lrlthpr <- rlvthres["pred"]
  ## --- sigma and threshold for standardized coefficient
  lsigma <- c(object$sigma, summary$sigma, 1)[1]
  lfamily <- object$family$family
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
      testtype <- "LRT"
      if (object$family %in% c("quasibinomial","quasipoisson")) testtype <- "F"
    }
    if (inherits(object, c("survreg"))) testtype <- "Chisq"
  }
  ## ---
  lterms <- terms(object)
  if(length(attr(lterms,"term.labels"))==0)
    return(data.frame(
      coef=c(object$coefficients,NA)[1], df = NA, se=NA, ciLow=NA, ciUp=NA,
      signif0=NA, stcoef=NA, stciLow = NA, stciUp = NA,
      testst=NA, p.value=NA, p.symbol="", R2.x=NA,
      stringsAsFactors=FALSE)
      )
  ## degrees of freedom
  ldfres <- df.residual(object)
  if (ldfres<1) {
    warning(":termtable: no degrees of freedom left.")
    return(data.frame(
      coef=c(object$coefficients,NA)[1], df = NA, se=NA, ciLow=NA, ciUp=NA,
      signif0=NA, stcoef=NA, stciLow = NA, stciUp = NA,
      testst=NA, p.value=NA, p.symbol="", R2.x=NA,
      stringsAsFactors=FALSE)
      )
  }
  ## --- coefficients
  lcoef <- object$coefficients
##-   lcoeftab <- object$coeftable
##-   if (length(lcoeftab)==0||ncol(lcoeftab)<10)
  lcoeftab <- rlstats(object)
  ljc <- match(c("effect","steffect"), names(lcoeftab), nomatch=0)
  names(lcoeftab)[ljc] <- c("coef", "stcoef")[ljc!=0]
  ## --- drop1
  ldr1 <-
    if (inherits(object, c("lm","lmrob"))&&!inherits(object, "glm")) {
      lcov <- summary$cov.unscaled
      try(drop1Wald(object, test=testtype, scope=lterms), silent=TRUE)
    } else {
      try(drop1(object, test=testtype, scope=lterms), silent=TRUE)
    }
  if (inherits(ldr1, "try-error")) {
    warning(":regr/termtable: drop1 did not work. I return the codfficient table")
##                  produced by ", object$fitfun
    return(list(termtable=lcoeftab))
  }
  ldr1 <- ldr1[-1,]
  ldr1$RSS <- NULL # same ncol for lm and glm
  if (inherits(object,"rlm"))  ldr1[,4] <- ldr1[,2]/ldr1[,1]
  ## -- critical value for test
  ltstq <- if (testtype=="F") qf(ltlev1,c(1,ldr1[,1]),ldfres) else {
    if (testtype=="Chisq") qchisq(ltlev1,c(1,ldr1[,1])) else NA }
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
  lnobs1 <- ldfres+sum(ldf)-1
  ## drop effect relevance
  ltst <- ldr1[, ifelse(inherits(object, "polr"), 3, 4)]
  ldrncci <- rbind(confintF(ltst, ldf, ldfres, testlevel))
  ldreff2 <- cbind(ltst*ldf, ldrncci)/lnobs1
  ldrrl <- cbind(NA,NA,NA, sqrt(ldreff2)/lrlthdr,
                 pmax(0.5*log((ldfres+ldreff2*lnobs1)/(ldfres+ldf))/lrlthpr, 0))
  dimnames(ldrrl) <-
    list(NULL, c(t(outer(c("","drop","pred"),c("Rle","Rls","Rlp"),paste, sep=""))))
  ltst <- ldr1[,lpvcol-1]
  ## table, filled partially
    ltb <- data.frame(coef=NA, df=ldr1[,1], se=NA, ciLow=NA, ciUp=NA,
                      signif0=sqrt(pmax(0,ltst)/ltstq),
                      stcoef=NA, st.Low=NA, st.Up=NA,
                      testst=ltst, p.value=lpv, p.symbol="", R2.x=lr2, ldrrl,
                      stringsAsFactors=FALSE)
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
            "Rle","Rls","Rlp")
    ltb[lcont1,lj] <- lcoeftab[ljc,lj]
    ltb[lcont1,"signif0"] <- sign(ltb[lcont1,"coef"])*ltb[lcont1,"signif0"]
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
  lpvs <- getopt(c("p.cutpoints","p.symbols"))
  ltb[,"p.symbol"] <-
    if (any(!is.na(ll <- ltb$p.value)))
      lpvs[[2]][as.numeric(cut(ll, lpvs[[1]]))] else NA
  ## Rls-symbol
  lrls <- getopt(c("rlv.cutpoints","rlv.symbols"))
##-   ltb[,"Rls.symbol"] <-
##-     if (any(!is.na(ll <- ltb$Rls)))
##-       lrls$symbol[as.numeric(cut(ll, lrls$cutpoint))] else NA
  ltb[,"dropRls.symbol"] <-
    if (any(!is.na(ll <- ltb$dropRls)))
      lrls[[2]][as.numeric(cut(ll, lrls[[1]]))] else NA
  ltb[,"predRls.symbol"] <-
    if (any(!is.na(ll <- ltb$predRls)))
      lrls[[2]][as.numeric(cut(ll, lrls[[1]]))] else NA
  ## ---
  structure(ltb, class=c("termtable", "data.frame"),
            testtype=testtype,
            fitclass=class(object), family=lfamily, dist=ldist,
            rlvthres=rlvthres,
            p.legend=getopt("p.legend"), rlv.legend=getopt("rlv.legend"))
}
## end termtable
## ===========================================================
confintF <-
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
  if (u.isnull(columns))
    columns <-
      if (i.getopt(printstyle)=="relevance")
        getopt("termcolumns.r") else getopt("termcolumns.p")
  columns[columns=="estimate"] <- "coef"
  names(x)[names(x)=="estimate"] <- "coef"
  if (length(columns)==1)
    columns <-
      switch (columns, all = names(x), relevance = getopt("termcolumns.r"),
              convent = getopt("termcolumns.p"), columns)
  if (length(ll <- setdiff(columns, names(x)))) {
    warning(":print.termtable: columns not found:  ", paste(ll, collapse=", "))
    if (length(lcol <- setdiff(columns, ll))) columns <- lcol else return()
  }
  lattr <- attributes(x)
  ## ---
  x <- data.frame(x)[,columns]
  lcnames <- colnames(x)
  ## --- round some columns to 3 digits
  ljrp <- lcnames[dropNA(pmatch(c("R2","signif0","p.v"), colnames(x)))]
  if (length(ljrp))
    x[,ljrp] <- round(as.matrix(x[,ljrp]),max(3,digits))
  if ("signif0" %in% ljrp) x$signif0 <- round(x$signif0,last(digits)-1)
  ## --- paste symbols to numbers
  ljsy <- grep(".symbol", lcnames)
  if (any(ll <- lcnames[ljsy]=="p.symbol"))
    lcnames[ljsy[ll]] <-
      if ("p.value" %in% lcnames) "p.value.symbol"
      else {
        if ("signif0" %in% lcnames) "signif0.symbol" else "coef.symbol"
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
  lprleg <- i.def(legend, getopt("show.symbolLegend"))
  if (lprleg) {
    if (length(grep(".symbol", columns))) {
      ## cat("---")
      if ("p.symbol" %in% columns)
        cat("\nSignificance codes:  ", lattr[["p.legend"]],"\n", sep = "")
      if (any(c("Rls.symbol", "dropRls.symbol", "predRls.symbol") %in% columns))
        cat("\nRelevance codes:    ", lattr[["rlv.legend"]],"\n", sep = "")
      cat("\n")
    }
  }
  invisible(xp)
}
## ===========================================================================
termeffects <-
  function (object, se = 2, df = df.residual(object), rlv = TRUE, ...)
    ## --------------------------------------------------------------
{
  if (is.atomic(object)||u.isnull(terms(object)))
      stop("!termeffects! inadequate first argument")
 ##  xl <- object$xlevels
  Terms <- delete.response(terms(object))
  tl <- attr(Terms, "term.labels")
  dcl <- attr(Terms,"dataClasses")[-1]
  if (all(dcl=="numeric")) ## ??? need coeftable. is it always available?
    return(as.list(coef(object))) ## !!! coeftable!
  ## result already available?
  allc <- object$termeffects
  if ((!u.isnull(allc))&&length(allc)==length(tl)&&
      (is.matrix(allc[[length(allc)]])|!se)) return(allc) ## !!! check!
  ## ---
  if (rlv) lsigma <- getscalepar(object)
  int <- attr(Terms, "intercept")
  facs <- attr(Terms, "factors")
  mf <- object$model  ##! d.c used all.vars
  if (u.isnull(mf)) mf <- model.frame(object)
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
    for (li in seq_along(imat))
      lctr <- c(lctr, list(diag(length(xtlv[[li]]))))
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
  lnna <- is.finite(coef)
  names(asgn) <- colnames(mm)
  if (any(!lnna)){
    coef <- coef[lnna]
    mm <- mm[,lnna]
    asgn <- asgn[lnna]
  }
  if (se) {
    cov <- vcov(object)
    if (u.isnull(cov)) {
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
    mmpart <- mm[mmr,mmc, drop=FALSE]
    rrj <- setNames(drop(mmpart %*% lcf), rnn[mmr]) ## coef[mmc]
    if (se) {
      sej <- sqrt(diag(mmpart %*% cov[mmc,mmc] %*% t(mmpart)))
      if (any(is.na(rrj))|any(!is.finite(sej))) {
##-         warning(":termeffects: missing coef or non-finite standard error for term '",
##-                 tl[j], "'. no standard errors etc")
        ljfail <- c(ljfail, tl[j])
      } else {
        rrj <- rlstats(rrj, sej, df, steff=rrj*0.5/lsigma, rlv=rlv)
      }
    }
    res[[j]] <- rrj
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
##-              list("(Intercepts)"=rlstats(lcfi, object$intercepts[,2], df=df, stcoef=lcfi) ))
##-   }
  ##  class(res) <- "termeffects" ## don't do that:
  ##                                 want to be able to print the whole table
  res
}
## ---------------------------------------
print.termeffects <-
  function (x, columns=NULL, printstyle=NULL, transpose=FALSE, single=FALSE, ...)
{
  if (u.isnull(columns))
    columns <-
      if (i.getopt(printstyle)=="relevance")
        getopt("termeffcolumns.r") else getopt("termeffcolumns.p")
  lnam <- names(x)
  for (li in seq_along(x)) {
    xi <- x[[li]]
    if (u.isnull(dim(xi))) next
    if (nrow(xi)==1) {
      if (!single) next
      cat("\n")
      }
    else cat("\n",lnam[li],":\n")
    print.termtable(xi, columns=columns, legend=FALSE)
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
  structure(lfac, sigma=lsigma, fitclass=lcls, family=lfamily, dist=ldist)
      ## model.matrix=lmmt, 
}
## ==========================================================================
drop1Wald <-
  function (object, scope=NULL, scale = 0, test = c("none", "Chisq", "F"),
           k = 2, ...)
{
    x <- model.matrix(object)
    offset <- model.offset(model.frame(object))
    n <- nrow(x)
    asgn <- attr(x, "assign")
    lterms <- terms(object)
    tl <- attr(lterms, "term.labels")
    attr(lterms, "order") <-  rep(1,length(tl))
    if (is.null(scope))
      scope <- tl # drop.scope(lterms)
    else {
        if (!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)),
                "term.labels")
        if (!all(match(scope, tl, 0L) > 0L))
            stop("scope is not a subset of term labels")
    }
    lsry <- summary(object)
    ndrop <- match(scope, tl)
    ns <- length(scope)
    rdf <- object$df.residual
    lsig <- c(lsry$sigma, lsry$scale, 1)[1]
    chisq <- lsig^2 * rdf
    ## sum(weighted.residuals(object)^2, na.rm = TRUE)
    ## deviance.lm(object)
    dfs <- numeric(ns)
    RSS <- numeric(ns)
    cov <- object$cov.unscaled
    if (is.null(cov)) cov <- summary(object)$cov.unscaled
    if (is.null(cov)) cov <- object$covariance/lsig^2
    if (length(cov)==0) stop("!drop1Wald! no covariance matrix found")
    cf <- object$coefficients
    ##- jj <- match(names(cf),colnames(cov), nomatch=0)
    ##-     if (!(any(jj==0)&&all(is.na(cf[jj==0]))))
    jj <- match(colnames(cov),names(cf), nomatch=0)
    if (any(jj==0))
      warning(":drop1Wald: coefficient(s) and cov. matrix may not correspond")
    coef <- cf[jj]
    asgn <- asgn[jj]
    if (any(names(coef[!is.na(coef)])%nin%names(coef)))
      stop("!drop1Wald! coefficient(s) not appearing in covariance matrix")
##-     y <- object$residuals + predict(object)
    for (i in 1:ns) {
        ii <- which(asgn==ndrop[i]) ## seq_along(asgn)[asgn == ndrop[i]]
        RSS[i] <- if (length(ii)==1) coef[ii]^2/cov[ii,ii] else
          coef[ii]%*%solve(cov[ii,ii])%*%coef[ii]  ## !!! REPLACE THIS
        dfs[i] <- length(ii)
##-         if (all.cols)
##-             jj <- setdiff(seq(ncol(x)), ii)
##-         else jj <- setdiff(na.coef, ii)
##-         z <- if (iswt)
##-             lm.wfit(x[, jj, drop = FALSE], y, wt, offset = offset)
##-         else lm.fit(x[, jj, drop = FALSE], y, offset = offset)
##-         dfs[i] <- z$rank
##-         oldClass(z) <- "lm"
##-         RSS[i] <- deviance(z)
    }
    scope <- c("<none>", scope)
    dfs <- c(c(object$rank,object$df)[1], dfs)
    RSS <- chisq + c(0, RSS)
    if (scale > 0)
        AIC <- RSS/scale - n + k * dfs
    else AIC <- n * log(RSS/n) + k * dfs
##-     dfs <- dfs[1] - dfs
##-     dfs[1] <- NA
    aod <- data.frame(Df = dfs, "Sum of Sq" = c(NA, RSS[-1] -
        RSS[1]), RSS = RSS, AIC = AIC, row.names = scope, check.names = FALSE)
    if (scale > 0)
        names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
    test <- match.arg(test)
    if (test == "Chisq") {
        dev <- aod$"Sum of Sq"
        if (scale == 0) {
            dev <- n * log(RSS/n)
            dev <- dev - dev[1]
            dev[1] <- NA
        }
        else dev <- dev/scale
        df <- aod$Df
        nas <- !is.na(df)
        dev[nas] <- pchisq(dev[nas], df[nas], lower.tail = FALSE)
        aod[, "Pr(Chi)"] <- dev
    }
    else if (test == "F") {
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
    head <- c("Single term deletions (Wald test)", "\nModel:",
              deparse(as.vector(formula(object))),
        if (scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
## ===================================================
vif.regr <- function (object, cov=NULL, mmat=NULL)
{
  ## Purpose:   vif.lm  of library  car
  ## ----------------------------------------------------------------------
  ## Author: objectified by Werner Stahel, Date: 11 Mar 2005, 09:18
  terms <- labels(terms(object))
  n.terms <- length(terms)
  if (n.terms < 2) {
    ##-         stop("model contains fewer than 2 terms")
    return(matrix(1,1,3))
  }
  if (length(cov)==0) {
    cov <- object$cov.unscaled
    if (is.null(cov)) cov <- summary(object)$cov.unscaled
    if (is.null(cov)) cov <- object$covariance ## /lsig^2 # no: a factor does not matter...
    if (is.null(cov)) cov <- object$var ## survreg
    if (length(cov)==0) {
      warning("!vif.regr! no covariance matrix found")
      return(NULL)
    }
  }
  if (length(mmat)==0) mmat <- model.matrix(object)
  if (length(mmat)==0) {
    warning("!vif.regr! no model matrix found")
    return(NULL)
  }
  cls <- dimnames(mmat)[[2]]%in%dimnames(cov)[[2]]
  ##-                                         # needed for singular cases
  assign <- attr(mmat,"assign")[cls]
  if (names(coefficients(object)[1]) == "(Intercept)") {
    cov <- cov[-1, -1]
    assign <- assign[-1]
  }
  else if (object$fitfun%nin%c("polr","coxph","survreg"))
    warning("No intercept: vifs may not be sensible.")
  sd <- 1/sqrt(diag(cov))
  if (any(!is.finite(sd))) {
    warning(":vif.regr: zero variances of predictors. no R2x")
    return(NULL)
  }
  R <- cov/outer(sd,sd)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) *
      det(as.matrix(R[-subs,-subs]))/detR
    result[term, 2] <- length(subs)
  }
  result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  result
}
## ==========================================================================
last <- function (data,n = NULL, ncol=NULL, drop=is.matrix(data))
{
  ldim <- dim(data)
  if (is.null(ldim)) {
    if (is.null(n)) n <- 1
    ldt <- length(data)
    return(data[sign(n)*((ldt-abs(n)+1):ldt)])
  }
  if (length(ldim)!=2)
    stop ("!last! not programmed for arrays of dimension >2")
  if (is.null(n)&is.null(ncol)) n <- 1
  if (is.null(n)) n <- ldim[1]
  if (is.null(ncol)) ncol <- ldim[2]
  data[sign(n)*((ldim[1]-abs(n)+1):ldim[1]),
       sign(ncol)*((ldim[2]-abs(ncol)+1):ldim[2]), drop=drop]
}
## -----------------------------------------------------------------
dropNA <- function (x, inf=TRUE) {
  if (length(dim(x))) {
    if (is.numeric(x)&inf) x[apply(is.finite(x),1,all)]
    else x[!apply(as.matrix(is.na(x)), 1, any),] ## as.matrix needed for Surv obj
  } else if (is.numeric(x)&inf) x[is.finite(x)] else x[!is.na(x)]
}
## ------------------------------------------------------------------
sumNA <- function (object, inf=TRUE)
{
  ## Purpose:   count NAs along columns
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   object   data.frame, matrix or vector
  ##   inf      treat Inf as NA
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 10 Oct 2007, 08:18
  ff <- if(inf) {
    function(x)
    if(is.numeric(x)) sum(!is.finite(x)) else sum(is.na(x)) }
      else function(x) sum(is.na(x))
  if (is.matrix(object)) apply(object,2,ff)  else {
    if (is.list(object)) sapply(object,ff)
    else if(is.atomic(object)) ff(object)
  }
}
## -----------------------------------------------------
getopt <- function(x, opt = NULL) {
  ## x is character, opt list or NULL
  if (is.null(opt))  opt <- options()
  llx <- length(x)
  if ((!is.atomic(x))||!is.character(x))
    stop("!getopt! Argument 'x' not suitable")
  if (llx>1)
    lx <- lapply(x, getopt, opt=opt)
  else {
    lx <- opt[[x]]
    lx <- if (length(lx)==0)  getOption(x) else lx <- checkopt(x, lx)
    if (length(lx)==0)  lx <- optdefault[[x]]
  }
  lx
}
## ---------------------------------------------------------------------------------
i.getopt <- function(x, opt = NULL) {
##  ldef <- optdefault
  if (is.null(opt))   opt <- options()
  lnam <- as.character(substitute(x))
  lx <- x
  if (is.null(lx)||(is.atomic(lx)&&all(is.na(lx))))
    lx <- opt[[lnam]]
  lx <- if (is.null(lx)||(is.atomic(lx)&&all(is.na(lx)))) getopt(lnam)
  else unlist(checkopt(lnam, lx))   ## check
  lx
}
## -----------------------------------------------------------------------
checkopt <- function(optname, value, list = NULL) {
  if (is.null(list)) list <- setNames(list(value), optname)
  lnl <- length(list)
  loptnames <- names(list)
  for (lil in seq_len(lnl)) {
    lnm <- loptnames[lil]
    lvalue <- list[[lnm]]
    lcheck <- optcheck[[lnm]]
    if (length(lcheck)) {
##    if (is.list(lcheck)) {
      if (!is.list(lcheck[[1]])) lcheck <- list(lcheck)
      lnopt <- length(lcheck)
      lmsg <- rep("", lnopt)
      for (lj in seq_len(lnopt)) {
        lch <- lcheck[[lj]]
        lfn <- get(lch[[1]])
        lmsg[lj] <- lmsgj <-
          switch(paste("v",length(lch),sep=""), v0="", v1=lfn(lvalue),
                 v2=lfn(lvalue, lch[[2]]), v3=lfn(lvalue, lch[[2]], lch[[3]]),
                 "")
        if (lmsgj=="") break
      }
##    }
      if (all(lmsg!="")) {
        warning(":check.option: argument '", lnm,
                "' not suitable. It should\n    ",
                paste(lmsg, collapse=" -- or \n  "),
                "\n  instead of (str())\n    ", format(str(lvalue)))
        list[lnm] <- NULL
      }
    }
  }
  list
}
## -----------------------------------------------------------------------
optcheck <- list()
## ========================================================================
u.isnull <- function(x)  length(x)==0
u.true <- function (x) length(x)>0 && is.logical(x) && (!is.na(x)) && all(x)
u.notfalse <- function (x) !(length(x)==1 && is.logical(x) && (!is.na(x)) && !x)
"%nin%" <- function (x,y) !x%in%y
DB <- function (on=TRUE) options(error=if(on) recover else NULL, warn=on)
## ---------------------------------------
i.def <- function(arg, value = TRUE, valuetrue = value, valuefalse = FALSE)
{
  rr <- arg
  if (length(arg)==0 ||
      (mode(arg)%in%c("numeric","character","logical","complex")&&
       all(is.na(arg)))
      )  rr <- value
  else {
    if (length(arg)==1 && is.logical(arg))
      rr <- if (arg) valuetrue else valuefalse
  }
  rr
}
## ===========================================================================
optdefault <- list(
  testlevel = 0.05,
  rlvthres = c(rel=0.1, stand=0.1, drop=0.1, pred=0.05),
  p.symbols = c("***", "**", "*", ".", " ", ""),
  p.cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
  rlv.symbols = c("-", " ", ".", "+", "++", "+++", ""),
  rlv.cutpoints = c(-Inf,-1,0,1,2,5,Inf),
  ## ----
  termtable = TRUE, vif = TRUE,
  show.termeffects = TRUE, show.coefcorr = FALSE,
  printstyle = "relevance",
  termcolumns.r = c("coef",  "df", "R2.x", "Rle", "Rls", # "Rlp",
                    "dropRle", "dropRls", "dropRls.symbol", "predRle"),
  termeffcolumns.r = c("coef","Rls.symbol"),
  coefcolumns.r = c("coef", "Rle", "Rls", "Rlp", "Rls.symbol"),
  termcolumns.p = c("coef",  "df", "ciLow","ciUp","R2.x", "signif0", "p.value",
                  "p.symbol"),
  termeffcolumns.p = c("coef","p.symbol"),
  coefcolumns.p = c("coef", "ciLow","ciUp", "signif0", "p.value", "p.symbol"),
  show.symbolLegend = TRUE
  )
optdefault$p.legend <-
  paste(rbind(optdefault$p.cutpoints, optdefault$p.symbols), collapse="  ")
optdefault$rlv.legend <-
  paste(rbind(optdefault$rlv.cutpoints, optdefault$rlv.symbols), collapse="  ")

## -----------------------------------------------------------------------
regrModelClasses <- c("regr","lm","lmrob","rlm","glm","survreg","coxph","rq","polr")
