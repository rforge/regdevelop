#### WS-tools : Some of Werner Stahel's R Tools, not related much to the package:
#### ========

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
## ---
stamp <- function(sure=TRUE, outer.margin = NULL,
                  project=getOption("project"), step=getOption("step"),
                  stamp=getOption("stamp"), ...)
{
  ## Purpose:   plot date and project information
  ## -------------------------------------------------------------------------
  ## Arguments:
  ##   sure     if F, the function only plots its thing if  getOption("stamp")>0
  ##   outer    if T, the date is written in the outer margin
  ##   project  project title
  ##   step     title of step of data analysis
  ##   ...      arguments to  mtext , e.g., line=3
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 13 Aug 96, 09:00
  if (length(stamp)==0) {
    message(":stamp: setting options(stamp=1)")
    options(stamp=1)
    stamp <- 1
  }
  if (length(outer.margin)==0) outer.margin <- par('oma')[4]>0
  t.txt <- date()
  t.txt <- paste(substring(t.txt,5,10),",",substring(t.txt,22,23),"/",
                 substring(t.txt,13,16),sep="")
  if (length(project)>0) t.txt <- paste(t.txt,project,sep=" | ")
  if (length(step)>0) t.txt <- paste(t.txt,step,sep=" | ")
  if (sure || stamp==2 || ( stamp==1 && (
                              ##     last figure on page
                              { t.tfg <- par("mfg") ; all(t.tfg[1:2]==t.tfg[3:4]) }
                              || isTRUE(outer.margin) ))  )
    mtext(t.txt, 4, cex = 0.6, adj = 0, outer = outer.margin, ...)
  invisible(t.txt)
}
## =======================================================================
factoreffects <- function (object, se = 2, # use.na = TRUE, 
                      df = df.residual(object), ...)  {
  if (is.atomic(object)||is.null(terms(object)))
      stop("!factoreffects! inadequate first argument")
 ##  xl <- object$xlevels
  Terms <- delete.response(terms(object))
  tl <- attr(Terms, "term.labels")
  dcl <- attr(Terms,"dataClasses")[-1]
  if (all(dcl=="numeric")) 
    return(as.list(coef(object)))
  ## result already available?
  faceff <- object$faceff
  if ((!is.null(faceff))&&length(faceff)==length(tl)&&
      (is.matrix(faceff[[length(faceff)]])|!se)) return(faceff) ## !!! check!
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
  mm <- model.matrix(Terms, df.dummy, contrasts.arg=lcontr, xlev=xtlv) ## 
  if (anyNA(mm)) {
    warning("some terms will have NAs due to the limits of the method")
    mm[is.na(mm)] <- NA
  }
  ## calculate dummy coefs
  coef <- object$coefficients ##!!! cf <-
##-   if (!use.na) 
##-     coef[is.na(coef)] <- 0
  if (se) {
    cov <- vcov(object)
    if (is.null(cov)) {
      warning(":factoreffects: no covariance matrix of coefficients found.",
              " Returning coefficients only")
      se <- FALSE
    }
  }
  asgn <- attr(mm, "assign")
  names(asgn) <- colnames(mm)
  asgn <- asgn[names(coef)] ## !!!
  res <- setNames(vector("list", length(tl)), tl)
  ljfail <- NULL
  for (j in seq_along(tl)) {
    mmr <- rn == tl[j]  ## rows corresponding to the term
    mmc <- names(asgn)[asgn == j & !is.na(coef)]  ## columns (logical fails for polr, vcov() too large) !!! was  which
    mmpart <- mm[mmr, mmc, drop=FALSE]
    rrj <- setNames(drop(mmpart %*% coef[mmc]), rnn[mmr])
    if (se) {
      if (any(is.na(rrj)))
        warning(":factoreffects: missing coef for term '", tl[j],
                "'. no standard errors etc") else {
##-      if (any(is.na(licov))) ljfail <- c(ljfail,j)  else {
          sej <- sqrt(diag(mmpart %*% cov[mmc,mmc] %*% t(mmpart)))
          rrj <- ciSignif(rrj, sej, df)
        }
    }
    res[[j]] <- rrj
  }
  if (length(ljfail))
    warning(":factoreffects: error calculating se for terms ",
            paste(ljfail, collapse=", "))
  if (int > 0) {
    res <- c(list(`(Intercept)` = coef[int]), res)
  }
  class(res) <- "faceff"
  res
}
## --------------------------------------------------------------------
print.faceff <- function(x, columns=NULL,  transpose=FALSE, ...)
{
  if (is.null(columns)) columns <- "all"
  columns[columns=="coef"] <- "estimate"
  csymb <- "coefsymb"%in%columns
  if ("all"%in%columns)  columns <-
      if(csymb)
        c("coefsymb", "se", "ciLow", "ciHigh", "testst",
          "signif", "p.value") else
        c("estimate", "se", "ciLow", "ciHigh", "testst",
          "signif", "p.value")
  for (li in seq_along(x)) {
    xi <- x[[li]]
    if (is.null(dim(xi))) next
    if (csymb)
      xi$coefsymb <-
        if ("p.symb"%in%names(xi))
          paste(format(xi[,1],...), format(xi[,"p.symb"])) else  xi[,1]
    xif <- format(xi[,intersect(columns,names(xi)), drop=FALSE],...)
    xif <- if (ncol(xif)==1 || (nrow(xif)>1 & transpose)) t(xif) else xif
    if (nrow(xif)==1) row.names(xif) <- " "  ## drop row name
    if (ncol(xif)==1) colnames(xif) <- " "  ## drop col name
    if (prod(dim(xif))==1) xif <- as.character(xif[1,1])
    x[li] <- list(xif)
  }
  class(x) <- NULL
  print.default(x, quote=FALSE, ...)
}
## -------------------------------------------------------------------------
ciSignif <- function(estimate, se=NULL, df=Inf, testlevel=0.05) {
  if (is.null(se))
    if (NCOL(estimate)>1) {
      se <- estimate[,2]
      estimate <- estimate[,1]
    } else
      stop("!ciSignif! no standard errors found")
  ltq <- qt(1-testlevel/2, df)
  lci <- estimate+outer(ltq*se, c(ciLow=-1,ciHigh=1))
  ltst <- estimate/se
  lsgf <- ltst/ltq
  lpv <- 2*pt(-abs(ltst), df)
  lipv <- as.numeric(cut(lpv, c(0, 0.001, 0.01, 0.05, 0.1, 1)))
  lsst <- c("***", "**", "*", ".", " ")[lipv]
  data.frame(estimate=estimate, se=se, lci, testst=ltst,
             signif=lsgf, p.value=lpv, p.symb=lsst)
}
