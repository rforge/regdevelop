## drop1Wald
## vif.regr

## ===========================================================================
drop1Wald <-
  function (object, scope=NULL, scale = 0, test = c("none", "Chisq", "F"),
           k = 2, ...)
{
    x <- model.matrix(object)
    ## offset <- model.offset(model.frame(object))
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
      if (length(ii)) {
        RSS[i] <- if (length(ii)==1) coef[ii]^2/cov[ii,ii] else
          coef[ii]%*%solve(cov[ii,ii])%*%coef[ii]  ## !!! REPLACE THIS
        dfs[i] <- length(ii)
      } else {
        dfs[i] <- 0
        RSS[i] <- NA
      }
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

## ==========================================================================
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
