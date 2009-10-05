if(getRversion() <= "2.7.1") {

add1.default <- function(object, scope, scale = 0, test=c("none", "Chisq"),
			 k = 2, trace = FALSE, ...)
{
    if(missing(scope) || is.null(scope)) stop("no terms in scope")
    if(!is.character(scope))
	scope <- add.scope(object, update.formula(object, scope))
    if(!length(scope))
	stop("no terms in scope for adding to object")
#     newform <- update.formula(object,
#                               paste(". ~ . +", paste(scope, collapse="+")))
#     data <- model.frame(update(object, newform)) # remove NAs
#     object <- update(object, data = data)
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2,
                  dimnames = list(c("<none>", scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
    n0 <- length(object$residuals)
    env <- environment(formula(object))
    for(i in seq(ns)) {
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying +", tt, "\n", sep='')
	    utils::flush.console()
	}
	nfit <- update(object, as.formula(paste("~ . +", tt)),
                       evaluate = FALSE)
	nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)
	ans[i+1, ] <- extractAIC(nfit, scale, k = k, ...)
        if(length(nfit$residuals) != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[,1] - ans[1,1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[,2])
    test <- match.arg(test)
    if(test == "Chisq") {
	dev <- ans[,2] - k*ans[, 1]
	dev <- dev[1] - dev; dev[1] <- NA
	nas <- !is.na(dev)
	P <- dev
	P[nas] <- pchisq(dev[nas], dfs[nas], lower.tail=FALSE)
	aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    head <- c("Single term additions", "\nModel:",
	      deparse(as.vector(formula(object))),
	      if(scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

drop1.default <- function(object, scope, scale = 0, test=c("none", "Chisq"),
			  k = 2, trace = FALSE, ...)
{
    tl <- attr(object$terms, "term.labels")
    if(missing(scope)) scope <- drop.scope(object)
    else {
	if(!is.character(scope))
	    scope <- attr(terms(update.formula(object, scope)), "term.labels")
	if(!all(match(scope, tl, 0L) > 0L))
	    stop("scope is not a subset of term labels")
    }
#    data <- model.frame(object) # remove NAs
#    object <- update(object, data = data)
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2,
                  dimnames =  list(c("<none>", scope), c("df", "AIC")))
    ans[1, ] <- extractAIC(object, scale, k = k, ...)
    n0 <- length(object$residuals)
    env <- environment(formula(object))
    for(i in seq(ns)) {
	tt <- scope[i]
	if(trace > 1) {
	    cat("trying -", tt, "\n", sep='')
	    utils::flush.console()
        }
        nfit <- update(object, as.formula(paste("~ . -", tt)),
                       evaluate = FALSE)
	nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)
	ans[i+1, ] <- extractAIC(nfit, scale, k = k, ...)
        if(length(nfit$residuals) != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[1,1] - ans[,1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[,2])
    test <- match.arg(test)
    if(test == "Chisq") {
        dev <- ans[, 2] - k*ans[, 1]
        dev <- dev - dev[1] ; dev[1] <- NA
        nas <- !is.na(dev)
        P <- dev
        P[nas] <- pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
        aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    head <- c("Single term deletions", "\nModel:",
	      deparse(as.vector(formula(object))),
	      if(scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

drop1.lm <-
  function (object, scope, scale = 0, all.cols = TRUE, test = c("none", 
    "Chisq", "F"), k = 2, ...) 
{
    x <- model.matrix(object)
    offset <- model.offset(model.frame(object))
    iswt <- !is.null(wt <- object$weights)
    n <- nrow(x)
    asgn <- attr(x, "assign")
    tl <- attr(object$terms, "term.labels")
    if (missing(scope)) 
        scope <- drop.scope(object)
    else {
        if (!is.character(scope)) 
            scope <- attr(terms(update.formula(object, scope)), 
                "term.labels")
        if (!all(match(scope, tl, 0L) > 0L)) 
            stop("scope is not a subset of term labels")
    }
    ndrop <- match(scope, tl)
    ns <- length(scope)
    rdf <- object$df.residual
    chisq <- deviance.lm(object)
    dfs <- numeric(ns)
    RSS <- numeric(ns)
##-     y <- object$residuals + predict(object)
    y <- object$residuals + object$fitted.values
    na.coef <- (1:length(object$coefficients))[!is.na(object$coefficients)]
    for (i in 1:ns) {
        ii <- seq_along(asgn)[asgn == ndrop[i]]
        if (all.cols) 
            jj <- setdiff(seq(ncol(x)), ii)
        else jj <- setdiff(na.coef, ii)
        z <- if (iswt) 
            lm.wfit(x[, jj, drop = FALSE], y, wt, offset = offset)
        else lm.fit(x[, jj, drop = FALSE], y, offset = offset)
        dfs[i] <- z$rank
        oldClass(z) <- "lm"
        RSS[i] <- deviance(z)
    }
    scope <- c("<none>", scope)
    dfs <- c(object$rank, dfs)
    RSS <- c(chisq, RSS)
    if (scale > 0) 
        aic <- RSS/scale - n + k * dfs
    else aic <- n * log(RSS/n) + k * dfs
    dfs <- dfs[1] - dfs
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, "Sum of Sq" = c(NA, RSS[-1] - 
        RSS[1]), RSS = RSS, AIC = aic, row.names = scope, check.names = FALSE)
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
    head <- c("Single term deletions", "\nModel:", deparse(as.vector(formula(object))), 
        if (scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
# <environment: namespace:stats>
deviance.lm <- function (object, ...) 
sum(weighted.residuals(object)^2, na.rm = TRUE)
}## only for R version <= 2.7.1
##- 
##- step <- 
##- function (object, scope, scale = 0, direction = c("both", "backward", 
##-     "forward"), trace = 1, keep = NULL, steps = 1000, k = 2, 
##-     ...) 
##- {
##-     mydeviance <- function(x, ...) {
##-         dev <- deviance(x)
##-         if (!is.null(dev)) 
##-             dev
##-         else extractAIC(x, k = 0)[2]
##-     }
##-     cut.string <- function(string) {
##-         if (length(string) > 1) 
##-             string[-1] <- paste("\n", string[-1], sep = "")
##-         string
##-     }
##-     re.arrange <- function(keep) {
##-         namr <- names(k1 <- keep[[1]])
##-         namc <- names(keep)
##-         nc <- length(keep)
##-         nr <- length(k1)
##-         array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, 
##-             namc))
##-     }
##-     step.results <- function(models, fit, object, usingCp = FALSE) {
##-         change <- sapply(models, "[[", "change")
##-         rd <- sapply(models, "[[", "deviance")
##-         dd <- c(NA, abs(diff(rd)))
##-         rdf <- sapply(models, "[[", "df.resid")
##-         ddf <- c(NA, diff(rdf))
##-         AIC <- sapply(models, "[[", "AIC")
##-         heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
##-             "\nInitial Model:", deparse(as.vector(formula(object))), 
##-             "\nFinal Model:", deparse(as.vector(formula(fit))), 
##-             "\n")
##-         aod <- data.frame(Step = I(change), Df = ddf, Deviance = dd, 
##-             "Resid. Df" = rdf, "Resid. Dev" = rd, AIC = AIC, 
##-             check.names = FALSE)
##-         if (usingCp) {
##-             cn <- colnames(aod)
##-             cn[cn == "AIC"] <- "Cp"
##-             colnames(aod) <- cn
##-         }
##-         attr(aod, "heading") <- heading
##-         fit$anova <- aod
##-         fit
##-     }
##-     Terms <- terms(object)
##-     object$call$formula <- object$formula <- Terms
##-     md <- missing(direction)
##-     direction <- match.arg(direction)
##-     backward <- direction == "both" | direction == "backward"
##-     forward <- direction == "both" | direction == "forward"
##-     if (missing(scope)) {
##-         fdrop <- numeric(0)
##-         fadd <- attr(Terms, "factors")
##-         if (md) 
##-             forward <- FALSE
##-     }
##-     else {
##-         if (is.list(scope)) {
##-             fdrop <- if (!is.null(fdrop <- scope$lower)) 
##-                 attr(terms(update.formula(object, fdrop)), "factors")
##-             else numeric(0)
##-             fadd <- if (!is.null(fadd <- scope$upper)) 
##-                 attr(terms(update.formula(object, fadd)), "factors")
##-         }
##-         else {
##-             fadd <- if (!is.null(fadd <- scope)) 
##-                 attr(terms(update.formula(object, scope)), "factors")
##-             fdrop <- numeric(0)
##-         }
##-     }
##-     models <- vector("list", steps)
##-     if (!is.null(keep)) 
##-         keep.list <- vector("list", steps)
##-     n <- length(object$residuals)
##-     fit <- object
##-     bAIC <- extractAIC(fit, scale, k = k, ...)
##-     edf <- bAIC[1]
##-     bAIC <- bAIC[2]
##-     if (is.na(bAIC)) 
##-         stop("AIC is not defined for this model, so 'step' cannot proceed")
##-     nm <- 1
##-     Terms <- fit$terms
##-     if (trace) 
##-         cat("Start:  AIC=", format(round(bAIC, 2)), "\n", cut.string(deparse(as.vector(formula(fit)))), 
##-             "\n\n", sep = "")
##-     models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
##-         edf, change = "", AIC = bAIC)
##-     if (!is.null(keep)) 
##-         keep.list[[nm]] <- keep(fit, bAIC)
##-     usingCp <- FALSE
##-     while (steps > 0) {
##-         steps <- steps - 1
##-         AIC <- bAIC
##-         ffac <- attr(Terms, "factors")
##-         scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
##-         aod <- NULL
##-         change <- NULL
##-         if (backward && length(scope$drop)) {
##-             aod <- drop1(fit, scope$drop, scale = scale, trace = trace, 
##-                 k = k, ...)
##-             rn <- row.names(aod)
##-             row.names(aod) <- c(rn[1], paste("-", rn[-1], sep = " "))
##-             if (any(aod$Df == 0, na.rm = TRUE)) {
##-                 zdf <- aod$Df == 0 & !is.na(aod$Df)
##-                 change <- rev(rownames(aod)[zdf])[1]
##-             }
##-         }
##-         if (is.null(change)) {
##-             if (forward && length(scope$add)) {
##-                 aodf <- add1(fit, scope$add, scale = scale, trace = trace, 
##-                   k = k, ...)
##-                 rn <- row.names(aodf)
##-                 row.names(aodf) <- c(rn[1], paste("+", rn[-1], 
##-                   sep = " "))
##-                 aod <- if (is.null(aod)) 
##-                   aodf
##-                 else rbind(aod, aodf[-1, , drop = FALSE])
##-             }
##-             attr(aod, "heading") <- NULL
##-             nzdf <- if (!is.null(aod$Df)) 
##-                 aod$Df != 0 | is.na(aod$Df)
##-             aod <- aod[nzdf, ]
##-             if (is.null(aod) || ncol(aod) == 0) 
##-                 break
##-             nc <- match(c("Cp", "AIC"), names(aod))
##-             nc <- nc[!is.na(nc)][1]
##-             o <- order(aod[, nc])
##-             if (trace) 
##-                 print(aod[o, ])
##-             if (o[1] == 1) 
##-                 break
##-             change <- rownames(aod)[o[1]]
##-         }
##-         usingCp <- match("Cp", names(aod), 0) > 0
##-         fit <- update(fit, paste("~ .", change), evaluate = FALSE)
##-         fit <- eval.parent(fit)
##-         if (length(fit$residuals) != n) 
##-             stop("number of rows in use has changed: remove missing values?")
##-         Terms <- terms(fit)
##-         bAIC <- extractAIC(fit, scale, k = k, ...)
##-         edf <- bAIC[1]
##-         bAIC <- bAIC[2]
##-         if (trace) 
##-             cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n", 
##-                 cut.string(deparse(as.vector(formula(fit)))), 
##-                 "\n\n", sep = "")
##-         BR()
##-         if (bAIC >= AIC + 1e-07) 
##-             break
##-         nm <- nm + 1
##-         models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - 
##-             edf, change = change, AIC = bAIC)
##-         if (!is.null(keep)) 
##-             keep.list[[nm]] <- keep(fit, bAIC)
##-     }
##-     if (!is.null(keep)) 
##-         fit$keep <- re.arrange(keep.list[seq(nm)])
##-     step.results(models = models[seq(nm)], fit, object, usingCp)
##- }
