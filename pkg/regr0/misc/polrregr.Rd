\name{polrregr}
\Rdversion{1.1}
\alias{polrregr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
polrregr(formula, data, weights, start, ..., subset, na.action, contrasts = NULL, Hess = FALSE, model = TRUE, method = c("logistic", "probit", "cloglog", "cauchit"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
}
  \item{contrasts}{
%%     ~~Describe \code{contrasts} here~~
}
  \item{Hess}{
%%     ~~Describe \code{Hess} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data, weights, start, ..., subset, na.action, 
    contrasts = NULL, Hess = FALSE, model = TRUE, method = c("logistic", 
        "probit", "cloglog", "cauchit"))
  ## copy of polr from MASS. 1 line added by WSt to keep eta in the result
{
    logit <- function(p) log(p/(1 - p))
    fmin <- function(beta) {
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 
            100)
        eta <- offset
        if (pc > 0) 
            eta <- eta + drop(x \%*\% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        if (all(pr > 0)) 
            -sum(wt * log(pr))
        else Inf
    }
    gmin <- function(beta) {
        jacobian <- function(theta) {
            k <- length(theta)
            etheta <- exp(theta)
            mat <- matrix(0, k, k)
            mat[, 1] <- rep(1, k)
            for (i in 2:k) mat[i:k, i] <- etheta[i]
            mat
        }
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 
            100)
        eta <- offset
        if (pc > 0) 
            eta <- eta + drop(x \%*\% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        p1 <- dfun(gamm[y + 1] - eta)
        p2 <- dfun(gamm[y] - eta)
        g1 <- if (pc > 0) 
            t(x) \%*\% (wt * (p1 - p2)/pr)
        else numeric(0)
        xx <- .polrY1 * p1 - .polrY2 * p2
        g2 <- -t(xx) \%*\% (wt/pr)
        g2 <- t(g2) \%*\% jacobian(theta)
        if (all(pr > 0)) 
            c(g1, g2)
        else rep(NA, pc + q)
    }
    m <- match.call(expand.dots = FALSE)
    method <- match.arg(method)
    pgumbel <- function(q) exp(pweibull(log(q))) # ???
    pfun <- switch(method, logistic = plogis, probit = pnorm, 
        cloglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logistic = dlogis, probit = dnorm, 
        cloglog = dgumbel, cauchit = dcauchy)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$start <- m$Hess <- m$method <- m$model <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    n <- nrow(x)
    pc <- ncol(x)
    cons <- attr(x, "contrasts")
    if (xint > 0) {
        x <- x[, -xint, drop = FALSE]
        pc <- pc - 1
    }
    else warning("an intercept is needed and assumed")
    wt <- model.weights(m)
    if (!length(wt)) 
        wt <- rep(1, n)
    offset <- model.offset(m)
    if (length(offset) <= 1) 
        offset <- rep(0, n)
    y <- model.response(m)
    if (!is.factor(y)) 
        stop("response must be a factor")
    lev <- levels(y)
    if (length(lev) <= 2) 
        stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- length(lev) - 1
    Y <- matrix(0, n, q)
    .polrY1 <- col(Y) == y
    .polrY2 <- col(Y) == y - 1
    if (missing(start)) {
        q1 <- length(lev)\%/\%2
        y1 <- (y > q1)
        X <- cbind(Intercept = rep(1, n), x)
        fit <- switch(method, logistic = glm.fit(X, y1, wt, family = binomial(), 
            offset = offset), probit = glm.fit(X, y1, wt, family = binomial("probit"), 
            offset = offset), cloglog = glm.fit(X, y1, wt, family = binomial("probit"), 
            offset = offset), cauchit = glm.fit(X, y1, wt, family = binomial("cauchit"), 
            offset = offset))
        if (!fit$converged) 
            stop("attempt to find suitable starting values failed")
        coefs <- fit$coefficients
        if (any(is.na(coefs))) {
            warning("design appears to be rank-deficient, so dropping some coefs")
            keep <- names(coefs)[!is.na(coefs)]
            coefs <- coefs[keep]
            x <- x[, keep[-1], drop = FALSE]
            pc <- ncol(x)
        }
        spacing <- logit((1:q)/(q + 1))
        if (method != "logistic") 
            spacing <- spacing/1.7
        gammas <- -coefs[1] + spacing - spacing[q1]
        thetas <- c(gammas[1], log(diff(gammas)))
        s0 <- c(coefs[-1], thetas)
    }
    else if (length(start) != pc + q) 
        stop("'start' is not of the correct length")
    else {
        s0 <- if (pc > 0) 
            c(start[seq_len(pc + 1)], diff(start[-seq_len(pc)]))
        else c(start[1], diff(start))
    }
    res <- optim(s0, fmin, gmin, method = "BFGS", hessian = Hess, 
        ...)
    beta <- res$par[seq_len(pc)]
    theta <- res$par[pc + 1:q]
    zeta <- cumsum(c(theta[1], exp(theta[-1])))
    deviance <- 2 * res$value
    niter <- c(f.evals = res$counts[1], g.evals = res$counts[2])
    names(zeta) <- paste(lev[-length(lev)], lev[-1], sep = "|")
    if (pc > 0) {
        names(beta) <- colnames(x)
        eta <- drop(x \%*\% beta)
    }
    else {
        eta <- rep(0, n)
    }
    cumpr <- matrix(pfun(matrix(zeta, n, q, byrow = TRUE) - eta), 
        , q)
    fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
    dimnames(fitted) <- list(row.names(m), lev)
    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance,
                linear.predictor = eta,  ## added by WSt
        fitted.values = fitted, lev = lev, terms = Terms, df.residual = sum(wt) - 
            pc - q, edf = pc + q, n = sum(wt), nobs = sum(wt), 
        call = match.call(), method = method, convergence = res$convergence, 
        niter = niter)
    if (Hess) {
        dn <- c(names(beta), names(zeta))
        H <- res$hessian
        dimnames(H) <- list(dn, dn)
        fit$Hessian <- H
    }
    if (model) 
        fit$model <- m
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    class(fit) <- "polr"
    fit
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
