library(lassogrp)

data(splice)
tol <- 5 * 10^-8
control <- lassoControl(tol = tol, trace = 0)

### ---- Check  that lambdamax()  *does* give the *maximal* lambda:
lambda.max <- lambdamax(y ~ ., data = splice, model = LogReg())
stopifnot(all.equal(lambda.max, 68.30413602469, tol = 1e-12))
## is lambdamax() correct?
fit. <- lasso(y ~ ., data = splice, model = LogReg(),
	      lambda = c(10, 1.1, 1.01, 1, 0.99, 0.5) * lambda.max,
	      control = control)
stopifnot(apply(coef(fit.)[-1, ], # -1 : drop intercept
		2, function(.) all(. == 0))
	  == c(TRUE, TRUE, TRUE,TRUE, FALSE,FALSE))
##  lambda = c( 10,  1.1,  1.01, 1,   0.99,  0.5) * lambda.max


##############################################
##                                          ##
## See whether correct solution is obtained ##
##                                          ##
##############################################

C.12.sum <- list(Pos.1 = "contr.sum",
                 Pos.2 = "contr.sum")
fit <- lasso(y ~ Pos.1 * Pos.2, data = splice, model = LogReg(), lambda = 25, control = control,
             center = FALSE, contrast = C.12.sum)

p0 <- function(...) paste(..., sep="")
m <- outer(p0("Pos.",1:2),1:3, p0)
rn <- c("(Intercept)", c(t(m), outer(m[1,], m[2,], paste, sep=":")))
sol <- structure(c(-0.13977233, 0.0226494585308358, 0.0897927902302861,
-0.0311859296853295, 0.519851153317134, -0.25546398114419, -0.209672224889056,
0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(16L, 1L), .Dimnames = list(
                                               rn, "l.1"))

stopifnot(all.equal(cf <- coef(fit), sol, tol = 1e-6),
          identical(cf == 0, sol == 0))

## here, 'center' is default --> TRUE, as we have intercept
fit1C <- lasso(y ~ Pos.1 * Pos.2, data = splice, model = LogReg(), lambda = 25, control = control,
               contrast = C.12.sum)
## FIXME: "compare" somehow with 'fit' above

##################################################################
##                                                              ##
## Check whether different contrasts lead to the same solutions ##
## I.e. check whether the (back-) transformations work right    ##
##                                                              ##
##################################################################

contr.A <- list(Pos.1 = "contr.sum",
                Pos.2 = "contr.sum",
                Pos.3 = "contr.sum")
contr.B <- list(Pos.1 = "contr.helmert",
                Pos.2 = "contr.helmert",
                Pos.3 = "contr.helmert")

fit.A <- lasso(y ~ Pos.1 * Pos.2 * Pos.3, nonpen = ~ 1, data = splice, model = LogReg(),
                  standardize = TRUE, lambda = 61:1,
                  contrasts = contr.A, control = control)
fit.B <- lasso(y ~ Pos.1 * Pos.2 * Pos.3, nonpen = ~ 1, data = splice, model = LogReg(),
                  standardize = TRUE, lambda = 61:1,
                  contrasts = contr.B, control = control)

op <- par(mfrow = c(1, 2), mgp = c(1.5, 0.6, 0), mar= .1+c(4,4,2,1))
plot(fit.A, log = "x")
plot(fit.B, log = "x")
par(op)

if(max(abs(fit.A$fn.val - fit.B$fn.val) / fit.A$fn.val) > tol)
  stop("Inconsistent result when changing the encoding scheme (fn.val)")

pred.A <- predict(fit.A, newdata = splice, type = "response")
pred.B <- predict(fit.B, newdata = splice, type = "response")
m      <- abs(pred.A - pred.B)

if(max(m) > sqrt(tol))
  stop("Inconsistent result when changing the encoding scheme (prediction)")

#range(pred.A - pred.B)
#range((pred.A - pred.B) / (1 + pred.A))

###############################################
##                                           ##
## Check whether offset is working correctly ##
##                                           ##
###############################################

## Fit an ordinary model, with unpenalized intercept
fit1 <- lasso(y ~ ., data = splice, model = LogReg(), lambda = 0.5 * lambda.max,
              standardize = TRUE, control = control)

## Plugging in the intercept as an offset should lead to an intercept
## which is close to zero
shift <- 4
intercept <- rep(shift, nrow(splice))

try( ## this currently fails  (in eval(.) inside model.frame.default()... AAARGH!
fit2. <- lasso(y ~ . + offset(intercept), lambda = 0.5 * lambda.max,
               data = splice, model = LogReg(), standardize = TRUE, control = control)
)
fit2 <- lasso(y ~ ., offset = intercept, lambda = 0.5 * lambda.max,
              data = splice, model = LogReg(), standardize = TRUE, control = control)
d.coef    <- coef(fit2) - coef(fit1)
d.coef[1] <- d.coef[1] + shift

if(max(abs(d.coef) / (1 + abs(coef(fit1)))) > 5 * 10^-4)
  stop("Inconsistent result when using offset (d.coef)")

#################################################
##                                             ##
## Check whether max.iter is working correctly ##
##                                             ##
#################################################

fit.maxiter <-
  lasso(y ~ Pos.1 * Pos.2, data = splice, model = LogReg(), lambda = c(1, 0.1),
        control = lassoControl(max.iter = 2, trace = 2, inner.loops = 0),
        contrast = list(Pos.1 = "contr.sum", Pos.2 = "contr.sum"))

stopifnot(all(!fit.maxiter$converged))
