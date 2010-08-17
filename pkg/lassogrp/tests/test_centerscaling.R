library(lassogrp)

rX <- function(n,p) runif(n * p, min = -2.5, max = 2.5)

set.seed(7198)

n <- 100
p <- 10
x <- matrix(rX(n,p), nrow = n, ncol = p)
fx <- 4 * sin(x[,1]) + x[,2]^2
y <- fx + rnorm(n)

x.new <- matrix(rX(n,p), nrow = n, ncol = p)

shift <- 10
scale <- 10

## Use Lasso shift, i.e. group-size = 1
fit <- lasso(cbind(1,x), y, index = c(NA, 1:p),
             lambda = c(50, 10), center = TRUE, standardize = TRUE,
             model = LinReg(), control = lassoControl(trace = 0))

## Rescale
x.resc <- shift + scale * x
fit.resc <- lasso(cbind(1,x.resc), y, index = c(NA, 1:10),
                  lambda = c(50, 10), center = TRUE, standardize = TRUE,
                  model = LinReg(), control = lassoControl(trace = 0))

## Compare estimators, without intercept
stopifnot(all.equal(coef(fit)[-1,], coef(fit.resc)[-1,] * scale))

## Check intercepts
mu.x <- colMeans(x)
int  <- mean(y) - colSums(coef(fit)[-1,] * mu.x)

mu.x.resc <- colMeans(x.resc)
int.resc  <- mean(y) - colSums(coef(fit.resc)[-1,] * mu.x.resc)

stopifnot(all.equal(int, coef(fit)[1,], tol = 1e-7))
stopifnot(all.equal(int.resc, coef(fit.resc)[1,], tol = 1e-7))

## Compare predictions
stopifnot(all.equal(predict(fit, newdata = cbind(1,x.new)),
                    predict(fit.resc,
                            newdata = cbind(1, shift + scale * x.new))))

## Check whether every case is running, including function lambda.max

## center = TRUE & unpenalized intercept

x.use <- cbind(1, x)
index <- c(0:10)

lambda1 <- c(1, 0.1) * lambdamax(x.use, y, index, model = LinReg())
fit1 <- lasso(x.use, y, index, model = LinReg(), lambda = lambda1,
              center = TRUE)
fit1

## center = TRUE & penalized intercept
x.use <- cbind(1, x)
index <- c(99, 1:10) # "99" !
lambda2 <- c(1, 0.1) * lambdamax(x.use, y, index, model = LinReg())
fit2 <- lasso(x.use, y, index, model = LinReg(), lambda = lambda2,
              center = TRUE)
##--> two warnings -- FIXME ??

## center = TRUE & *no* intercept -- gives a warning
index <- 1:10
(lambda3 <- c(1, 0.1) * lambdamax(x, y, index, model = LinReg()))
fit3 <- lasso(x, y, index, model = LinReg(), lambda = lambda3,
              center = TRUE)

## center = FALSE & *no* intercept
fit4 <- lasso(x, y, index, model = LinReg(), lambda = lambda3,
              center = FALSE)

## NOTE:  fit3 and fit4 are NO LONGER the same
##     {with Lukas' grplasso, they were, as he had
##		if(!has.interc) center <- FALSE  }
## stopifnot(all.equal(coef(fit3), coef(fit4)))


## center = FALSE & standardize = TRUE -- + intercept
x.use <- cbind(1,x)
index <- c(0:10)
lambda5 <- c(1, 0.1) * lambdamax(x.use, y, index, model = LinReg())
fit5 <- lasso(x.use, y, index, model = LinReg(), lambda = lambda5,
              center = FALSE)
## -> warning (which is ok)
## should be close  to 'fit1' as 'x' is typically quite similar to centered x
fit5
