\name{lambdamax}
\Rdversion{1.1}
\alias{lambdamax}
\alias{lambdamax.formula}
\alias{lambdamax.default}

\title{Lasso: get the Maximal lambda}
\description{
  Calculates the maximal value of the weight lambda of the L1 penalty term
  in a Lasso regression. For values >= this value, the null model will
  be obtained as the result of the penalized regression.
}
\usage{
lambdamax(x, ...)

\method{lambdamax}{formula}(formula, nonpen = ~1, data,
  weights, subset, na.action, offset, coef.init,
  penscale = sqrt, model = LogReg(),
  center = NA, standardize = TRUE, contrasts = NULL,
  nlminb.opt = list(), ...)

\method{lambdamax}{default}(x, y, index, weights = NULL,
  offset = rep(0, length(y)), coef.init = rep(0, ncol(x)),
  penscale = sqrt, model = LogReg(),
  center = NA, standardize = TRUE,
  nlminb.opt = list(), ...)
}
\arguments{
  \item{x}{design matrix (including intercept)}
  \item{y}{response vector}
  \item{formula}{\code{formula} of the penalized variables.  The response
    has to be on the left hand side of \code{~}.}
  \item{nonpen}{\code{formula} of the nonpenalized variables. This will
    be added to the \code{formula} argument above and doesn't need to have the
    response on the left hand side.}
  \item{data}{\code{data.frame} containing the variables in the model.}
  \item{index}{vector which defines the grouping of the
    variables. Components sharing the same
    number build a group. Non-penalized coefficients are marked with
    \code{NA}.}
  \item{weights}{vector of observation weights.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
  \item{offset}{vector of offset values.}
  \item{coef.init}{initial parameter vector. Penalized groups are discarded.}
  \item{penscale}{rescaling function to adjust the value of the penalty
    parameter to the degrees of freedom of the parameter group. See the
    reference below.}
  \item{model}{an object of class \code{\link{lassoModel}} implementing
    the negative log-likelihood, gradient, hessian etc. See
    \code{\link{lassoModel}} for more details.}
  \item{standardize, center}{logical; see \code{\link{lasso}}.}
  \item{contrasts}{an (optional) list with the contrasts for the factors
    in the model.}
  \item{nlminb.opt}{arguments to be supplied to \code{\link{nlminb}}.}
  \item{...}{additional arguments to be passed to the functions defined
    in \code{model}.}
}
\details{Uses \code{\link{nlminb}} to optimize the non-penalized parameters.}
\value{Numerical value of the maximal lambda}
\author{Lukas Meier, Seminar f. Statistik, ETH Zurich}
\examples{
data(splice)
lambdamax(y ~ ., data = splice, model = LogReg(), center = TRUE,
          standardize = TRUE)
}
\keyword{misc}
