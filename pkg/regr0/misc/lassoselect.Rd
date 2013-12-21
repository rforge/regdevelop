\name{lassoselect}
\alias{lassoselect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model Selection by Lasso}
\description{
  Calculate the Lasso estimates for the coefficients of a linear
  regression model and select the model with smallest cross validated
  mean square prediction error
}
\usage{
lassoselect(object, data, s = NA, cv = is.na(s), cv.k=10,
  cv.s = seq(0,1,length=100), adaptive = TRUE,
  coefficients = object$coefficients, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{either a formula or a result of an earlier call to lassoselect}
  \item{data}{data frame in which the formula is evaluated}
  \item{s}{importance of L1 penalty given as in
    \code{predict.lars (mode="fraction")}\cr
    if NA, \code{s} is determined by cross validation}
  \item{cv}{calculate cross validated mean squared prediction error
    (see  \code{\link[lars]{cv.lars}} )}
  \item{cv.k}{k-fold cross-validation is used, where k=\code{cv.k}}
  \item{cv.s}{fractions for which cross-validation ist done}
  \item{adaptive}{if TRUE, the L1 penalty is \code{sum(abs(beta/coefficients))}
    where  \code{coefficients}  is given as the next argument}
  \item{coefficients}{named vector of initial coefficients used as
    inverse weights in the L1 penalty. Only needed if adaptive is
    TRUE.\crIf \code{NULL}, \code{lassoselect} is called with
    \code{adaptive=FALSE} to get these preliminary coefficients.
  }
  \item{\dots}{arguments passed to \code{cv.lars} }
}
%%- \details{
%%-   adaptive ...
%%- }
\value{
  List consisting of the components
  \item{lars}{result of the call to \code{lars}, see \code{\link[lars]{lars}}}
  \item{data}{the data used}
  \item{formula}{the formula used}
  \item{s}{the fraction value for which the optimal model results}
  \item{coefficients}{the coefficients of the optimal model}
  \item{fitted}{the fitted values for the optimal model}
  \item{adaptive}{indicates if adaptive lasso has been used}
  \item{cv}{data.frame containing the results of cross validation\cr
    with columns \code{fraction, cv, cv.error} }
}
\references{see \code{?lars} and ... for adaptive lasso}
\author{Werner A. Stahel, ETH Zurich}
\note{
}
\seealso{\code{\link[lars]{lars}}}
\examples{
## data(d.blast)
## t.r <- lassoselect(log10(tremor)~location+log10(distance)+log10(charge),
##                    data=d.blast)
## plot(t.r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
