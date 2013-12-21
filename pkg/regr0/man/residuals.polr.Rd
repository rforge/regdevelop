\name{residuals.polr}
\alias{residuals.polr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residuals of a Binary or Ordered Regression}
\description{
  Calculates quartiles and random numbers according to the
  conditional distribution of residuals for the latent variable of a
  binary or ordinal regression, given the observed response value.
  See Details for an explanation.
}
\usage{
\method{residuals}{polr}(object, na.action = object, ...)
}
\arguments{
  \item{object}{the result of \code{polr} or of
    \code{glm(,family=binomial)} with binary data.}
  \item{na.action}{treatment of NAs}
  \item{...}{unused}
}
\details{
  For binary and ordinal regression, the regression models can be
  described by introducing a latent response variable Z of which the
  observed response Y is a classified version, and for which a linear
  regression applies. The errors of this "latent regression" have a
  logistic distribution. Given the linearly predicted value eta[i],
  which is the fitted value for the latent variable, the residual for
  Z[i] can therefore be assumed to have a logistic distribution.

  This function calculates quantiles and random numbers according to the
  conditional distribution of residuals for Z[i], given the observed
  y[i].
}
\value{
  a data.frame with the variables 
  \item{median}{medians of the conditional distributions}
  \item{lowq}{lower quartiles}
  \item{uppq}{upper quartiles}
  \item{random}{random numbers, drawn according to the conditional
    distributions}
  \item{fit}{linear predictor values}
  \item{y}{observed response values}
}
\references{See http://stat.ethz.ch/~stahel/regression}
\author{Werner A. Stahel, ETH Zurich}
%%- \note{ ~~further notes~~ 
%%- 
%%-  ~Make other sections like Warning with \section{Warning }{....} ~
%%- }
\seealso{\code{\link{condquant}}, \code{\link{plot.regr}}}
\examples{
house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
house.resid <- residuals(house.plr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
