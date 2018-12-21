\name{residuals.regrpolr}
\alias{residuals.regrpolr}
\alias{residuals.regrsurvreg}
\alias{residuals.regrcoxph}
% \alias{residuals.polr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residuals of a Binary or Ordered Regression}
\description{
  Calculates quartiles and random numbers according to the
  conditional distribution of residuals for the latent variable of a
  binary or ordinal regression, given the observed response value.
  See Details for an explanation.
}
\usage{
\method{residuals}{regrpolr}(object, type="condquant", ...)
%% \method{residuals}{polr}(object, type=c("condquant","deviance"), ...)
\method{residuals}{regrsurvreg}(object, type="condquant", ...)
\method{residuals}{regrcoxph}(object, type="CoxSnellMod", ...)
}
\arguments{
  \item{object}{the result of \code{polr} or of
    \code{glm(,family=binomial)} with binary data.}
%%-   \item{na.action}{a list with component \code{na.action}, which
%%-     determines the na.action for the residual component of
%%\code{object} }
\item{type}{type of residuals:
  \code{"condquant"} for response type residuals augmented, for censored
  observations, by conditional quantiles (see
  \code{\link{condquant}}). 
%% or \code{"deviance"} for deviance residuals\cr
  For \code{residuals.regrsurvreg}, type \code{CoxSnellMod} yields
  a modified version of Cox-Snell residuals, also including
  a \code{condquant} attribute, see Details.
  }
  \item{...}{arguments passed to standard methods of \code{residuals}}
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

  Modified Cox-Snell residuals:
  Cox-Snell residuals are defined in a way that they always follow an
  exponential distribution. Since this is an unususal law for residuals,
  it is convenient to transform them such that they then obey a standard
  normal distribution. See the vignette for more detail.
}
\value{
  Vector of residual values. If conditional quantiles are requested,
  it carries an attribute \code{"condquant"}, which is 
  a data.frame with the variables
  \item{median}{median of the conditional distributions}
  \item{lowq}{lower quartile}
  \item{uppq}{upper quartile}
  \item{random}{random number, drawn according to the conditional
    distribution}
  \item{prob}{probability of the condition being true}
  \item{index}{index of the observation in the sequence of the result
    (residuals)}
  \item{fit}{linear predictor value}
  \item{y}{observed response value}
}
\references{See http://stat.ethz.ch/~stahel/regression}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{condquant}}, \code{\link{plot.regr}}}
\examples{
require(MASS)
data(housing, package="MASS")
rr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
t.res <- residuals.regrpolr(rr)
head   (t.res)
summary(t.res)
}
\keyword{regression}
