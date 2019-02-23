\name{residuals.regrpolr}
\alias{residuals.polr}
\alias{residuals.regrpolr}
\alias{residuals.regrsurvreg}
\alias{residuals.regrcoxph}
% \alias{residuals.polr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residuals of a Binary, Ordered, or Censored Regression}
\description{Methods of \code{residuals} for classes
  \code{polr, survreg} and \code{coxph},
  calculating quartiles and random numbers according to the
  conditional distribution of residuals for the latent variable of a
  binary or ordinal regression or a regression with censored response,
  given the observed response value.
  See Details for an explanation.
}
\usage{
\method{residuals}{polr}(object, type="condquant", ...)
\method{residuals}{regrpolr}(object, type="condquant", ...)
%% \method{residuals}{polr}(object, type=c("condquant","deviance"), ...)
\method{residuals}{regrsurvreg}(object, type="condquant", ...)
\method{residuals}{regrcoxph}(object, type="CoxSnellMod", ...)
}
\arguments{
  \item{object}{the result of \code{polr}, of
    \code{glm(,family=binomial)} with binary data
    for the \code{regrpolr} method,
    or of \code{survreg} or \code{coxph} for the respective methods.}
  \item{type}{type of residuals:
    \code{"condquant"} requires conditional quantiles (and more) of the
    residuals of the model, see Details.\cr
%%-  for the latent variable of binary or ordered regression,
%%-  given the response.
%%-  For regression with censored response, conditional quantiles
%%-  are given for the observations with censored values.
%%-  See \code{\link{condquant}} for more information.  
%% or \code{"deviance"} for deviance residuals\cr
  For \code{residuals.regrsurvreg}, type \code{CoxSnellMod} yields
  a modified version of Cox-Snell residuals, also including
  a \code{condquant} attribute, see Details.\cr
  Other \code{type}s are available according to the standard
  methods of the \code{residuals} function.
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
\note{\code{residuals.polr} and \code{residuals.regrpolr} are identical
  for the time being. Only \code{type="condquant"} is available now.
}
\value{
  Vector of residual values. If conditional quantiles are requested,
  the residuals for censored observations are replaced by conditional
  medians, and an attribute \code{"condquant"} is attached, which is 
  a data.frame with the variables
  \item{median}{median of the conditional distributions}
  \item{lowq}{lower quartile}
  \item{uppq}{upper quartile}
  \item{random}{random number, drawn according to the conditional
    distribution}
  \item{prob}{probability of the condition being true}
  \item{limlow, limup}{lower and upper limits of the intervals}
  \item{index}{index of the observation in the sequence of the result
    (residuals)}
  \item{fit}{linear predictor value}
  \item{y}{observed response value}
}
\references{See http://stat.ethz.ch/~stahel/regression}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{condquant}}, \code{\link{plregr}}}
\examples{
require(MASS)
data(housing, package="MASS")
rr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
t.res <- residuals.regrpolr(rr)
head   (t.res)
summary(t.res)
}
\keyword{regression}
