\name{residuals.survreg}
\Rdversion{1.1}
\alias{residuals.survreg}
\alias{residuals.coxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residuals of a regression for survival data
}
\description{
  Calculates residuals for uncensored observations and
  quartiles and random numbers according to the
  conditional distribution of residuals for censored observations.
  See Details for an explanation.
}
\usage{
\method{residuals}{survreg}(object, type = "response", na.action=object, ...)
\method{residuals}{coxph}(object, type = "CoxSnellMod", na.action=object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object inheriting from classes Survreg or coxph
}
  \item{type}{type of residuals required
}
  \item{na.action}{a list with component \code{na.action}, which
    determines the na.action for the residual component of \code{object}
}
  \item{\dots}{unused
}
}
\details{
For censored observations, the conditional distribution of the survival
time, given the censoring, is used to calculate quantiles, see
\code{\link{condquant}}
}
\value{
  matrix with class \code{condquant}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{See http://stat.ethz.ch/~stahel/regression
}
\author{Werner A. Stahel, ETH Zurich
}
\note{BUG: \code{residuals.coxph} is not yet implemented in this way.
}

\seealso{\code{\link{condquant}}, \code{\link{plot.regr}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(survival)
t.r <- regr(Surv(time,status)~x, data=aml, family="weibull")
residuals(t.r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
