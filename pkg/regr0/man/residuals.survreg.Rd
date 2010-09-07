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
\method{residuals}{survreg}(object, type = "response", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object inheriting from classes Survreg or coxph
}
  \item{type}{
}
  \item{\dots}{unused
}
}
\details{
For censored observations, the conditional distribution of the survival
time, given the censoring, is used to calculate quantiles ...\cr
}
\value{
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
