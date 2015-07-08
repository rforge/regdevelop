\name{Tobit}
\alias{Tobit}

\title{
Prepare a Response for a Tobit Model
}
\description{
Returns a \code{Surv} object that allows for setting up a Tobit
regression model by calling \code{survreg} 
}
\usage{
Tobit(data, limit = 0, transform = NULL, log = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    the variable to be used as the response in the Tobit regression
}
  \item{limit}{
    Lower limit which censors the observations.
    If \code{log} is \code{TRUE}, then the default is the minimal value
    of \code{\link{logst}}(data), and if \code{limit>0}, it refers to the
    untransformed data.
}
  \item{transform}{
    if data should be transformed, specify the function to be used.
  }
  \item{log}{
    logical. If \code{TRUE}, \code{data} will be log transformed
    by calling \code{\link{logst}}.
  }
  \item{...}{
    any additional arguments to the \code{transform} function
  }
}
\details{
Tobit regression is a special case of regression with left censored
response data. The function \code{survreg} is suitable for fitting.
In \code{regr}, this is done automatically.
}
\value{
A \code{Surv} object.
}
\author{
Werner A. Stahel
}

\seealso{
  \code{Surv}
}
\examples{
data("tobin", package="survival")

Tobit(tobin$durable) 
(t.r <- regr(Tobit(durable) ~ age + quant, data=tobin))
plot(t.r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}

