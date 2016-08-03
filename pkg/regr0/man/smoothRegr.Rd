\name{smoothRegr}
\alias{smoothRegr}

\title{Smoothing function used as a default in regr
}
\description{This function wraps the \code{loess} smoothing function in
  order to meet the argument conventions used in the \code{regr0} package.
}
\usage{
smoothRegr(x, y, weights, par = 5 * length(x)^log10(1/2),
  iterations = 50, minobs=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of x values
}
  \item{y}{vector of y values to be smoothed
}
  \item{weights}{vector of weigths used for fitting the smooth 
%%     ~~Describe \code{weights} here~~
}
\item{par}{value for the \code{span} argument of \code{loess}.
}
\item{iterations}{number of iterations for the \code{loess} algorithm.
  If \code{==1}, the non-robust, least squares version is applied.
}
\item{minobs}{minimal number of observations. If less valid observations
  are provided, the result is \code{NULL}.
}
}
\value{vector of smoothed values. If \code{loess} fails, NAs will be
  returned without issuing a warning.}
\author{Werner A. Stahel, ETH Zurich}

\seealso{\code{\link{loess}}, \code{\link{smoothM}}
}
\examples{
t.x <- (1:50)^1.5
t.y <- log10(t.x) + rnorm(length(t.x),0,0.3)
t.y[40] <- 5
r.sm <- smoothRegr(t.x, t.y, par=0.5)
r.sm1 <- smoothRegr(t.x, t.y, iterations=1, par=0.5)

plot(t.x,t.y)
lines(t.x,r.sm, col=2)
lines(t.x,r.sm1, col=3)
}
\keyword{regression}

