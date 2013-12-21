\name{condquant}
\alias{condquant}
\title{Quantiles of a Conditional Distribution}
\description{
  Calculates quantiles of a conditional distribution, as well
  as corresponding random numbers.  The condtion is simply to restrict
  the distribution (given by \code{dist}) to a range (given by
  \code{x})
}
\usage{
condquant(x, dist = "normal", sig = 1)
}
\arguments{
  \item{x}{matrix with 2 columns or vector of length 2 giving the
    limits for the conditional distribution}
  \item{dist}{(unconditional) distribution. Currently, only normal and
    logistic are implemented.
}
  \item{sig}{standard deviation for normal distribution
}
}
\value{Matrix consisting of a row for each row of \code{x} and the
  following columns:
  \item{median}{Median}
  \item{lowq, uppq}{lower and upper quartiles}
  \item{random}{random number according to the conditional distribution
  (one for each row)}
  \item{prob}{probability of the condition being true}
}
\author{Werner A. Stahel, Seminar for Statistics, ETH Zurich}
\note{The function is needed for \code{\link{residuals.polr}} and
  \code{\link{residuals.survreg}}
}

\examples{
condquant(cbind(seq(-2,0),c(0,1,Inf)))
}
\keyword{distribution}
