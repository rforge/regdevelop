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
condquant(x, dist = "normal", mu = 0, sig = 1, randomrange = 0.9)
}
\arguments{
  \item{x}{matrix with 2 columns or vector of length 2 giving the
    limits for the conditional distribution}
  \item{dist}{(unconditional) distribution. Currently, only
    \code{normal} (or \code{gaussian}),
    \code{logistic} and \code{revgumbel} (reverse-Gumbel,
    distribution of the logarithm of a Weibull variable) are implemented.
  }
  \item{mu, sig}{locarion and scale parameter of the distribution
  }
  \item{randomrange}{random numbers from the conditional distribution
    are drawn for the inner \code{100*randomrange} percent of the
    suitable p-range. This avoids random extreme outliers and points
    close to the limit of the intervals on which they are conditioned.}
}
\value{Matrix consisting of a row for each row of \code{x} for which
  x[,1] differs from x[,2] and the following columns:
  \item{median}{Median}
  \item{lowq, uppq}{lower and upper quartiles}
  \item{random}{random number according to the conditional distribution
    (one for each row)}
  \item{prob}{probability of the condition being true}
  \item{index}{(row) index of the corresponding entry in the input 'x'}
}
\author{Werner A. Stahel, Seminar for Statistics, ETH Zurich}
%%- \note{The function is needed for \code{\link{residuals.polr}} and
%%-   \code{\link{residuals.survreg}}
%%- }

\examples{
condquant(cbind(seq(-2,1),c(0,1,Inf,1)))
}
\keyword{distribution}
