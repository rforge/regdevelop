\name{logst}
\alias{logst}
\title{Started Logarithmic Transformation}
\description{
Transforms the data by a log10 transformation, modifying small and zero
observations such that the transformation yields finite values.
}
\usage{
logst(data, calib=data, threshold=NULL, mult = 1)
}
\arguments{
  \item{data}{a vector or matrix of data, which is to be transformed}
  \item{calib}{a vector or matrix of data used to calibrate the
    transformation(s),
    i.e., to determine the constant \code{c} needed}
  \item{threshold}{constant c that determines the transformation}
  \item{mult}{a tuning constant affecting the transformation of small
    values, see Details}
}
\details{
  Small values are determined by the threshold c. If not given by the
  argument \code{threshold}, then it is determined by the quartiles
  \eqn{q_1}{q1} and \eqn{q_3}{q3} of the non-zero data as those
  smaller than \eqn{c=q_1 / (q_3/q_1)^{mult}}{q1 / (q3/q1)^{mult}}.
  The rationale is that for lognormal data, this constant identifies
  2 percent of the data as small.
  Beyond this limit, the transformation continues linear with the
  derivative of the log curve at this point. See code for the formula.

  The function chooses log10 rather than natural logs because they can
  be backtransformed relatively easily in the mind.
}
\value{
  the transformed data. The value c needed for the transformation is
  returned as \code{attr(.,"threshold")}. 
}
\author{Werner A. Stahel, ETH Zurich}
\note{
  The names of the function alludes to Tudey's idea of "started logs".
}
\examples{
dd <- c(seq(0,1,0.1),5*10^rnorm(100,0,0.2))
dd <- sort(dd)
r.dl <- logst(dd)
plot(dd, r.dl, type="l")
abline(v=attr(r.dl,"threshold"),lty=2)
}
\keyword{manip}
