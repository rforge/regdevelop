\name{xdistResdiff}
\alias{xdistResdiff}
\alias{xdistResscale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residual Differences for Near Replicates: Tabulate and Test
}
\description{
A test for the completeness of a linear regression model can be
performed based on comparing the differences of residuals for
pairs of observations that are close to each other to the
estimated standard deviation of the 
}
\usage{
xdistResdiff(object, perc = c(3, 10, 80), trim = 0.1,
  nmax = 100, nsim = 100, out = "aggregate")
xdistResscale(x, perc = c(3, 10, 90), trim = 1/6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object containing the result of fitting a linear
    model by \code{regr}}
  \item{x}{an object produced by \code{xdistResdiff}}
  \item{perc}{Percentage points to define distance classes}
  \item{trim}{Trimming proportion for calculating means of
    absolute residual differences}
  \item{nmax}{maximal number of observations to form pairs}
  \item{nsim}{number of simulations for determining the null
    distribution of test statistics}
  \item{out}{determines the value of \code{xdistResdiff}:
    if \code{=="aggregate"} (the default), the value will be
    produced by calling \code{xdistResscale}, otherwise,
    all x distances and respective residual differences will be
    returned.}
}
\details{
  See package vignette.
}
\value{
  For \code{xdistResdiff} with \code{out="aggregate"} and
  \code{xdistResscale}, a matrix is returned with a row for
  each class of x distances and the columns
  \item{xdist}{mean x distance}
  \item{rdiff.mean}{absolute differences of residuals for pairs of
    observations in the distance class, averaged over the class}
  \item{rdiff.simmean}{mean of (trimmed) means for simulated data}
  \item{rdiff.se}{standard error of (trimmed) means as obtained from
    simulation}
  The matrix carries along the following attributes:
  \item{perc}{given argument \code{perc}}
  \item{xd.classlim}{the actual class limits corresponding to
    \code{perc}}
  \item{trim}{given argument \code{trim}}
  \item{rdiff.grandmean}{overall mean of absolute residual differences}
  \item{p-values}{p values for the classes as obtained from simulation,
    and p-value for the sum of squares statistic}
  \item{class}{The value has S3 class \code{xdistResscale} and
    \code{matrix}}.

  If \code{xdistResdiff} with \code{out} different from
  \code{"aggregate"}, then a data.frame is returned containing a row for
  each pair of observations and the columns
  \item{id1, id2}{the labels of the two observations}
  \item{xdist}{the x distance between the two observations}
  \item{resdiff}{the difference of residuals for the two observations}
  The value has S3 class \code{xdistResdiff} and \code{data.frame}.
}
\references{
  See package vignette.
}
\author{Werner A. Stahel, ETH Zurich}

\examples{
data(d.blast)
rr <- regr(tremor~distance+charge, data=d.blast)
## an inadequate model!
xdistResdiff(rr)
xdrd <- xdistResdiff(rr, out="all")
showd(xdrd[,1:10])
xdistResscale(xdrd)
## same as first call of  xdiffResdiff
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{test}% __ONLY ONE__ keyword per line
