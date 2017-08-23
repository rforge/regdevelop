\name{plcoord}
\alias{plcoord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determines Values for Plotting with Limited "Inner" Plot Range }
\description{
  For plots with an "inner plot range" (see Details) this function
  converts the data values to the coordinates in the plot 
}
\usage{
plcoord(x, range = NULL, limfac = 3, limext = 0.1, plext = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ data to be represented }
  \item{range}{ vector of 2 elements giving the inner plot range. Data
    beyond the given interval will be non-linearly transformed to fit
    within the (outer) plot margins. Defaults to
    \code{\link{robrange}(x, fac=fac)}.
  }
  \item{limfac}{ factor used to determine the default of \code{range} }
  \item{limext}{ factor for extending the \code{range} to determine the
    outer plot range }
  \item{plext}{ vector of 1 or 2 elements setting the extension factor
    for the plotting range}
}
\details{
  When plotting data that contain outliers, the non-outlying data is
  represented poorly. Rather than simply clipping outliers, one can
  split the plotting area into an inner region, where the (non-outlying)
  data is plotted as usual, and a plot area margin, in which outliers
  are represented on a highly non-linear scale that allows to display
  them all.

  This function converts the data to the coordinates used in the
  graphical display, and also returns the inner and outer ranges for
  plotting.
}
\value{
  vector of coordinates used for plotting, that is, unchanged \code{x}
  values for   those within the \code{range} and transformed values
  for those outside.

  Attributes:
  \item{attr(,"plrange")}{the range to be used when plotting}
  \item{attr(,"range")}{the "inner" plot range, either the argument
  \code{range} or the values determined by default.}
  \item{attr(,"nmod")}{the number of modified observations}
}
% \references{ ~put references to the literature/web site here ~ }
\author{Werner Stahel}
\seealso{ \code{\link{robrange}} }
\examples{
  x <- c(rnorm(20),rnorm(3,5,20))
  xmod <- plcoord(x)
  plot(x,xmod)
  plot(xmod)
  abline(h=attr(xmod,"range"),lty=3, lwd=2)
}
\keyword{manip}
\keyword{aplot}
