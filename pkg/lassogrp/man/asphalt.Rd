\name{asphalt}
\alias{asphalt}
\docType{data}
\title{Asphalt Rutting Data}
\description{
  Measurement series on the rutting of asphalt in function of some of
  its material properties.
}
\usage{data(asphalt)}
\format{
  A data frame with 31 observations on the following 7 variables
  \describe{
    \item{\code{VISC}}{numeric; the viscosity of the asphalt.}
    \item{\code{ASPH}}{numeric; the percentage of asphalt in the surface
      course.}
    \item{\code{BASE}}{numeric; the percentage of asphalt in the base
      course.}
    \item{\code{RUN}}{a vector of 0:1 values, indicator of the two runs
      of measurements.}
    \item{\code{FINES}}{numeric; 10* percentage of fines in surface course.}
    \item{\code{VOIDS}}{numeric; percentage of voids in surface course.}

    \item{\code{RUT}}{numeric, the \dQuote{rate of rutting}, i.e. the
      change of rut depth in inches per million wheel passes.}
  }
}
\details{
  Typical model for regression teaching, including variable transformations.
}
\source{
  Gorman, J.W. and Toman, R.J. (1966).
  Selection of Variables for Fitting Equations to Data;
  \emph{Technometrics} \bold{8}, 27--51.

  Daniel, C. and Wood, F.S. (1980).
  \emph{Fitting Equations to Data}, 2nd ed., Wiley, N.Y.
}
\references{
  Detailed applications in Daniel and Wood, 6.5, p. 95~ff.
}
\examples{
data(asphalt)
pairs(asphalt, gap = .1)
pairs(asphalt[, -4], gap = .1,
      main = "asphalt -- color coded  'RUN'",
      col = 1+ asphalt[, 4])
}
\keyword{datasets}
