\name{d.hail}
\alias{d.hail}
\docType{data}
\title{
d.hail: hail prevention data
}
\description{
data from an experiment to test the suitability of a hail
prevention method
}
\usage{data("d.hail")}
\format{
  A data frame with 253 observations on the following 23 variables.
  \describe{
    \item{\code{TGNR}}{a numeric vector}
    \item{\code{SEED}}{a numeric vector}
    \item{\code{HV}}{a numeric vector}
    \item{\code{TV}}{a numeric vector}
    \item{\code{G0}}{a numeric vector}
    \item{\code{G2}}{a numeric vector}
    \item{\code{DT46}}{a numeric vector}
    \item{\code{E0}}{a numeric vector}
    \item{\code{A46}}{a numeric vector}
    \item{\code{ZM0}}{a numeric vector}
    \item{\code{H0}}{a numeric vector}
    \item{\code{TB}}{a numeric vector}
    \item{\code{HB}}{a numeric vector}
    \item{\code{TAU0}}{a numeric vector}
    \item{\code{T0}}{a numeric vector}
    \item{\code{EGR}}{a numeric vector}
    \item{\code{EGR0}}{a numeric vector}
    \item{\code{E56}}{a numeric vector}
    \item{\code{E61}}{a numeric vector}
    \item{\code{RGR}}{a numeric vector}
    \item{\code{RGR0}}{a numeric vector}
    \item{\code{R56}}{a numeric vector}
    \item{\code{R61}}{a numeric vector}
  }
}
\details{
This dataset was used to assess the suitability of a method for hail
suppression based on seeding clouds with silver iodide.
}
\source{
B. Federer and A. Waldvogel and W. Schmid and 
                  H. H. Schiesser and F. Hampel and M. Schweingruber and 
                  W. Stahel and J. Bader and J.-F. Mezeix and N. Doras and 
                  G. D'Aubigny and G. DerMegreditchian and D. Vento (1986).
Main results of Grossversuch IV.
J. Climate Appl. Meteorology 25, 917-957
}
\examples{
data(d.hail)
## maybe str(d.hail) ; plot(d.hail) ...
}
\keyword{datasets}
