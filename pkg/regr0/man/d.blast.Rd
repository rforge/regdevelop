\name{d.blast}
\alias{d.blast}
\docType{data}
\title{tremor caused by blasting for excavation of a tunnel, with
  explanatory information}
\description{
  When excavating a tunnel, blasting may cause damage in houses 
  in the neighborhood. It is important to keep the tremor caused by the
  blasting under control. To this end, tremor is measured.
  There were 4 devices to measure tremor. They were moved to different
  locations while the excavation advanced. 
}
\usage{data(d.blast)}
\format{
  A data frame with 388 observations on the following 7 variables.
  \describe{
    \item{\code{no}}{number of blasting}
    \item{\code{datetime}}{date and time of the blasting}
    \item{\code{device}}{measurement device}
    \item{\code{charge}}{charge (loading) in kg}
    \item{\code{distance}}{between place of blasting and location}
    \item{\code{tremor}}{measured in mm/s}
    \item{\code{location}}{house in which tremor was measured}
  }
}
\source{
  Basler \& Hofmann, Zurich
}
%- \references{
%-   ``` possibly secondary sources and usages ```
%- }
\examples{
data(d.blast)
doc(d.blast)
coplot(log10(tremor)~distance|charge*location, data=d.blast)
}
\keyword{datasets}
