\name{d.blast}
\alias{d.blast}
\docType{data}
\title{Blasting for a tunnel}
\description{Blasting causes tremor in buildings, which can lead to
  damages. This dataset shows the relation between tremor and
  distance and charge of blasting.
}
\usage{data("d.blast")}
\format{
  A data frame with 388 observations on the following 7 variables.
  \describe{
    \item{\code{date}}{date in Date format}
    \item{\code{location}}{Code for location of the building,
      \code{loc1} to \code{loc8}}
    \item{\code{device}}{Number of measuring device, 1 to 4}
    \item{\code{distance}}{Distance between blasting and location
      of measurement}
    \item{\code{charge}}{Charge of blast}
    \item{\code{tremor}}{Tremor energy (target variable)}
  }
}
\details{The charge of the blasting should be controled in order to
  avoid tremors that exceed a threshold.
  This dataset can be used to establish the suitable rule:
  For a given \code{distance}, how large can \code{charge} be in order
  to avoid exceedance of the threshold?
}
\source{Basler and Hoffmann AG, Zurich
}
%%- \references{
%%- %%  ~~ possibly secondary sources and usages ~~
%%- }
\examples{
data(d.blast)

summary(lm(log10(tremor)~location+log10(distance)+log10(charge),
           data=d.blast))
}
\keyword{datasets}
