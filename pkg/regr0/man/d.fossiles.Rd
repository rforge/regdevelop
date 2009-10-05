\name{d.fossiles}
\alias{d.fossiles}
\docType{data}
\title{Shapes of Shells and Environmental Variables}
\description{
  This dataset relates shapes of Gephyrocapsa shells with environmental
  variables
}
\usage{data(d.fossiles)}
\format{
  The format is:
 chr "d.fossiles"
}
\details{
  NA
}
\source{
  Joerg Bollmann, Jorijntje Henderiks and Bernhard Brabec,
  Geological Institute, ETH Zurich,
}
\references{
  Joerg Bollmann, Jorijntje Henderiks and Bernhard Brabec,
  Geological Institute, ETH Zurich,\cr
  Marine Micropaleontology 29, 319-350 (1997)
}
\examples{
data(d.fossiles)
showd(d.fossiles)
regr(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+lChlorophyll+region+N,
                data=d.fossiles)
}
\keyword{datasets}
