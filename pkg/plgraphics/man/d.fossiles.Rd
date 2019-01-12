\name{d.fossiles}
\alias{d.fossiles}
\docType{data}
\title{Coccolith Abundance and Environmental Variables
}
\description{The abundance of cocolith shells can be used to infer
  environmental conditions in epochs corresponding to earlier epochs.
  This data set contains the core location, the relative abundance of
Gephyrocapsa morphotypes and the sea surface temperatures from all
deep see cores used in this study.
}
\usage{data("d.fossiles")}
\format{
  A data frame with 108 observations on the following 39 variables.

  Identification and location of the core:
  \describe{
    \item{\code{Sample}}{Identification number}
    \item{\code{Sname}}{Identification code}
    \item{\code{SDepth}}{a numeric vector}
    \item{\code{Latitude ,Longitude}}{Coordinates of the location}
    \item{\code{WDepth}}{Water depth}
  }
  Environment:
  \describe{
    \item{\code{SST.Mean}}{a numeric vector}
    \item{\code{SST.Spring}}{a numeric vector}
    \item{\code{SST.Summer}}{a numeric vector}
    \item{\code{SST.Fall}}{a numeric vector}
    \item{\code{SST.Winter}}{a numeric vector}
    \item{\code{Chlorophyll}}{a numeric vector}
    \item{\code{SST.Max}}{a numeric vector}
    \item{\code{SST.Min}}{a numeric vector}
    \item{\code{Salinity}}{a numeric vector}
    \item{\code{N}}{number of specimen measured}
  }
  Abundances of the classes \code{T, E, C, L, O, M}:
  \describe{
    \item{\code{GT}}{a numeric vector}
    \item{\code{GE}}{a numeric vector}
    \item{\code{GC}}{a numeric vector}
    \item{\code{GL}}{a numeric vector}
    \item{\code{GO}}{a numeric vector}
    \item{\code{GM}}{a numeric vector}
  }
  \describe{
    \item{\code{SST.Mean.est}}{a numeric vector}
    \item{\code{SST.Spring.est}}{a numeric vector}
    \item{\code{SST.Fall.est}}{a numeric vector}
    \item{\code{SST.Winter.est}}{a numeric vector}
    \item{\code{SST.Summer.est}}{a numeric vector}
    \item{\code{SST.Max.est}}{a numeric vector}
    \item{\code{SST.Min.est}}{a numeric vector}
  }
  Shape variables:
  \describe{
    \item{\code{Angle}}{a numeric vector}
    \item{\code{Length}}{a numeric vector}
    \item{\code{rWidth}}{a numeric vector}
    \item{\code{rCLength}}{a numeric vector}
    \item{\code{Cratio}}{a numeric vector}
    \item{\code{sAngle}}{a numeric vector}
    \item{\code{lLength}}{a numeric vector}
    \item{\code{lChlorophyll}}{a numeric vector}
    \item{\code{region}}{Region, levels \code{Atlantik, Ind.Ozean, Pazifik}}
    \item{\code{region.num}}{a numeric vector}
  }
}
\details{The paradigm of research associated with this dataset is the
  following:
  Datasets of this kind are used to establish the relationship
  between the shell shapes of cocoliths (species Gephyrocapsa)
  from the most recent sediment
  layer with actual environmental conditions. This relationship is then
  used to infer environmental conditions of earlier epochs from
  the shell shapes from the corresponding layers.

  The analysis presented in the paper cited below consisted of first
  introducing classes of shells based on the shapes and then use the
  relative abundance of the classes to predict the environmental
  conditions. 
}
\source{
  J\"org Bollmann, Jorijntje Henderiks and Bernhard Brabec (2002).
  Global calibration of Gephyrocapsa coccolith abundance in Holocene
  sediments for paleotemperature assessment.
  Paleoceanography, 17(3), 1035
}
%%- \references{
%%- }
\examples{
data(d.fossiles)
names(d.fossiles)

}
\keyword{datasets}
