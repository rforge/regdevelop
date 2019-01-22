\name{d.fossileShapes}
\alias{d.fossileShapes}
\alias{d.fossileSamples}
\docType{data}
\title{Coccolith Abundance and Environmental Variables
}
\description{The abundance of cocolith shells can be used to infer
  environmental conditions in epochs corresponding to earlier epochs.
  This data set contains the core location, the relative abundance of
  Gephyrocapsa morphotypes and the sea surface temperatures from all
  deep see cores used in this study.
}
\usage{data("d.fossileShapes")
data("d.fossileSamples")}
\format{
  \code{d.fossilShapes}:
    A data frame with 5864 observations on the following 15
    variables:\cr
    Identification and location of the sample:
  \describe{
    \item{\code{Sample}}{Identification number of the sample}
    \item{\code{Sname}}{Identification code}
    \item{\code{Magnification}}{(technical)}
  }
  Shape features and recommended transformations:
  \describe{
    \item{\code{Angle}}{bridge angle}
    \item{\code{Length, Width}}{lengtha and width of the shell}
    \item{\code{CLength, CWidth}}{length and width of the 'central area'}
    \item{\code{Cratio}}{ratio between width and length of the central area}
    \item{\code{sAngle}}{sqrt of \code{Angle}}
    \item{\code{lLength}}{\code{log10(Length)}}
    \item{\code{rWidth, rCLength, rCWidth}}{relative measures,
      percentage of \code{Length}}
    \item{\code{Cratio}}{\code{CWidth/Clength}}
    \item{\code{ShapeClass}}{shape class as defined in the cited paper,
      classes ar
      \code{CM} < \code{CC} < \code{CT} < \code{CO} < \code{CE} < \code{CL}}
  }
  
  \code{d.fossilSamples}:
    A data frame with 108 observations on the following 32
    variables:\cr
    Identification and location:
    \describe{
      \item{\code{Sample}}{Identification number of the sample (as above)}
      \item{\code{Sname}}{Identification code}
      \item{\code{Latitude, Longitude}}{Coordinates of the location}
      \item{\code{Region}}{Ocean: \code{Pacific, Atlantic, Indian.Ocean}}
      \item{\code{SDepth}}{sample depth below soil surface [cm]}
      \item{\code{WDepth}}{Water depth [m]}
      \item{\code{N}}{number of specimen measured}
    }
    Shape features as above, averaged.
    (This is the reason for introducing transformed variables above:
    The transformed values are averaged.)\cr
    \describe{
      \item{\code{CM, CC, CT, CO, CE, CL}}{percentages of shape classes
	in the sample}
    }
    Environment:
    \describe{
      \item{\code{SST}}{Sea Surface Temperature, mean, [deg C]}
      \item{\code{SST.Spring}, \code{SST.Summer}, \code{SST.Fall},
	\code{SST.Winter}}{... in each season}
      \item{\code{Chlorophyll, lChlorophyll}}{Chlorophyll content
	[microgram/L] and \code{log10} of it}
      \item{\code{Salinity}}{Salinity of the sea water}
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
\references{
  J\"org Bollmann (1997).
  Morphology and biogeography of Gephyrocapsa coccoliths in Holocene
  sediments.
  Marine Micropaleontology, 29, 319-350
}
\examples{
data(d.fossileShapes)
names(d.fossileShapes)

data(d.fossileSamples)
plyx(sqrt(Angle) ~ SST, data=d.fossileSamples)

}
\keyword{datasets}
