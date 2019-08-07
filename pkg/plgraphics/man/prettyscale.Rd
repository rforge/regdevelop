\name{prettyscale}
\alias{prettyscale}

\title{Pretty Tickmark Locations for Transformed Scales
}
\description{
  Compute about \code{n} 'round' values that are about equally spaced
  in a transformed (plotting) scale and cover the range of the values
  in \code{x}.
}
\usage{
prettyscale(x, plscale = "log10", inverse = NULL, range = NULL,
    range.transformed = NULL, n = NULL, logscale = NULL)
}

\arguments{
  \item{x}{numeric vector of data (original scale)
  }
  \item{plscale}{name of the transformation defining the plotting scale
  }
  \item{inverse}{back (or inverse) back transformation
  }
  \item{range, range.transformed}{possible range of original and
    transformed values, respectively. Usually not needed, cf. Details
  }
  \item{n}{approximate number of tickmark locations.
    If of length \code{>=2}, \code{n[2]} can be varied to obtain more
    adequate locations. See Details.
  }
  \item{logscale}{if \code{NULL}, R's function \code{axTicks}
    will be called if the plscale is a log function.
  }
}
\details{\code{prettyscale} generates \code{n+2} "anchor" values in the
  transformed scale which cover the range of the transformed \code{x}
  values and are equidistant within the range.
  It then back-transforms these anchor values. For each one of them,
  say \code{c}, 
  it seeks a pretty value near to it by the following construction:
  it calls the R function \code{pretty} on the range given by the
  back-transformed neighboring anchor values, asking for \code{n[2]}
  pretty values. From these, it chooses the one for which the
  transformed value is closest to the transformed \code{c}.

  Therefore, if \code{n[2]} is large, the pretty values may be less
  pretty, whereas small \code{n[2]} may lead to equal pretty values
  for neighboring anchors and thus to too few resulting pretty values.
  The default value for \code{n[2]} is 3.

  The ranges are needed to get the limits as pretty values when
  appropriate (and to avoid warning messages).
  They are generated in the function for the commonly used plscales
  and may be given as attributes of the \code{plscale} function,
  see Examples.
}
\value{Numeric vector of tick mark locations in transformed scale,
  with an attribute \code{ticklabels} containing the appropriate
  tick marks and labels (in original scale)
}
\note{The function does not always lead to consistent results.
  Increasing \code{n} sometimes leads to fewer resulting values.
}
\author{W. A. Stahel
}

\seealso{
\code{axTicks}
}
\examples{
  prettyscale(10^rnorm(10))
  prettyscale(c(0.5, 2, 10, 90), plscale="sqrt")
  prettyscale(c(50,90,95,99), plscale="asinp", n=10)
  ## asinp has the useful attributes:
  asinp
}
\keyword{dplot}
\keyword{attribute}
