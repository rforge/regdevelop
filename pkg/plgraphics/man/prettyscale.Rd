\name{prettyscale}
\alias{prettyscale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pretty Tickmark Locations for Transformed Scales
}
\description{
  Compute about \code{n} 'round' values that are about equally spaced
  in a transformed (plotting) scale and cover the range of the values
  in \code{x}.
}
\usage{
prettyscale(x, plscale = "log10", inverse = NULL, n = NULL,
    logscale = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric vector of data
  }
  \item{plscale}{name of the transformation defining the plotting scale
  }
  \item{inverse}{back (or inverse) back transformation
  }
  \item{n}{approximate number of tickmark locations.
    If of length \code{>=2}, \code{n[2]} can be varied to obtain more
    adequate locations. See Details    
  }
  \item{logscale}{if \code{NULL}, R's function \code{axTicks}
    will be called if the plscale is a log function.
  }
}
\details{
!!!
}
\value{Numeric vector of tick mark locations in transformed scale,
  with an attribute \code{ticklabels} containing the appropriate
  tick marks and labels (in original scale)
}
\author{W. A. Stahel
}

\seealso{
\code{axTicks}
}
\examples{
  prettyscale(10^rnorm(10))
}
\keyword{dplot}
\keyword{attribute}
