\name{plticks}
\alias{plticks}

\title{Ticks for plotting
}
\description{Find ticks locations and labels
}
\usage{
plticks(range, plscale = NULL, transformed = FALSE, nouter = 0,
  tickintervals = NULL, ploptions = NULL)
}

\arguments{
  \item{range}{range of values that the ticks should cover
  }
  \item{plscale}{function defining the scale of the axis.
    Either the name of the function or a function, see Details.
  }
  \item{transformed}{logical: Is \code{range} scaled according to
    \code{plscale} rather than in original scale?
  }
  \item{nouter}{
  }
  \item{tickintervals}{approximate number of tick intervals desired.
    Default is taken from \code{ploptions('tickintervals')}.
  }
  \item{ploptions}{pl options
  }
}
\details{\code{plticks} calls \code{\link{pretty}} for getting
  tick locations if \code{plscale} is not specified and
  \code{\link{prettyscale}} if it is.
  It generates another set for locations of tick labels if
  \code{tickintervals} has 2 elements, such that not all ticks
  are labelled. %% see \code{\link{prettyscale}}

  The scaling function \code{plscale} can be given by its name
  if that name is one of
  \code{log, log10, logst, sqrt, asinp, logit, qnorm}.
  Otherwise, it must be a function with an attribute
  \code{inverse} that defines the inverse function.
  It should also have an attribute \code{range} and an
  attribute \code{range.transformed} if the possible
  range for its argument or its values are restricted,
  like \code{\link{asinp}} that is defined for values between 0 and 100
  and has values in the interval from 0 to 1.
}
\value{
  A list with components
  \item{ticksat}{locations of ticks}
  \item{ticklabelsat}{locations of tick labels}
  \item{ticklabels}{tick labels, if \code{plscale} is given}
}
\author{Werner A. Stahel
}

\seealso{
\code{\link{pretty}}, \code{\link{prettyscale}}, \code{\link{plaxis}}
}
\examples{
plticks(c(23,87))
plticks(c(23,91), plscale="asinp", transformed=FALSE,
  tickintervals=c(10,2))

asinp ## shows the attributes 'inverse', 'range' and 'range.transformed'
}
\keyword{dplot}
\keyword{attribute}
