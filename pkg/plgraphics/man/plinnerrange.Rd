\name{plinnerrange}
\alias{plinnerrange}

\title{Inner Plotting Limits
}
\description{Calculates inner limits for plotting, based on a robust
  estimate of the range.
}
\usage{
plinnerrange(innerrange, data, factor = 4, FUNC = robrange)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{innerrange}{logical: Should range be calculated?
    If \code{FALSE}, the result will contain only the values
    \code{FALSE}.
    If it is a list or matrix of the approriate size, it will be
    returned as is.
  }
  \item{data}{vector or data.frame for which the range(s) will be calculated
  }
  \item{factor}{expansion of the calculated robust range to yield the
    plotting range
  }
  \item{FUNC}{function used to calculate the robust range.
    The \code{factor} will be handed over to \code{FNC} as the argument
    \code{fac}.
  }
}
\value{Matrix of 2 rows giving the ranges to be used as
  inner plotting ranges for the variables.
  If \code{innerrange} is such a matrix or data.frame, it will be returned
  as is.
}
\author{Werner A. Stahel
}

\seealso{
\code{\link{robrange}}, \code{\link{plcoord}}
}
\examples{
data(d.blast)
dd <- d.blast[,c("charge","distance","tremor")]
( t.ipl <- plinnerrange(TRUE, dd) )
plot(dd[,"tremor"], plcoord(dd[,"tremor"], t.ipl[,"tremor"]))
abline(h=t.ipl[,"tremor"])
}
\keyword{utilities}

