\name{plsubset}
\alias{plsubset}

\title{Subsetting a Data.Frame with pl Attributes
}
\description{Select rows of data.frames keeping the variable attributes
  that drive pl graphics
}
\usage{
plsubset(x, subset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data.frame from which the subset is to be generated
  }
  \item{subset}{logical vector or vector of indices of rows
  }
}
\details{\code{plsubset} maintains the 'pl' attributes of the variables
  of the data.frame (if there are), such as 'col', 'lty', ..., and
  subsets the two attributes 'numvalues' and 'plcoord'.
  This is useful if the way of displaying the axis is to be kept when a
  new plot is drawn.
}
\value{Data.frame with the selected rows and all attributes as described
  above. 
}
\author{Werner A. Stahel
}

\seealso{
  Argument \code{subset} of the high level 'pl' functions
  \code{\link{plyx}, \link{plmatrix}}
}
\examples{
data(d.river)
dd <- d.river[seq(1,1000,4),]
dd$time <- gentimeaxis("date",hour="hour", data=dd)
dsubs <- plsubset(dd, 1:50)
attr(dsubs$time, "ticksat")

plyx(O2~time, data=dsubs)
## same as
## plyx(O2~time, data=dd, subset=1:50)
}
\keyword{manip}
\keyword{attribute}

