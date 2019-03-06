\name{plsubset}
\alias{plsubset}

\title{Subsetting a Data.Frame with pl Attributes
}
\description{Select rows of data.frames keeping the variable attributes
  that drive pl graphics
}
\usage{
plsubset(x, subset = NULL, omit = NULL, select = NULL, drop = FALSE,
         keeprange = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data.frame from which the subset is to be generated
  }
  \item{subset, omit}{logical vector or vector of indices of rows or or
    rownames of \code{x}.
    If \code{subset} is used, \code{omit} is ignored.
  }
  \item{select}{vector of indices or names of variables to be selected
  }
  \item{drop}{logical: if only one variable remains, should the
    data.frame be converted into a vector?
  }
  \item{keeprange}{logical: should ranges
    (\code{inner.range} and \code{plrange}) be maintained?
  }
}
\details{\code{plsubset} maintains the 'pl' attributes of the variables
  of the data.frame (if there are), such as 'col', 'lty', ..., and
  subsets the two attributes 'numvalues' and 'plcoord'.
  This is useful if the way of displaying the axis is to be kept when a
  new plot is drawn.
}
\value{Data.frame with the selected rows (or without the \code{omit}ted
  rows, respectively) and all attributes as described above. 
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
attr(dd$time, "ticksat")

dsubs <- plsubset(dd, subset=1:50)
attr(dsubs$time, "ticksat")

plyx(O2~time, data=dsubs)
## same as
## plyx(O2~time, data=dd, subset=1:50)
}
\keyword{manip}
\keyword{attribute}

