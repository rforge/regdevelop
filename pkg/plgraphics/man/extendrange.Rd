\name{extendrange}
\alias{extendrange}
\alias{clipat}
\title{extend a range, clip data outside a range
}
\description{Extend a Range on Both Ends by a Proportion
}
\usage{
extendrange(range, ext=0.05)
clipat(x, range=NULL, clipped=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{range}{range, a numerical vector of 2 elements
  }
  \item{ext}{proportion by which the range is to be extended.
    Can also be a vector of 2 proportions, for the respective
    ends of the range.
  }
  \item{x}{vector of data to be clipped at \code{range}
  }
  \item{clipped}{if \code{NULL}, the clipped data will be dropped.
    Otherwise, they will be replaced by \code{clipped}, which is
    typically set to \code{NA}
  }
}
\value{\code{extendrange}:
  Extended range, a numerical vector of length 2 \cr
  \code{clippedat}: as the input \code{x}, with pertinent elements
    dropped or replaced 
}
\author{Werner A, Stahel
}
\examples{
## extendrange is just
## rg + c(-1, 1) * ext * diff(rg)
extendrange(c(10,20), ext=0.08)

## clipat
clipat(rnorm(10,8,2), c(10,20), clipped=-1)
}
\keyword{arith}

