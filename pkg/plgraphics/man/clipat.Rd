\name{clipat}
\alias{clipat}
\title{Clip Data Outside a Range
}
\description{Drop values outside a given range
}
\usage{
clipat(x, range=NULL, clipped=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of data to be clipped at \code{range}
  }
  \item{range}{range, a numerical vector of 2 elements
  }
  \item{clipped}{if \code{NULL}, the clipped data will be dropped.
    Otherwise, they will be replaced by \code{clipped}, which is
    typically set to \code{NA}.
    If \code{clipped} is numerical of length 2, the elements of
    \code{x} clipped below are set to \code{clipped[1]},
    those clipped by \code{range[2]}, by \code{clipped[2]}.
    Therefore, if \code{clipped} equals \code{range},
    \code{x} will be "Winsorized".
  }
}
\value{As the input \code{x}, with pertinent elements dropped or replaced 
}
\author{Werner A, Stahel
}
\examples{
clipat(rnorm(10,8,2), c(10,20), clipped=NA)
}
\keyword{arith}

