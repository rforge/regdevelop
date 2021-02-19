\name{shift}
\alias{shift}

\title{Shift a Vector, Matrix or Data.frame
}
\description{Shifts a vector by introducing NA's at the beginning and
  dropping elements at the end or vice versa. Analogous for rows of
  matrices and data.frames
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
shift(x, k = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector or matrix to be shifted
  }
  \item{k}{number of elements or rows to be dropped in the end (if
    \code{k} is positive) or at the beginning (otherwise)
  }
}
\value{shifted vector or matrix
}
\note{This function helps to regress a target variable on lagged
  (=shifted) versions of itself (autocorrelation) and/or explanatory
  variables (time series regression)
}
\author{Werner A. Stahel
}
\seealso{\code{\link{lag}} for \code{ts} objects
}
\examples{
shift(1:5,k=2)
shift(c("a","b","c"),-1)

## need an example for time series regression!!!
}
\keyword{manip}
