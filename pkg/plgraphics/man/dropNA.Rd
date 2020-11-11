\name{dropNA}
\alias{dropNA}
\alias{replaceNA}

\title{drop or replace NA values
}
\description{\code{dropNA} returns the vector 'x', without elements that are NA or NaN
    or, if 'inf' is TRUE, equal to Inf or -Inf.
    \code{replaceNA} replaces these values by values from the second argument
}
\usage{
dropNA(x, inf = TRUE)
replaceNA(x, na, inf = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector from which the non-real values should be dropped
    or replaced
  }
  \item{na}{replacement or vector from which the replacing values
    are taken.
  }
  \item{inf}{logical: should 'Inf' and '-Inf' be considered "non-real"?
  }
}
\value{For \code{dropNA}: Vector containing the 'real' values
  of 'x' only \cr
  For \code{replaceNA}: Vector with 'non-real' values replaced by
  the respective elements of \code{na}.
}
\author{
Werner A. Stahel
}
\note{The differences to 'na.omit(x)' are:
  'Inf' and '-Inf' are also dropped, unless 'inf==FALSE'.\\
  no attribute 'na.action' is appended.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{na.omit}, \link{sumNA}, \link{ifelse}}
}
\examples{
dd <- c(1, NA, 0/0, 4, -1/0, 6)
dropNA(dd)
na.omit(dd)

replaceNA(dd, 99)
replaceNA(dd, 100+1:6)
}
\keyword{manip}
\keyword{NA}
