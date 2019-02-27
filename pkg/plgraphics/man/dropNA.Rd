\name{dropNA}
\alias{dropNA}

\title{drop NA values from vector
}
\description{Return the vector 'x', without elements that are NA or NaN
  or, if 'inf' is TRUE, equal to Inf or -Inf  
}
\usage{
dropNA(x, inf = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector from which the non-real values should be dropped
}
  \item{inf}{Should 'Inf' and '-Inf' be considered "non-real"?
}
}
\value{Vector containing the 'real' values of 'x' only
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
\code{\link{na.omit}, \link{sumNA}}
}
\examples{
dd <- c(1, NA, 0/0, 4, -1/0, 6)
dropNA(dd)
na.omit(dd)
}
\keyword{manip}
\keyword{NA}
