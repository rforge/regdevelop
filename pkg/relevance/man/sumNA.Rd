\name{sumNA}
\alias{sumNA}

\title{Count NA Values
}
\description{Count NA Values
}
\usage{
sumNA(object, inf = TRUE)
}

\arguments{
  \item{object}{a vector, matrix, or data.frame
}
  \item{inf}{logical: Shall infinite values be treated as NAs?
}
}
\value{Number of \code{NA}s -- for each column if \code{object} is a
  matrix or data.frame
}
\author{Werner A. Stahel
}

\seealso{
  \code{\link{dropNA}}
}
\examples{
 sumNA(c(0,3,NA,5,Inf,NA))  ## 3
 sumNA(c(0,3,NA,5,Inf,NA), inf=FALSE)  ## 2

 data(d.blast)
 sumNA(d.blast)  ##  26 NAs in variable 'tremor'
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{attribute}% use one of  RShowDoc("KEYWORDS")

