\name{transferAttributes}
\alias{transferAttributes}

\title{Transfer Attributes
}
\description{Attach the attributes of an object to another object
}
\usage{
transferAttributes(x, xbefore)
}

\arguments{
  \item{x}{the object to which the attributes should be transferred
  }
  \item{xbefore}{the object which delvers the attributes
  }
}
\value{Object \code{x} with attributes from \code{xbefore} (and possibly
  some that it already had)
}
\author{W. A. Stahel
}
\note{This function would not be needed if \code{structure}
  allowed for a list of attributes.
}
\seealso{\code{\link{structure}}
}

\examples{
a <- structure(1:10, title="sequence")
transferAttributes(31:40, a)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}% use one of  RShowDoc("KEYWORDS")

