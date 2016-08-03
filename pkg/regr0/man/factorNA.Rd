\name{factorNA}
\alias{factorNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
For a factor, onvert NA to a separate level 
}
\description{
Introduces a new level for missing values and assigns it to all missing
values of a factor
}
\usage{
factorNA(x, na.label = ".NA.", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a factor -- or a variable that will be converted to a factor
}
  \item{na.label}{label to be used for the new level of the factor
}
  \item{\dots}{other arguments, to be passed to the function \code{factor}
}
}
\value{
  The factor with the new level. There will be no NAs left.
}
\author{
Werner A. Stahel
}

\seealso{
\code{\link{xNA}}, \code{\link{factor}}, \code{\link{ordered}}
}
\examples{
factorNA(c("a","c",NA,"a",NA), na.label="missing")
}
\keyword{manip}

