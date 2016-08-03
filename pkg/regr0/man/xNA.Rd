\name{xNA}
\alias{xNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
For quantitative variables, make NAs available
}
\description{
  Introduces a new indicator (binary) variable for missing values and
  assigns a fixed, non-missing value for them in the (renamed) original
  variable 
}
\usage{
  xNA(data, na.values = NULL, name.suffix=c(".X",".NA"))
}

\arguments{
  \item{data}{a data.frame of quantitative variables. A vector will be
  converted to a one column data.frame
}
  \item{na.values}{label to be used for the new level of the factor
}
  \item{name.suffix}{suffix to be used for the new quantitative
  variables and the new indicator variables
}
}
\value{
  data.frame, containing, for each variable \code{Vx}, a variable
  \code{Vx.X} and a variable \code{Vx.NA} (if the default
  \code{name.suffix} is used). Variables that do not contain any
  NAs will be tranferred as they are to the result.
}
\author{
Werner A. Stahel
}

\seealso{
\code{\link{factorNA}}
}
\examples{
xNA(c(3,NA,NA,6,5,9,1,11))
df <- data.frame(a=c(3,NA,NA,6,5,9,1,11), b=c(1,2,1,0,3,4,2,0))
( dfn <- xNA(df) )
attr(dfn, "xNA.values")
}
\keyword{manip}

