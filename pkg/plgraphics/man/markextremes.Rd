\name{markextremes}
\alias{markextremes}

\title{
  Adjust the default proportion of extreme points to be labeled to the number
  of observations
}
\description{
  Adjusts the proportion of extreme points to be labeled to the number of
  observations. It is the default of the \code{ploption markextremes}. 
}
\usage{
markextremes(n)
}

\arguments{
  \item{n}{number of observations
  }
}
\details{
  The function simply applies \code{ceiling(sqrt(n)/2)/n}.
}
\value{A scalar between 0 and 0.5
}
\author{Werner A. Stahel
}
\examples{
markextremes(20)
for (n in c(10,20,50,100,1000))  print(c(n,markextremes(n)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{auxiliary}% use one of  RShowDoc("KEYWORDS")

