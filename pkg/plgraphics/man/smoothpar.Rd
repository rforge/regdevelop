\name{smoothpar}
\alias{smoothpar}

\title{
  Adjust the smoothing parameter to number of observations
}
\description{Adjust the smoothing parameter to number of observations
}
\usage{
smoothpar(n)
}

\arguments{
  \item{n}{number of observations
  }
}
\value{smoothing parameter
}
\author{Werner A. Stahel
}
\examples{
  smoothpar(50)
  t.n <- c(5,10,20,100,1000)
  smoothpar(t.n)
}
\keyword{auxiliary} % use one of  RShowDoc("KEYWORDS")
