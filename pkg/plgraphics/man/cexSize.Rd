\name{cexSize}
\alias{cexSize}

\title{
Adjust character size to number of observations
}
\description{
Adjusts the character size \code{cex} to number of observations
}
\usage{
cexSize(n)
}

\arguments{
  \item{n}{number of observations
}
}
\details{
  The function simply applies \code{min(1.5/log10(n), 2)}
}
\value{A scalar, defining \code{cex}
}
\author{Werner A. Stahel
}
\examples{
cexSize(20)
for (n in c(10,20,50,100,1000))  print(c(n,cexSize(n)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
