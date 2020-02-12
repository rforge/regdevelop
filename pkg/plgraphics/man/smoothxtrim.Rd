\name{smoothxtrim}
\alias{smoothxtrim}

\title{
  Adjust range for smooth lines to number of observations
}
\description{
  The range in which smooth lines are drawn should be restricted in
  order to avoid the ill determined parts at both ends.
  The proportion of suppressed values is determined as a function of
  the number of observations.
}
\usage{
smoothxtrim(n, c=2)
}

\arguments{
  \item{n}{number of observations
  }
  \item{c}{tuning parameter: how rapidly should the result decrease
    with \code{n}?
  }
}
\value{proportion of x values for which the smoothline will not be shown
  on both ends. Equals \ 1.6^(log10(n)*c) / n
}
\author{W. Stahel
}
\examples{
  smoothxtrim(50)
  t.n <- c(5,10,20,100,1000)
  t.n * smoothxtrim(t.n)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{auxiliary}% use one of  RShowDoc("KEYWORDS")

