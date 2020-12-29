\name{twogroups}
\alias{twogroups}
\title{Compare Two Groups: Relevance and Significance
}
\description{Test for a difference between two independent groups
  based on the t-test:
  Calculate Relevance and Significance measures
}
\usage{
twogroups(x, y, var.equal = TRUE, testlevel = 0.05, RelTh = 0.25)
}

\arguments{
  \item{x, y}{data of the two groups
  }
  \item{var.equal}{logical: Should the variances be assumed equal?
  }
  \item{testlevel}{level for the test used for confidence intervals
  }
  \item{RelTh}{Relevance threshold
  }
}
\value{List with the components:
  \item{estimate}{estimated difference between expectations of the two groups}
  \item{rlv}{relevance measures: }
  \item{V}{single observation variance}
  \item{sig0}{significance measure for test on 0 difference}
  \item{sigth}{significance measure for test on difference=rlevance
    threshold}
}
\references{W. A. Stahel (2021)
}
\author{Werner A. Stahel}

\seealso{\code{\link{t.test}}
}
\examples{
  t.test(sleep[sleep$group == 1, "extra"], sleep[sleep$group == 2, "extra"])
  twogroups(sleep[sleep$group == 1, "extra"], sleep[sleep$group == 2, "extra"])
}
\keyword{htest}
\concept{relevance}
