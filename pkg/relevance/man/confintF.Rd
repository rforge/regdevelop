\name{confintF}
\alias{confintF}

\title{
Confidence Interval for the Non-Central F and Chisquare Distribution
}
\description{
Confidence Interval for the Non-Central F and Chisquare Distribution
}
\usage{
confintF(f, df1, df2, testlevel = 0.05)
}

\arguments{
  \item{f}{observed F value(s)
  }
  \item{df1}{degrees of freedom for the numerator of the F distribution
  }
  \item{df2}{degrees of freedom for the denominator of the F distribution
  }
  \item{testlevel}{level of the (two-sided) test that determines the
    confidence interval, 1 - confidence level
  }
}
\details{
  The confidence interval is calculated by solving the two implicit
  equations  \code{qf(f, df1, df2, x) = testlevel/2}  and
  \code{... = 1 - testlevel/2}.
  For \code{f>100}, the usual \code{f +- standard error} interval is
  used as a rather crude approximation.

  A confidence interval for the non-centrality of the Chisquare
  distribution is obtained by setting \code{df2} to \code{Inf}
  (the default) and \code{f=x2/df1} if \code{x2} is the observed
  Chisquare value.  
}
\value{vector of lower and upper limit of the confidence interval,
  or, if any of the arguments has length \code{>1}, matrix containing
  the intervals as rows.
}
\author{\packageAuthor{relevance}}

\seealso{\code{\link{qf}}
}
\examples{
confintF(5, 3, 200)
## [1] 2.107 31.95
confintF(1:5, 5, 20)   ## lower limit is 0 for the first 3 f values
}
\keyword{distribution}
\keyword{univar}
