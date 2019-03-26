\name{nobs.survreg}
\alias{nobs.survreg}
\alias{nobs.coxph}

\title{Number of Observations
}
\description{Methods of \code{nobs} for classes
  \code{survreg} and \code{coxph} (since they are lacking presently)
}
\usage{
\method{nobs}{survreg}(object, ...)
\method{nobs}{coxph}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{result of calling \code{survreg} or \code{coxph}
  }
  \item{\dots}{not used
  }
}
\value{scalar: number of non-missing observations
}
\seealso{\code{\link{nobs}}
}
\examples{
data(ovarian, package="survival")
rr <- survival::survreg(survival::Surv(futime, fustat) ~ ecog.ps + rx, ovarian,
             dist="exponential")
nobs(rr)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}% use one of  RShowDoc("KEYWORDS")

