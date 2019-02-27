\name{ciSignif}
\alias{ciSignif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create Confidence Intervals and Significance Values for a Termtable
}
\description{Calculates confidence intervals and significance values for
  the coefficients of a fitted model to be used in a
  \code{termtable} as produced by \code{regr}.
}
\usage{
ciSignif(estimate, se = NULL, df = Inf, testlevel = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estimate}{vector of coefficients or matrix containing
    coefficients and their standard errors. The latter is needed if
    \code{se} is not given.
  }
  \item{se}{vector of standard errors of the coefficients
  }
  \item{df}{degrees of freedom of the residuals
  }
  \item{testlevel}{test level
  }
}
\details{The coefficients divided by standard errors are assumed to be
  t-distributed with \code{df} degrees of freedom.
  For \code{df==Inf}, this is the standard normal distribution.
}
\value{
  A data.frame with the variables
  \item{estimate, se}{coefficients and standard errors}
  \item{ciLow, ciUp}{lower and upper limit of the confidence interval}
  \item{teststatistic}{t-test statistic}
  \item{signif}{significance value, i.e., test statistic divided by
    critical value, which in turn is the \code{1-testlevel/2}-quantile
    of the t-distribution.
  }
  \item{p.value}{p value}
  \item{p.symbol}{the conventional symbol corresponding to the p value}
}
\author{Werner A. Stahel
}
%%- \seealso{
%%- %% ~~objects to See Also as \code{\link{help}}, ~~~
%%- }
\examples{
data(d.blast, package="plgraphics")
rr <-
  lm(logst(tremor)~location+log10(distance)+log10(charge), data=d.blast) 
ciSignif(rr) %% same as
ciSignif(summary(rr)$coef)
}
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")

