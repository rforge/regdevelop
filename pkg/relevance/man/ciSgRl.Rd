\name{ciSgRl}
\alias{ciSgRl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create Confidence Intervals and Significance Values for a Termtable
}
\description{Calculates confidence intervals and significance values for
  the coefficients of a fitted model to be used in a
  \code{termtable} as produced by \code{regr}.
}
\usage{
ciSgRl(coef, se = NULL, df = Inf, testlevel = 0.05, stcoef=TRUE, rlv=TRUE, 
            rlvThres=0.1, object=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coef}{vector of coefficients or matrix containing
    coefficients and their standard errors. The latter is needed if
    \code{se} is not given.
  }
  \item{se}{vector of standard errors of the coefficients
  }
  \item{df}{degrees of freedom of the residuals
  }
  \item{testlevel}{test level
  }
  \item{stcoef}{standardized coefficients.
    If \code{NULL}, these will be calculated from \code{object}.
  }
  \item{rlv}{logical: Should relevances be calculated?
  }
  \item{rlvThres}{Relevance threshold
  }
  \item{object}{model fit, needed for calculating \code{stcoef}.
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
  If \code{rlv}
}
\author{Werner A. Stahel
}
%%- \seealso{
%%- %% ~~objects to See Also as \code{\link{help}}, ~~~
%%- }
\examples{
data(d.blast)
rr <-
  lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast) 
ciSgRl(rr) %% same as
ciSgRl(summary(rr)$coef)
}
\keyword{ regression }% use one of  RShowDoc("KEYWORDS")

