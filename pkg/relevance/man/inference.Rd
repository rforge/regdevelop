\name{inference}
\alias{inference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate Confidence Intervals and Relevance and Significance Values
}
\description{Calculates confidence intervals and relevance and
  significance values given estimates, standard errors and, for
  relevance, additional quantities. 
}
\usage{
inference(estimate, se = NULL, df = NULL, testlevel = 0.05, stcoef=TRUE, rlv=TRUE, 
            rlv.threshold=0.1, object=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estimate}{vector of estimates or matrix containing
    estimates and their standard errors. The latter is needed if
    \code{se} is not given.
  }
  \item{se}{vector of standard errors of the estimates
  }
  \item{df}{degrees of freedom of the residuals
  }
  \item{testlevel}{test level
  }
  \item{stcoef}{standardized coefficients, suitable for 
    If \code{NULL}, these will be calculated from \code{object}.
  }
  \item{rlv}{logical: Should relevances be calculated?
  }
  \item{rlv.threshold}{Relevance threshold
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
  \item{testst}{t-test statistic}
  \item{Sig0}{significance value, i.e., test statistic divided by
    critical value, which in turn is the \code{1-testlevel/2}-quantile
    of the t-distribution.
  }
  \item{p.value}{p value}
  \item{p.symbol}{the conventional symbol corresponding to the p value}
  \cr
  If \code{rlv} is \code{TRUE},
  \item{stcoef}{standardized coefficient}
  \item{st.Low, st.Up}{confidence interval for \code{stcoef}}
  \item{Rle}{estimated relevance of \code{coef}}
  \item{Rls}{secured relevance, lower end of confidence interval
    for the relevance of \code{coef}}
  \item{Rlp}{potential relevance, upper end of confidence interval
    ...}
  \item{Rls.symbol}{symbols for the secured relevance}
}
\author{Werner A. Stahel}

\examples{
data(d.blast)
rr <-
  lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast) 
inference(rr) ## same as
inference(summary(rr)$coef)
}
\keyword{ regression }
\concept{relevance}
