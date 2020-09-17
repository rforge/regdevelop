\name{rlstats}
\alias{rlstats}

\title{Calculate Relevance and Significance Statistics
}
\description{Calculate relevance and significance statistics
  for samples or two group comparisons
}
\usage{
rlstats(x = NULL, se = NULL, df = Inf, testlevel = NULL,
  steff = NULL, rlv = TRUE, rlvthres = NULL, object = NULL)
}

\arguments{
  \item{x}{estimated effect(s) 
  }
  \item{se}{their standard errors
  }
  \item{df}{degrees of freedom for the t test
  }
  \item{testlevel}{level of the statistical test, usually 0.05
  }
  \item{steff}{vector of standardized effects, or \cr
    logical: should standardized coefficients be calculated?
    This only works if \code{object} is a regression result
  }
  \item{rlv}{logical: should relevances be calculated?
  }
  \item{rlvthres}{relevance threshold. May be a vector containing
    a component names \code{steff}
  }
  \item{object}{a regression fit, needed in case standardized effects
    (coefficients) need to be calculated.
  }
}
\details{
to be described  
}
\value{List with components ...
}
\references{
  Stahel, Werner A. (2020). 
  Measuring Significance and Relevance instead of p-values.
  In preparation.
}
\author{\packageAuthor{relevance}

Maintainer: \packageMaintainer{relevance}
}

\examples{
  rlstats(5, se=3, rlv=FALSE)
}
\keyword{univariate}

