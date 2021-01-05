\name{twosamples}
\alias{twosamples}
\alias{twosamples.default}
\alias{twosamples.formula}
\title{Compare Two Groups, or Inference for one Sample: Relevance and Significance
}
\description{Inference for a difference between two independent groups
  or for a single sample
  based on the t-test:
  Calculate Relevance and Significance measures
}
\usage{
twosamples(x, ...)
\method{twosamples}{default}(x, y = NULL, paired = FALSE, var.equal = TRUE, 
    testlevel = 0.05, rlvThres = NULL, ...)
\method{twosamples}{formula}(formula, data, subset, na.action, ...)
}

\arguments{
  \item{x, y}{data of the two groups
  }
  \item{paired}{logical: In case \code{x} and \code{y} are given.
    are their values paired?
  }
  \item{var.equal}{logical: Should the variances be assumed equal?
  }
  \item{testlevel}{level for the test used for confidence intervals
  }
  \item{rlvThres}{Relevance threshold
  }
  For the \code{formula} method:
  \item{formula, data}{formula and data for the \code{formula} method
  }
  \item{subset, na.action}{subset and na.action to be applied to
    \code{data}
  }
  \item{...}{further arguments. For the \code{formula} method,
    this includes the arguments \code{var.equal, testlevel} and
    code{rlvThres} of the \code{default} method.
  }
}
\value{List with the components:
  \item{estimate}{estimated difference between expectations of the two
    groups, or mean in case of a single sample, with confidence interval
  }
  \item{rlv}{relevance measures: estimated, secured, potential}
  \item{V}{single observation variance}
  \item{sig0}{significance measure for test on 0 difference}
  \item{sigth}{significance measure for test on difference=rlevance
    threshold}
  \item{p.value}{p value for test against 0}
  \item{n}{number of observations}
  \item{df}{degrees of freedom for the t distribution}
  \item{means}{in the case of 2 independent groups: their means}
}
\references{W. A. Stahel (2021)
}
\author{Werner A. Stahel}

\seealso{\code{\link{t.test}}
}
\examples{
  t.test(sleep[sleep$group == 1, "extra"], sleep[sleep$group == 2, "extra"])
  twosamples(sleep[sleep$group == 1, "extra"], sleep[sleep$group == 2, "extra"])

## Two-sample test, wilcox.test example,  Hollander & Wolfe (1973), 69f.
## Permeability constants of the human chorioamnion (a placental membrane)
## at term and between 12 to 26 weeks gestational age
d.permeabililty <-
  data.frame(perm = c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46,
                      1.15, 0.88, 0.90, 0.74, 1.21), atterm = rep(1:0, c(10,5))
             )
t.test(perm~atterm, data=d.permeabililty)
twosamples(perm~atterm, data=d.permeabililty)
}
\keyword{htest}
\concept{relevance}
