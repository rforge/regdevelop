\name{twosamples}
\alias{twosamples}
\alias{onesample}
\alias{twosamples.default}
\alias{twosamples.formula}

\title{Relevance and Significance for One or Two Samples}
\description{Inference for a difference between two independent samples
  or for a single sample: Collect quantities for inference, 
  including Relevance and Significance measures
}
\usage{
twosamples(x, ...)
onesample(x, ...)

\method{twosamples}{default}(x, y = NULL, paired = FALSE,
  hypothesis = 0,var.equal = TRUE,
  testlevel=getOption("testlevel"), rlv.threshold=getOption("rlv.threshold"), ...)
\method{twosamples}{formula}(formula, data, subset, na.action, ...)
}
\arguments{
  \item{x}{a formula or the data for the first or the single sample
  }
  \item{y}{data for the second sample
  }
  \item{paired}{logical: In case \code{x} and \code{y} are given.
    are their values paired?
  }
  \item{hypothesis}{the null effect to be tested, and anchor for the
    relevance
  }
  \item{var.equal}{logical: In case of two samples,
    should the variances be assumed equal?
    Only applies for quantitative data.
  }
  \item{testlevel}{level for the test,
    also determining the confidence level
  }
  \item{rlv.threshold}{Relevance threshold, or a vector of thresholds
    from which the element \code{stand} is taken for quantitative data
    and the element \code{prop}, for binary data.
  }

  For the \code{formula} method:

  \item{formula, data}{formula and data for the \code{formula} method
  }
  \item{subset, na.action}{subset and na.action to be applied to
    \code{data}
  }
  \item{...}{further arguments. For the \code{formula} method,
    this includes the arguments \code{var.equal, testlevel} and
    code{rlv.threshold} of the \code{default} method.
  }
}
\value{an object of \code{\link{class}} \code{'effecttable'}, a
  \code{\link{data.frame}} with columns (components)
  \describe{
    \item{\code{effect}: }{for quantitative data:
      estimated difference between expectations of the two
      samples, or mean in case of a single sample.

      For binary data: log odds (for one sample or paired samples) or
      log odds ratio (for two samples)}
    \item{\code{ciLow, ciUp}: }{confidence interval for the effect}
    \item{\code{Rle, Rls, Rlp}: }{relevance measures: estimated, secured, potential}
    \item{\code{Sig0}: }{significance measure for test or 0 effect}
    \item{\code{Sigth}: }{significance measure for test of
      \code{effect}  == relevance threshold}
    \item{\code{p.value}: }{p value for test against 0}
  }
  In addition to the columns/components, it has \code{\link{attributes}}
  \describe{
    \item{\code{method}: }{problem and inference method}
    \item{\code{effectname}: }{label for the effect}
    \item{\code{hypothesis}: }{the null effect}
    \item{\code{n}: }{number(s) of observations}
    \item{\code{means}: }{in the case of 2 independent samples: their means}
    \item{\code{statistic}: }{test statistic}
    \item{\code{V}: }{single observation variance}
    \item{\code{df}: }{degrees of freedom for the t distribution}
    \item{\code{rlv.threshold}: }{relevance threshold}
  }
}
\note{\code{onesample} and \code{twosamples} are identical.
  }
\references{see those in \code{\link{relevance-package}}.}

\author{Werner A. Stahel}

\seealso{\code{\link{t.test}, \link{binom.test}}
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

## one sample
onesample(sleep[sleep$group == 2, "extra"])
}
\keyword{htest}
\concept{relevance}
