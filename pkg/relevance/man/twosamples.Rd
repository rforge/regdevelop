\name{twosamples}
%zu lang:\title{Compare Two Groups, or Inference for one Sample: Relevance and Significance}
\title{Relevance and Significance for One or Two Samples}
\alias{twosamples}
\alias{twosamples.default}
\alias{twosamples.formula}
\alias{print.effecttable}
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

\method{print}{effecttable}(x, columns=NULL, ...)
}
\arguments{
  \item{x, y}{for \code{twotables()}: data of the two groups.
    \cr for the \code{"effecttable"} \code{print()} method: the effect
      table to be printed.}
  \item{paired}{logical: In case \code{x} and \code{y} are given.
    are their values paired?
  }
  \item{var.equal}{logical: Should the variances be assumed equal?
  }
  \item{testlevel}{level for the test used for confidence intervals
  }
  \item{rlvThres}{Relevance threshold}

  %% [MM] Auskommentiert: Man bringt das nicht gut hin (Text|HTML|PDF) :
  %% \cr
  %% For the \code{formula} method:

  \item{formula, data}{formula and data for the \code{formula} method
  }
  \item{subset, na.action}{subset and na.action to be applied to
    \code{data}
  }
  \item{...}{further arguments. For the \code{formula} method,
    this includes the arguments \code{var.equal, testlevel} and
    code{rlvThres} of the \code{default} method.
  }

  %% for print.effecttable()
  \item{columns}{.. ..}% <<<<< TODO !
}
\value{an object of \code{\link{class}} \code{'effecttable'}, a
  \code{\link{data.frame}} with columns (components)
  \describe{
    \item{\code{effect}: }{estimated difference between expectations of the two
      groups, or mean in case of a single sample, with confidence interval}
    \item{\code{ciLow, ciUp}: }{confidence interval for the effect}
    \item{\code{Rle, Rls, Rlp}: }{relevance measures: estimated, secured, potential}
    \item{\code{Rls.symbol}: }{symbol describing \code{Rls}}
    \item{\code{Sig0}: }{significance measure for test or 0 effect}
    \item{\code{Sigth}: }{significance measure for test of \code{effect}  == relevance threshold}
    \item{\code{p.value}: }{p value for test against 0}
    \item{\code{p.symbol}: }{symbol for \code{p.value}}
  }
  In addition to the columns/components, it has \code{\link{attributes}}
  \describe{
    \item{\code{rlvThres}: }{relevance threshold}
    \item{\code{n}: }{number(s) of observations}
    \item{\code{df}: }{degrees of freedom for the t distribution}
    \item{\code{means}: }{in the case of 2 independent groups: their means}
    \item{\code{V}: }{single observation variance}
    \item{\code{pSymbols, rlvSymbols}: }{lists defining the symbols used for
      p-values and relevances}
  }
}
\references{see those in \code{\link{relevance-package}}.}

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
