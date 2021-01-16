\name{termtable}
\alias{termtable}
\alias{print.termtable}

\title{Statistics for Linear Models, Including Relevance Statistics}
\description{Calculate a table of statistics for (multiple) regression
  mdels with a linear predictor and respective \code{print} method
}
\usage{
termtable(object, summary = NULL, testtype = NULL, r2x = TRUE, rlv = TRUE,
          rlvThres = c(rel = 0.1, coef = 0.1, drop = 0.1, pred = 0.05),
          testlevel = 0.05)

\method{print}{termtable}(x, columns = NULL, printstyle = NULL, legend = NULL,
      digits = NULL, na.print = "  ", ...)
}

\arguments{
  \item{object}{result of a model fitting function like \code{\link{lm}}
  }
  \item{summary}{result of \code{summary(object)}.  If \code{NULL}, the
    \code{summary} will be called.
  }
  \item{testtype}{type of test to be applied for dropping each term in
    turn. If \code{NULL}, it is selected according to the class of the
    object, see Details.
  }
  \item{r2x}{logical: should the collinearity measures \code{"R2.x"}
    for the terms be calculated?
  }
  \item{rlv}{logical: Should relevances be calculated?
  }
  \item{rlvThres}{Relevance thresholds, vector with the three components
    \cr{\code{rel}:}  {threshold for relative effects,}
    \cr{\code{coef}:}  {for standardized coefficients,}
    \cr{\code{drop}:}  {for drop effects,}
    \cr{\code{pred}:}  {for prediction intervals.}
  }
  \item{testlevel}{1 - confidence level}

  %%-- print method:
  \item{x}{the \code{termtable} object to be printed.}
  \item{columns}{.. ..}%--- TODO: replace all '.. ..' <<<<<<<<<<<<<<<<<
  \item{printstyle}{.. ..}
  \item{legend}{.. ..}
  \item{digits}{.. ..}
  \item{na.print}{.. ..}
  \item{...}{potentially further arguments to \code{\link{print}()} methods.}
}
\details{
  If \code{testtype} is not specified, it is determined by the class of
  \code{object} and its attribute \code{family} as follows:\cr
  --- {F}{or t for objects of class \code{lm, lmrob} and \code{glm}
    with families \code{quasibinomial} and \code{quasipoisson},}\cr
  --- {Chi-squared}{for other \code{glm}s and \code{survreg}}\cr
}
\value{\code{\link{data.frame}} with columns
  \describe{
  \item{\code{coef}: }{coefficients for terms with a single degree of freedom}
  \item{\code{df}: }{degrees of freedom}
  \item{\code{se}: }{standard error of \code{coef}}
  \item{\code{ciLow, ciUp}: }{confidence interval for \code{coef}}
  \item{\code{Sig0}: }{significance value for the test of \code{coef==0}}
  \item{\code{stcoef}: }{standardized coefficient (standardized using
    the standard deviation of the 'error' term, \code{sigma},
    instead of the response's standard deviation)}
  \item{\code{st.Low, st.Up}: }{confidence interval for \code{stcoef}}
  \item{\code{testst}: }{test statistic}
  \item{\code{p.value, p.symbol}: }{p value and symbol for it}
  \item{\code{R2.x}: }{collinearity measure
    (\eqn{= 1 - 1 / vif}, where \eqn{vif} is the variance inflation factor)}
  \item{\code{coefRle}: }{estimated relevance of \code{coef}}
  \item{\code{coefRls}: }{secured relevance, lower end of confidence interval
    for the relevance of \code{coef}}
  \item{\code{coefRlp}: }{potential relevance, the upper end of the confidence interval.}
  \item{\code{dropRle, dropRls, dropRlp}: }{analogous values for drop effect}
  \item{\code{predRle, predRls, predRlp}: }{analogous values for prediction effect}
  \item{\code{coefRls.symbol, dropRls.symbol, predRls.symbol}: }{symbols for
    the secured relevances}
  }

  In addition, it has \code{\link{attributes}}
  \describe{
    \item{\code{testtype}: }{as determined by the argument \code{testtype} or
      the class and attributes of \code{object}.}
    \item{\code{fitclass, family, dist}: }{class and attributes of \code{object}.}
  }
}
\references{Werner A. Stahel (2020).
  Measuring Significance and Relevance instead of p-values.
  Submitted
}
\author{Werner A. Stahel}
\seealso{
  \code{summary(object)$coefficients}
}
\examples{
  data(swiss)
  rr <- lm(Fertility ~ . , data = swiss)
  rt <- termtable(rr)
  rt
  names(rt)
  print(rt[2:4,], columns="all")
}
\keyword{regression}
\concept{relevance}

