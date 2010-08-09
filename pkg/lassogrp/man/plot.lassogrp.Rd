\name{plot.lassogrp}
\Rdversion{1.1}
\alias{plot.lassogrp}
\title{Plotting Method for 'lassogrp' Objects}
\description{
  Plots either the traces of the coefficients, the penalty of the
  coefficient groups, or the penalty terms of the whole model obtained
  by lasso regression as functions of the penalty weight lambda.
}
\usage{
\method{plot}{lassogrp}(x, type = c("norms", "coefficients", "criteria"),
    col = NULL, lty = NULL, mar = NULL, main = NULL,
    cv = NULL, se = TRUE, ylim = NULL, legend = TRUE, ...)
}
\arguments{
  \item{x}{object of class 'lassogrp' containing the results of fitting
    a model by lasso.}
  \item{type}{character string defining the type of results to be plotted:
    \describe{
      \item{"norms"}{For simple (non-grouped) regressors, show the
	absolute value of the estimated coefficient.  For grouped
	variables, show the (L2) norms of the coefficient vector.  In
	both cases, this is the contribution to the L1 penalty term of
	the (grouped) lasso.}
      \item{"coefficients"}{the estimated coefficients}
      \item{"criteria"}{the (average) log-likelihood and the L1 penalty
	term}
    }
  }
  \item{col, lty}{colors and line types for the curves.
    If \code{type} is \code{"criteria"}, they are used as follows
    \describe{
      \item{[1]}{for the log-likelihood,}
      \item{[2]}{for the cross-validated root mse (!!!)}
      \item{[3]}{for the cross-validated standard error band
	of the root mse (!!!)}
      \item{[4]}{for the L1 penalty term.}
    }
  }
  \item{mar}{plot margins}
  \item{main}{main title of the plot}
  \item{cv}{results of cross validation (see \code{\link{cv.lasso}}) to
    be plotted, or logical: if true, \code{cv.lasso} will be called and
    cross validation results will be plotted.}
  \item{se}{logical: Should the cross validation standard error band
    be plotted?}
  \item{ylim}{limits of the vertical plotting axis.
    Defaults to the range of the values to be plotted except for the
    values for lambda==0 (if these are contained in \code{x})}
  \item{legend}{logical: should a legend be displayed?}
  \item{\dots}{further arguments passed to \code{matplot}}
}
\value{none}
\author{Werner A. Stahel, ETH Zurich, \email{stahel@stat.math.ethz.ch}}
\note{For model other than ordinary regression, the cross validation
  results are not yet treated in a sensible way.
}
\examples{
data(asphalt)
rr <- lasso(log10(RUT) ~ log10(VISC)+ASPH+BASE+FINES+VOIDS+RUN,
            data=asphalt)
plot(rr)
plot(rr, type="criteria")
}
\keyword{hplot}
\keyword{regression}

