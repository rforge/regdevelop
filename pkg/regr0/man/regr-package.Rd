\name{regr-package}
\alias{regr-package}
\docType{package}
\title{
  Tools for Building Regression Models
}
\description{
  This package helps to perform efficient model building for regression
  models. It also contains some useful utility functions.

}
\details{
\tabular{ll}{
Package: \tab regr\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-08-20\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
The function \code{regr} is a wrapper function for fitting regression
models of many types by calling well-known R functions:
\itemize{
\item
Ordinary linear regression,
\item
logistic, Poisson, other generalized linear models, also for
\item
ordinal or nominal response variable.
}

The output of the function is somewhat different from the output of the
base R functions \code{lm}, \code{glm}, \code{polr} etc. -- it should
be more informative for the task of model building.

The next most important function is the \code{plot} method for results
of \code{regr}, which gives more complete regression diagnostic plots
than the specific plot methods of the mentioned incorporated R functions.

An important part of residual analysis which is often neglected is
plotting residuals versus explanatory variables, often called term
plots. This is part of \code{plot.regr}, but can also be called
independently by calling \code{plresx}.

If a \code{regr} result is printed, the output is more informative than
the classical tables.

The package contains the methods for the usual extractor functions, like
\code{residuals}, \code{predict}, ..., as well as for \code{drop1},
\code{add1}, ...

The package contains a number of additional functions that should be
useful for efficient data analysis. They include
\itemize{
\item
\code{plresx} to plot residuals against explanatory variables,
\item
\code{plres2x} to produce a plot showing interactions between two
continuous regressors,
\item
\code{plmatrix}, a version of \code{pairs} that displays any rectangular
scatterplot matrix and can also be used to produce the full usual
quadratic scatterplot matrix for many variables in appropriate junks
with panels that are still large enough to be useful.
\item
\code{lassoselect} to do model selection.
}

The package also includes some additional utility functions, like
\itemize{
\item
\code{logst} for log transformation of variables with zeros,
\item
\code{showd} for printing a "representative" part of a dataset.
}

For documenting the steps of the data analysis, there are several
documenting functions:

\code{tit} and \code{doc} add attributes to data and results, which
should allow for tracing the history of these objects (not yet
thoroughly implemented).

\code{stamp} adds tracking information to plots, consisting of a project
title, a step label, and the date when the plot is produced.

Some default options are added to \code{options} to reflect user
preferences for output, see the help entries for
\code{\link{stamp}}. % and \code{\link{print.regr}}.

}
\author{
Werner A. Stahel, ETH Zurich, \email{stahel@stat.math.ethz.ch}
}
\references{
There will be a website on \url{http://stat.ethz.ch/~stahel/regression}
}
\keyword{ package }
\seealso{\code{library(car)}
}
\examples{
  r.savings <- regr(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
  plot(r.savings)
}
