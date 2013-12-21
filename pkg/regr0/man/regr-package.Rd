\name{regr-package}
\alias{regr-package}
\docType{package}
\title{
  Tools for Building Many Kinds of Regression Models
}
\description{
  This package helps to perform efficient model building for many kinds
  of regression models, including (generalized) linear, robust linear,
  ordinal (aka \code{\link[MASS]{polr}}) etc.
  It also contains some useful utility functions.
}
\details{
\tabular{ll}{
Package: \tab regr\cr
Type: \tab Package\cr
Date: \tab >= 2008-08-20\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

The principal function \code{\link{regr}()} allows fitting regression
models of many types by calling well-known \R functions for
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
of \code{regr} (see \code{\link{plot.regr}}), which gives more complete
regression diagnostic plots than the specific plot methods of the
mentioned \dQuote{base} \R functions.

An important part of residual analysis which is often neglected is
plotting residuals versus explanatory variables, often called \emph{term
plots}.  This is part of \code{\link{plot.regr}}, but can also be called
independently by calling \code{\link{plresx}}.

If a \code{regr} result is printed, the output is more informative than
the classical tables.

The package contains the methods for the usual extractor functions, like
\code{residuals}, \code{predict}, ..., as well as for \code{drop1},
\code{add1}, etc, ...

The package contains a number of additional functions that should be
useful for efficient data analysis. They include
\itemize{
  \item \code{\link{plresx}()} to plot residuals against explanatory variables,
  \item \code{\link{plres2x}()} to produce a plot showing interactions between two
  continuous regressors,
  \item \code{\link{plmatrix}()}, a version of \code{\link{pairs}} that
  displays any rectangular scatterplot matrix and can also be used to
  produce the full usual quadratic scatterplot matrix for many variables
  in appropriate junks with panels that are still large enough to be useful.
  %% Gibt's nicht: \item \code{\link{lassoselect}()} to do model selection.
}

The package also includes some additional utility functions, like
\itemize{
\item \code{\link{logst}()} for log transformation of variables with zeros,
\item \code{\link{showd}()} for printing a "representative" part of a dataset.
}

For documenting the steps of the data analysis, there are several
documenting functions:

\code{\link{tit}} and \code{\link{doc}} add attributes to data and results, which
should allow for tracing the history of these objects (not yet
thoroughly implemented).

\code{\link{stamp}} adds tracking information to plots, consisting of a project
title, a step label, and the date when the plot is produced.

Some default options are added to \code{\link{options}} to reflect user
preferences for output, see the help entries for
\code{\link{stamp}}. % and \code{\link{print.regr}}.

}
\author{
Werner A. Stahel, ETH Zurich, \email{stahel@stat.math.ethz.ch}
}
\references{
  More information is available at \url{http://stat.ethz.ch/~stahel/regression}
}
\keyword{ package }
\seealso{
  The \pkg{car} package.
}
\examples{
  r.savings <- regr(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
  plot(r.savings)
}
