\name{extract.lassogrp}
\Rdversion{1.1}
\alias{extract.lassogrp}
\alias{[.lassogrp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Regression Results from a Lasso Fit
}
\description{
Extract several fits from a lasso fit table.}
\usage{
extract.lassogrp(object, i = NULL, lambda = NULL, fitfun = "lm", ...)

object[i]
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object with class \code{lassogrp}}
  \item{i}{a single index for \code{extract.lassogrp}, or a vector of
    indices for the simple subset operator '\code{[i]}'}
  \item{lambda}{alternatively to specifying \code{i}, a lambda value may
    be specified. }
  \item{fitfun}{Fitting function that determines the structure of the
    return value. Note that the coefficients and ...
    will be taken from \code{object}.}
  \item{\dots}{additional arguments passed to ...}
}
\details{
  \code{extract.lassogrp} generates an object of a regression class
  like \code{lm} or \code{regr}. This is useful for applying the
  respective plot and print methods to the lasso fit.
}
\value{
  \code{extract.lassogrp}: object of class 'lassofit' inheriting from
  the class specified by
  \code{fitfun}.

  \code{x[i]}: an object of class 'lassogrp' containing only the specified
  fits, i.e. all the information corresponding to these fits.
}
\author{Werner Stahel, \email{stahel@stat.math.ethz.ch}}
\examples{
data(asphalt)
rr <- lasso(log10(RUT)~log10(VISC)+ASPH+BASE+FINES+VOIDS+RUN,
  data=asphalt)
rr[c(1,10,15)]
extract.lassogrp(rr, lambda=2.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
