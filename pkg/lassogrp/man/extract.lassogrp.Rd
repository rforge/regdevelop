\name{extract.lassogrp}
\title{Extract Regression Results from a Lasso Fit}
\Rdversion{1.1}
\alias{extract.lassogrp}
\alias{[}% <- needed for R CMD check {MM thinks this is bogous}
\alias{[.lassogrp}
\description{
  Extract a regression fits from a lasso fit (table).
}
\usage{
extract.lassogrp(object, i = NULL, lambda = NULL, data=NULL, fitfun = "lm", ...)

\S3method{[}{lassogrp}(x, i)
%] << for emacs
}
\arguments{
  \item{object, x}{an object of class \code{"lassogrp"}.}
  \item{i}{a single index for \code{extract.lassogrp}, or a vector of
    indices for the simple subset operator '\code{[i]}'}
  \item{lambda}{alternatively to specifying \code{i}, a single (numeric)
    lambda value may be specified.}
  \item{data}{the data originally used which must still be available.
    (The latter restriction will possibly be relaxed in the future.).
    By default it is found in the environment.}
  \item{fitfun}{fitting function that determines the structure of the
    return value.  Note that coefficients and more will be taken from
    \code{object}.}
  \item{\dots}{additional arguments passed to \code{fitfun}.}
}
\details{
  \code{extract.lassogrp} generates an object of a regression class
  like \code{lm} or \code{regr}.  This is useful for applying the
  respective plot and print methods to the lasso fit.
  The result of an unpenalized fit to the ``selected'' model
  (terms with non-zero coefficients) is available as 'fit.unpen'
  component of the result.
}
\value{
  \code{extract.lassogrp}: object of class \code{lassofit} inheriting from
  the class specified by \code{fitfun}.
  \code{$fit.unpen}: The result of fitting the model (by code{fitfun})
  to the reduced model.

  \code{x[i]}: an object of class \code{lassogrp} containing only the specified
  fits, i.e. all the information corresponding to these fits.
}
\author{Werner Stahel, \email{stahel@stat.math.ethz.ch}}
\examples{
data(asphalt)
rr <- lasso(log10(RUT) ~ log10(VISC) + ASPH+BASE+FINES+VOIDS+RUN,
            data=asphalt, adaptive=TRUE)
## Extract results for three lambda's:
rr[c(1,19,20)]
extract.lassogrp(rr, 19)

## The above relies on finding the original data;
## it does not work otherwise
d.a <- asphalt
rm(asphalt)
try(extract.lassogrp(rr, lambda=2.5)) # -> error: cannot find 'asphalt'
%% MM thinks this is a design bug! (FIXME!)
## it works if you can specify the data :
extract.lassogrp(rr, lambda=2.5, data=d.a)
}
\keyword{utilities}
