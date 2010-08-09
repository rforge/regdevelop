\name{print.lassogrp}
\Rdversion{1.1}
\alias{print.lassogrp}
\alias{print.lassofit}
\title{Print lassogrp or lassofit Objects}
\description{
  Methods for \code{\link{print}}ing of \code{"lassogrp"} or
  \code{"lassofit"} objects.
}

\usage{
\method{print}{lassogrp}(x, coefficients = TRUE, doc = options("doc")[[1]], ...)
\method{print}{lassofit}(x, residuals = FALSE, doc = options("doc")[[1]], ...)
}
\arguments{
  \item{x}{Object to be printed}
  \item{coefficients}{logical: should coefficients be shown?}
  \item{residuals}{logical: should a summary of residuals be shown?}
  \item{doc}{integer specifying if documentation of \code{x} should be
    printed (if any):
    \describe{
      \item{0:}{do not show documentation}
      \item{1:}{show \code{attr(x,'tit')}}
      \item{2:}{show \code{attr(x,'descr')} in addition}
    }
  }
  \item{...}{further arguments, passed to \code{\link{print}}.}
}
\value{\code{x}, invisibly.}% requirement of *every* print() method in R

\author{Werner A. Stahel, ETH Zurich, \email{stahel@stat.math.ethz.ch}}

\examples{
data(asphalt)
rr <- lasso(log10(RUT)~log10(VISC)+ASPH+BASE+FINES+VOIDS+RUN,
            data= asphalt)
print(rr, coefficients=FALSE)

rrf <- extract.lassogrp(rr, 13)
print(rrf, residuals=TRUE)
}
\keyword{misc}
