\name{cv.lasso}
\Rdversion{1.1}
\alias{cv.lasso}
\title{Cross Validation of a Lasso Fit}
\description{
Calculation of cross validation for a \code{lassogrp} object.
}
\usage{
cv.lasso(object, blocks = 10, trace = FALSE,
  control = lassoControl(trace = 0), env = globalenv(),
  plot.it = NULL, se = TRUE, ...)
}
\arguments{
  \item{object}{a \code{lassogrp} object}
  \item{blocks}{number of blocks to be used in block crossvalidation,
    or name of variable defining the blocks}
  \item{trace}{logical: Should completion of each cross validation cycle
    be indicated?}
  \item{control}{see \code{?lassoControl}}
  \item{env}{environment for getting blocks}
  \item{plot.it}{logical: should results be plotted?}
  \item{se}{logical: should cross validation standard errors be plotted?}
  \item{\dots}{further arguments passed to \code{plot}}
}
\details{If the blocks are determined by a factor which appears in the
  model, ...
}
\value{
\item{rmse}{}
\item{rmse.see}{}
\item{rmse.blocks}{}
\item{fitted}{}
\item{lambda}{}
\item{blocksinmodel}{}
}
\author{Werner Stahel, \email{stahel@stat.math.ethz.ch}}

\seealso{
\code{cv.lars} in package \code{lars}
}
\examples{
data(asphalt)
rr <- lasso(log10(RUT)~log10(VISC)+ASPH+BASE+FINES+VOIDS+RUN, data=asphalt)
cv.lasso(rr)
}
\keyword{regression}
