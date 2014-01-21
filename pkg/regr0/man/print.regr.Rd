\name{print.regr}
\alias{print.regr}
\title{Print Method for 'regr' objects}
\description{
  Printing \code{"regr"} objects, i.e. the results of
  \code{\link{regr}(*)}.
}
\usage{
\method{print}{regr}(x, correlation = FALSE, dummycoef = NULL,
      digits = max(3, getUserOption("digits") - 2), symbolic.cor = p > 4,
      signif.stars = getOption("show.signif.stars"),
      residuals = FALSE, niterations = FALSE, ...)
}
\arguments{
  \item{x}{'regr' object}
  \item{correlation}{logical: should correlation among estimated
    coefficients be shown?}
  \item{dummycoef}{logical: should the dummy coefficients for the
    factors be shown?}
  \item{digits}{number of significant digits to be shown}
  \item{symbolic.cor}{logical: should correlations be shown by symbols?}
  \item{signif.stars}{logical: should significance asterisks be shown?
    If >1, a legend will also be shown.}
  \item{residuals}{logical: should summary(residuals) be shown?}
  \item{niterations}{???}
  \item{\dots}{further arguments, passed to the print methods used.}
}
\details{
  The \code{tit} attribute of \code{x} will be printed if available and
  \code{getUserOption("doc") > 0}, and any \code{doc} attribute,
  if \code{getUserOption("doc") >= 2} (see \code{\link{tit}}).
}
\value{Print functions return the object to be printed invisibly.
}
\author{Werner A. Stahel, ETH Zurich}
\examples{
  data(d.blast)
  r.blast <-
       regr(log10(tremor)~location+log10(distance)+log10(charge),
            data=d.blast)
  print(r.blast, correlation=TRUE, symbolic.cor=TRUE, dummycoef=FALSE)
}
\keyword{print}
