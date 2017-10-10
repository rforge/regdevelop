\name{print.regr}
\alias{print.regr}
\title{Print Method for 'regr' objects}
\description{
  Printing \code{"regr"} objects, i.e. the results of
  \code{\link{regr}(*)}.
}
\usage{
\method{print}{regr}(x, call=TRUE, correlation = FALSE,
      dummycoef = getUserOption("show.dummycoef"),
      termcolumns = getUserOption("termcolumns"),
      allcoefcolumns = getUserOption("allcoefcolumns"),
      digits = max(3, getUserOption("digits") - 2), symbolic.cor = p > 4,
      signif.stars = getOption("show.signif.stars"),
      na.print = getUserOption("na.print"),
      residuals = FALSE, niterations = FALSE, ...)
}
\arguments{
  \item{x}{'regr' object}
  \item{call}{logical: should the call to \code{regr} be shown?}
  \item{correlation}{logical: should correlation among estimated
    coefficients be shown?}
  \item{dummycoef}{logical: should the dummy coefficients for the
    factors be shown? By default, coefficients corresponding to factors
    are shown}
  \item{termcolumns}{columns of the term table to be printed.
    If \code{NULL}, all columns will be printed.}
  \item{allcoefcolumns}{columns of the (dummy) coefficent tables to be
    printed.} % See \code{\link{print.allcoef}}.
  \item{digits}{number of significant digits to be shown}
  \item{symbolic.cor}{logical: should correlations be shown by symbols?}
  \item{signif.stars}{logical: should significance asterisks be shown?
    A legend will also be shown after the term table.
    To avoid the legend, set \code{signif.stars} false and include
    \code{"p.symb"} in \code{allcoefcolumns}.}
  \item{na.print}{character string used to display NAs}
  \item{residuals}{logical: should summary(residuals) be shown?}
  \item{niterations}{???}
  \item{\dots }{further arguments, passed to the print methods used.}
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
  print(r.blast,
    termcolumns = c("coef","ciLow","ciHigh","df","signif"),
    allcoefcolumns = c("coef","ciLow","ciHigh"),
    correlation=TRUE, symbolic.cor=TRUE)
}
\keyword{print}
