\name{print.regr}
\alias{print.regr}
\title{Print Method for 'regr' objects}
\description{
  Printing \code{"regr"} objects, i.e. the results of
  \code{\link{regr}(*)}.
}
\usage{
\method{print}{regr}(x, call=TRUE, correlation = FALSE,
      termeffects = getUseroption("show.termeffects"),
      termcolumns = getUseroption("termcolumns"),
      termeffcolumns = getUseroption("termeffcolumns"),
      digits = getUseroption("digits"), symbolic.cor = p > 4,
      signif.stars = getOption("show.signif.stars"),
      na.print = getUseroption("na.print"),
      residuals = FALSE, niterations = FALSE, ...)
}
\arguments{
  \item{x}{'regr' object}
  \item{call}{logical: should the call to \code{regr} be shown?}
  \item{correlation}{logical: should correlation among estimated
    coefficients be shown?}
  \item{termeffects}{logical: should the dummy coefficients for the
    factors be shown? if TRUE (default), coefficients corresponding
    to factors are shown. If $>1$, all termeffects are shown}
  \item{termcolumns}{columns of the term table to be printed.
    If \code{NULL}, all columns will be printed.}
  \item{termeffcolumns}{columns of the (dummy) coefficent tables to be
    printed.} % See \code{\link{print.termeffects}}.
  \item{digits}{number of significant digits to be shown}
  \item{symbolic.cor}{logical: should correlations be shown by symbols?}
  \item{signif.stars}{logical: should significance asterisks be shown?
    A legend will also be shown after the term table.
    To avoid the legend, set \code{signif.stars} false and include
    \code{"p.symb"} in \code{termeffcolumns}.}
  \item{na.print}{character string used to display NAs}
  \item{residuals}{logical: should summary(residuals) be shown?}
  \item{niterations}{logical: should number of iterations be shown
    (if available)?}
  \item{\dots }{further arguments, passed to the print methods used.}
}
\details{
  The \code{tit} attribute of \code{x} will be printed if available and
  \code{getUseroption("doc") > 0}, and any \code{doc} attribute,
  if \code{getUseroption("doc") >= 2} (see \code{\link{tit}}).
}
\value{Print functions return the object to be printed invisibly.
}
\author{Werner A. Stahel, ETH Zurich}
\examples{
  data(d.blast, package="plgraphics")
  r.blast <-
       regr(log10(tremor)~location+log10(distance)+log10(charge),
            data=d.blast)
  print(r.blast,
    termcolumns = c("coef","ciLow","ciHigh","df","signif"),
    termeffbcolumns = c("coef","ciLow","ciHigh"),
    correlation=TRUE, symbolic.cor=TRUE)
}
\keyword{print}
