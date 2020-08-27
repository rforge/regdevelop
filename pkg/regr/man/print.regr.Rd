\name{print.regr}
\alias{print.regr}
\title{Print Method for 'regr' objects}
\description{
  Printing \code{"regr"} objects, i.e. the results of
  \code{\link{regr}(*)}.
}
\usage{
\method{print}{regr}(x, call=TRUE, residuals = FALSE, 
  termeffects = TRUE, coefcorr = FALSE, niterations = FALSE,
  printstyle = NULL, digits = getRegrOption("digits"),
  na.print = NULL, symbolic.cor = p > 4, ...)
}
\arguments{
  \item{x}{'regr' object}
  \item{call}{logical: should the call to \code{regr} be shown?}
  \item{residuals}{logical: should summary(residuals) be shown?}
  \item{termeffects}{logical: should the dummy coefficients for the
    factors be shown? if TRUE (default), coefficients corresponding
    to factors are shown. If $>1$, all termeffects are shown}
%%   \item{termcolumns}{columns of the term table to be printed.
%%     If \code{NULL}, all columns will be printed.}
%%   \item{termeffcolumns}{columns of the term effects tables to be
%%     printed.} % See \code{\link{print.termeffects}}.
%%   \item{coefcolumns}{columns of the (dummy) coefficent table to be
%%     printed. This table is printed only if there is no termtable}
  % See \code{\link{print.coeftable}}.
  \item{coefcorr}{logical: should correlation among estimated
    coefficients be shown?}
  \item{niterations}{logical: should number of iterations be shown
    (if available)?}
  \item{printstyle}{determines the columns to be printed}
  \item{digits}{number of significant digits to be shown}
  \item{na.print}{character string used to display NAs}
  \item{symbolic.cor}{logical: should correlations be shown by symbols?}
  \item{\dots }{further arguments, passed to the print methods used.}
}
\details{
  The \code{tit} attribute of \code{x} will be printed if available and
  \code{getRegrOption("doc") > 0}, and any \code{doc} attribute,
  if \code{getRegrOption("doc") >= 2} (see \code{\link{tit}}).
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
    termcolumns = c("coef","ciLow","ciUp","df","signif0"),
    termeffbcolumns = c("coef","ciLow","ciUp"),
    correlation=TRUE, symbolic.cor=TRUE)
}
\keyword{print}
