\name{drop1.mlm}
\title{Drop All Possible Single Terms from or Add Terms to a Multivariate Model}
\alias{drop1.mlm}
\alias{add1.mlm}
\description{
  \code{drop1.mlm} / \code{add1.mlm} performs a test for the suitability
  of dropping / adding each term in scope from / to a multivariate
  regression or manova model.
}
\usage{
\method{drop1}{mlm}(object, scope = NULL,
      test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"),
      total = TRUE, add = FALSE, ...)
\method{add1}{mlm}(object, scope = NULL,
     test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"), ...)
}
\arguments{
  \item{object}{a multivariate regression or manova object, possibly
    with class \code{"regr"}, typically resulting from \code{\link{regr}()}.}
  \item{scope}{see \code{\link{drop1}} (in \R's package \pkg{stats}).}
  \item{test}{a \code{\link{character}} string specifying the
    multivariate test to be used.}
  \item{total}{logical specifying if the test for null model should also
    be performed (and a \code{"<total>"} row be added to the resulting
    table).}
  \item{add}{logical, not to be set typically;  \code{add=TRUE}
    is basically equivalent to \code{add1.mlm()}.}
  \item{\dots}{potentially further arguments passed from or to methods;
    for unused for \code{drop1()}.  For \code{add1()}: Arguments passed
    to \code{\link{drop1}(..)}.}
}
%% \details{ %% FIXME: explain the four different tests
%% }
%% \references{
%% }
\value{
  a \code{\link{data.frame}} matrix containing statistics for each term.
}
\author{modification of \R's \code{\link{summary.manova}} by
  Werner Stahel
}
\note{Note that \code{regr} performs such a (Wilks) test automatically.
  Bug: extract.AIC is not yet available
}
\seealso{\code{\link{drop1.regr}}, \code{\link{drop1.lm}}}
\examples{
data(d.fossileSamples, package="plgraphics")
r.lm <-
  lm(cbind(sAngle,lLength,rWidth) ~ SST+Salinity+lChlorophyll+Region+N,
     data=d.fossileSamples)
drop1(r.lm)  ##  should be the same as summary(manova(...)) for its last row.
r.mregr <-
  regr(cbind(sAngle,lLength,rWidth) ~ SST+Salinity+lChlorophyll+Region+N,
       data=d.fossileSamples)
r.mregr$drop1
drop1(r.mregr, test="Pillai")
}
\keyword{models}
