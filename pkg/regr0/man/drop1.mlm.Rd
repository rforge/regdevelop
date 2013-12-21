\name{drop1.mlm}
\alias{drop1.mlm}
\alias{add1.mlm}
\title{Drop All Possible Single Terms from or Add Terms to a Multivariate Model}
\description{
  \code{drop1.mlm} / \code{add1.mlm} performs a test for the suitability
  of dropping / adding each term in scope from / to a multivariate
  regression or manova model
}
\usage{
\method{drop1}{mlm}(object, scope = NULL,
  test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"), total = TRUE,
    add=FALSE, ...)
\method{add1}{mlm}(object, scope = NULL,
  test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a multivariate regression or manova object, possibly
    with class \code{regr}}
  \item{scope}{as in general \code{drop1 / add1} }
  \item{test}{multivariate test to be used}
%%  \item{intercept}{as in general \code{drop1 / add1} }
\item{total}{should the test for null model also be performed?}
\item{add}{if TRUE, the function performs \code{add1} instead of
  dropping terms}
\item{\dots}{not used} % additional arguments passed to \code{drop1}
}
%- \details{
%-   ~~ If necessary, more details than the description above ~~
%- }
\value{
  matrix containing statistics for each term
}
%- \references{ ~put references to the literature/web site here ~ }
\author{modification of summary.manova by W. Stahel}
\note{Note that \code{regr} performs such a (Wilks) test automatically.
Bug: extract.AIC is not yet available
}
\seealso{\code{\link{drop1.regr}}, \code{\link{drop1.lm}}}
\examples{
data(d.fossiles)
r.lm <-
  lm(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+lChlorophyll+region+N,
                data=d.fossiles)
drop1(r.lm)  ##  should be the same as summary(manova(...)) for its last row.
r.mregr <-
  regr(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+lChlorophyll+region+N,
                data=d.fossiles)
r.mregr$drop1
drop1(r.mregr, test="Pillai")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%- \keyword{}
%- \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
