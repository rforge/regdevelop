\name{stdresiduals}
\alias{stdresiduals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get Standardized Residuals}
\description{
  Calculates standardized residuals and leverage values
}
\usage{
stdresiduals(x, residuals=NULL, sigma=x$sigma, weights=NULL,
    leveragelimit = c(0.99, 0.5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    an object containing a component \code{fit$leverage} and possibly
    a component \code{fit$na.action}
  }
  \item{residuals}{unstandardized residuals}
  \item{sigma}{error standard deviation or other scale}
  \item{weights}{weights}
  \item{leveragelimit}{limit on leverage values}
}
\value{
  vector or matrix of standardized residuals,
  with attributes\cr
  \code{attr(.,"stdresratio")}:
  ratio of standardized / unstandardized residuals,\cr
  \code{attr(.,"leverage")}: leverage (hat) values,\cr
  \code{attr(.,"weighted")}: weights used in the standardization,\cr
  \code{attr(.,"stddev")}: error standard deviation or scale parameter.
}
\details{
  The difference to \code{stdres} of package \code{MASS} 
  is that \code{stdresiduals} also applies to multivariate regression
  and can be used with other regression model fits.

  The function uses the \code{qr} decomposition of \code{object}.
  If necessary, it generates it.
}

\author{Werner A. Stahel, ETH Zurich}

\seealso{\code{\link{stdres}; \link{hat}; \link{hatvalues}; \link{influence}}}

\examples{
data(d.blast)
r.blast <-
  lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
t.stdr <- stdresiduals(r.blast)
showd(t.stdr)
showd(attr(t.stdr, "leverage"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
% \keyword{diagnostics}
