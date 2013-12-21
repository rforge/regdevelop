\name{legendr}
\alias{legendr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add a Legend to a Plot}
\description{
  Adds a legend to a plot as does \code{\link{legend}}. This function
  just expresses the position relative to the range of the coordinates
}
\usage{
legendr(x = 0.05, y = 0.95, legend, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{position in horizontal direction, between 0 for left margin
    and 1 for right margin}
  \item{y}{position in vertical direction, between 0 for bottom margin
    and 1 for top margin}
  \item{legend}{text of the legend}
  \item{\dots}{arguments passed to \code{\link{legend}}}
}
\value{
  See \code{\link{legend}}
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{legend}}
  }
  \examples{
     ts.plot(ldeaths, mdeaths, fdeaths,xlab="year", ylab="deaths", lty=c(1:3))
     legendr(0.7,0.95, c("total","female","male"), lty=1:3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{aplot}
