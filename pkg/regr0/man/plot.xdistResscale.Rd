\name{plot.xdistResscale}
\alias{plot.xdistResscale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Illustrate Completeness Test for Regression Model}
\description{
  Plots mean absolute differences of residuals against mean
  x-distance of pairs of observations, classified on the basis
  of their x-distance. The distance is proportional to the
  Mahalanobis distance in the design space.
}
\usage{
\method{plot}{xdistResscale}(x, lwd=2, cex=2, xlab = "distance in x space",
  ylab = "average abs. residual difference", col.aux = "grey30", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{xdistResscale}}
  \item{lwd, cex}{line width and character expansion for plotting}
  \item{xlab, ylab}{label for axes}
  \item{col.aux}{color to be used for auxiliary lines}
  \item{...}{arguments passed to \code{plot.default}}
  }
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{none.
}
\references{See package vignette.
}
\author{Werner A. Stahel, ETH Zurich}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{xdistResdiff}}
  }

\examples{
data(d.blast)
rr <- regr(tremor~distance+charge, data=d.blast)
## an inadequate model!
xdrd <- xdistResdiff(rr)
plot(xdrd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{hplot}% __ONLY ONE__ keyword per line
