\name{plot.lassoselect}
\alias{plot.lassoselect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a Lassoselect Object}
\description{
  Plots the coefficients as a function of the tuning parameter "s" of
  the L1 penalty term in the lasso model selection method and the cross
  validated mean square prediction error.
}
\usage{
\method{plot}{lassoselect}(x, plot.s = TRUE, plot.cv = TRUE,
  plot.cvse = TRUE, main = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the result of \code{\link{lassoselect}}}
  \item{plot.s}{if TRUE, the coefficients will be plotted}
  \item{plot.cv}{if TRUE, the cross validated mean square prediction
    error will be plotted} 
  \item{plot.cvse}{if TRUE, error bars for these values will be plotted}
  \item{main}{main title of plot}
  \item{\ldots}{further arguments, passed to \code{plot.lars}}
}
%%- \details{
%%-   ~~ If necessary, more details than the description above ~~
%%- }
\value{
  none.
}
% \references{ ~put references to the literature/web site here ~ }
\author{Werner A. Stahel, ETH Zurich}
\note{\code{plot.lars} has been slightly modified
}
\seealso{\code{\link[lars]{plot.lars}}}
\examples{
## data(d.blast)
## t.r <- lassoselect(log10(tremor)~location+log10(distance)+log10(charge),
##                    data=d.blast)
## plot(t.r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{regression}% __ONLY ONE__ keyword per line
