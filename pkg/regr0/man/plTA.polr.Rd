\name{plTA.polr}
\alias{plTA.polr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Residuals against Fitted Values for polr objects}
\description{
  Displays the "Tukey-Anscombe" plot for ordered regression models.
}
\usage{
plTA.polr(object, colbars = grey(0.7), colref = grey(0.7), ploty = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{result of a call to \code{polr}}
  \item{colbars}{color to be used for plotting the bars representing
    the residuals}
  \item{colref}{color for reference line}
  \item{ploty}{if TRUE, the latent response will be plotted instead of
    the residuals}
}
\details{
  The plot intends to show the residuals for the latent variable Z that
  determines the response values Y. Since the precise residuals of the
  latent variable cannot be known, the plot shows their conditional
  distribution, given the observed response values. This conditional
  distribution is characerized by its median (by -) and the two other
  quartiles (by vertical bars).
}
\value{
  none.
}
\references{See \code{http://stat.ethz.ch/~stahel/regression}}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{plot.regr}}, \code{\link{residuals.polr}} }
\examples{
house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
plTA.polr(house.plr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{regression}% __ONLY ONE__ keyword per line
