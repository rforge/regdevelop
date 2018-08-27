\name{colorpale}
\alias{colorpale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{determine more pale colors for given colors
%%  ~~function to do ... ~~
}
\description{
Finds colors that are `equivalent' to the colors given as the first
argument, but more pale.
}
\usage{
colorpale(col = NA, pale = 0.3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{col}{
  a color or a vector of colors for which the pale version should be found
}
  \item{pale}{
  number between 0 and 1 determining how much paler the result should be.
}
  \item{...}{
  further arguments passed on to \code{\link{rgb}}
}
}
\details{
  The function increases rgb coordinates of colors `proportionally':
    crgb <- t(col2rgb(col)/255);
    rgb(1 - pale * (1 - crgb))
}
\value{character vector: names of colors to be used as color argument
  for graphical functions.
}
\author{Werner A. Stahel, ETH Zurich
}

\seealso{
\code{\link{rgb}}
}
\examples{
( t.col <- colorpale(c("red","blue")) )
plot(0:6, type="h", col=c("black","red","blue",t.col, colorpale(t.col)), lwd=5)
}
\keyword{utility}

