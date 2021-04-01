\name{colorpale}
\alias{colorpale}

\title{determine more pale colors for given colors
}
\description{
Finds colors that are `equivalent' to the colors given as the first
argument, but more pale or less pale
}
\usage{
colorpale(col = NA, pale = NULL, rgb = FALSE, ...)
}

\arguments{
  \item{col}{
    a color or a vector of colors for which the pale version should be found
  }
  \item{pale}{
    number between -1 and 1 determining how much paler the result should
    be. If \code{=0}, the original color, \code{col} will be returned
    unchanged (but in the 'rgb' or 'hexadecimal' form). If \code{=1} or
    \code{-1},     the result is white (\code{#FFFFFF}) or black, respectively. 
  }
  \item{rgb}{should result be expressed in 'rgb' form? If \code{FALSE},
    it will be in hexadecimal form.
  }
  \item{...}{
    further arguments passed on to \code{\link{rgb}}
    if \code{rgb} is \code{FALSE}
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

