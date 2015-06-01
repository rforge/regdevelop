\name{quinterpol}
\alias{quinterpol}

\title{Interpolated Quantiles
}
\description{This function implements a version of empirical quantiles
  based on interpolation 
}
\usage{
quinterpol(x, probs = c(0.25, 0.5, 0.75), extend = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of data determining the quantiles 
}
  \item{probs}{vector of probabilities defining which quantiles should
    be produced
}
  \item{extend}{logical: Should quantiled be calculated outside the
    range of the data by linear extrapolation?
}
}
\details{The empirical quantile function jumps at the data values
  according to the usual definition.
  The version of quantiles calculated by 'quinterpol' avoids jumps.
  It is based on linear interpolation of the step version of the
  empirical cumulative distribution function, using as the given points
  the midpoints of both vertical and horizontal pieces of the latter.
  See 'examples' for a visualization.
}
\value{vector of quantiles
}
\author{Werner A. Stahel
}

\seealso{
quantile
}
\examples{
## This example illustrates the definition of the "interpolated quantiles"

set.seed(2)
t.x <- sort(round(2*rchisq(20,2)))
table(t.x)
t.p <- ppoints(100)
plot(quinterpol(t.x,t.p),t.p, type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
