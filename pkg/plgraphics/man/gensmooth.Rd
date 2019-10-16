\name{gensmooth}
\alias{gensmooth}

\title{Smooth: wrapper function
}
\description{Generate fits of a smoothing function for multiple y's.
  Smooths can be calculated within given groups.
}
\usage{
gensmooth(x, y, band = FALSE, power = 1, resid = "difference",
  plargs=NULL, ploptions=NULL, ...)
}
\arguments{
  \item{x}{vector of x values.
  }
  \item{y}{vector or matrix of y values.
  }
  \item{band}{logical: Should a band consisting of low and high smooth
    be calculated? It will only be calculated for the first column of \code{y}.
  }
  \item{power}{\code{y} will be raised to \code{power} before smoothing.
    Results will be back-transformed. (Useful for smoothing absolute
    values for a 'scale plot', for which \code{power=0.5} is recommended.)    
  }
  \item{resid}{Which residuals be calculated?
    \code{resid=1} or \code{="difference"} means usual residuals;
    \code{resid=2} or \code{="ratio"} means $y_i/\\hat y_i$, which is useful
    to get scaled y's (regression residuals) according to a smooth fit in the
    scale plot.
  }
  \item{plargs, ploptions}{result of calling \code{pl.control}.
    The component \code{plargs$pdata} may contain \code{smooth.weight}
    and \code{smooth.group}, and \code{ploptions} specifies
    \code{smoothPar} and \code{smoothIter}.
    All of these may be used by the smoothing function.
  }
  \item{...}{Further arguments, passed to the smoothing function.
  }
}
\details{This function is useful for generating the smooths enhancing
  residual plots.
  It generates a smooth for a single x
  variable and multiple y's.
  It is also used to draw smooths from simulated residuals.

  NA's in either \code{x} or any column of \code{y} cause dropping the
  observation (equivalent to \code{na.omit}).

  The smoothing function used to produce the smooth is
  \code{\link{smoothRegr}}, which relies \code{\link{loess}}, by default.
  This may be changed via \code{ploptions(smooth.function = func)} where
  \code{func} is a smoothing function with the same arguments as
  \code{\link{smoothRegr}}.

  The result of the smoothing function may carry an attribute
  \code{xtrim}. This regulates if the fitted values corresponding to
  extreme x values will be suppressed when plotting:
  The number of extreme x values corresponding to
  \code{ploptions("smooth.xtrim")} will be multiplied by
  this attribute to obtain the number of extreme points suppressed at
  each end. If the smoothing function is \code{smoothLm}, which fits a
  straight line, then trimming is suppressed since this function returns
  0 as the \code{xtrim} attribute.

  If \code{band} is \code{TRUE}, a vector of "low" and a vector of
  "high" smooth values will be calculated for the first column of
  \code{y} in the following way:
  Residuals are calculated as the diference
  between the observations and the respective smoothed values hat.$s_i$.
  Then a smooth is calculated for the square roots of the positive residuals,
  and the squared fitted values are added to the hat.$s_i$.
  (The transformation by square roots makes the distribution of the residuals 
  more symmetric.)
  This defines the ``high'' smooth values. 
  The construction of the ``low'' one is analogous.
  The resulting values of the two are stored in the list component
  \code{yband}, and \code{ybandindex} contains the information to which
  group ("low" or "high") the value belongs.
}
\value{
  A list with components:  
  \item{x}{vector of x values, sorted, within levels of \code{group}
    if grouping is actif.}
  \item{y}{matrix with 1 or more columns of corresponding fitted values
    of the smoothing.}
  \item{group}{grouping factor, sorted, if actif. \code{NULL} otherwise.}
  \item{index}{vector of indices of the argument \code{x} used for
    sorting. This is useful to relate the results to the input.
    Use \code{ysmoothed[value$index,] <- value$y} to get values
    corresponding to input \code{y}.}
  \item{xorig}{original \code{x} values}
  \item{ysmorig}{corresponding fitted values}
  \item{residuals}{if required by the argument \code{resid},
    residuals from the smooth fit are provided in the original order,
    i.e. \code{value$resid[i,j]} corresponds to the input
    \code{value$y[i,j]}.}
  If \code{band==TRUE},
  \item{yband}{vector of low and high smoothed values (for the first
    column of \code{y})}
  \item{ybandindex}{Indicator if \code{yband} is a high value}
}
\note{
  This function is called by \code{\link{plyx}} and
  \code{\link{plmatrix}} when \code{smooth=T} is set,
  as well as by
  \code{\link{plregr}} applied to model objects.
  It is rarely needed to call it directly.\cr
  A band is generated only for the first columnn of \code{y} since the
  others are supposed to be simulated versions of the first one
  and do not need a band. 
}

\author{Werner A. Stahel, ETH Zurich}

\seealso{\code{\link{smoothRegr}},
  \code{\link{plsmooth}}, \code{\link{plsmoothline}}
}
\examples{
data(d.blast)
r.blast <-
  lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast,
    na.action=na.exclude)
r.smooth <- gensmooth( fitted(r.blast), residuals(r.blast))
showd(r.smooth$y)
plot(fitted(r.blast), resid(r.blast), main="Tukey-Anscombe Plot")
abline(h=0)
lines(r.smooth$x,r.smooth$y, col="red")

## grouped data
t.plargs <- list(pdata=data.frame(".smooth.group."=d.blast$location))

r.smx <- gensmooth( d.blast$dist, residuals(r.blast), plargs=t.plargs)

plot(d.blast$dist, residuals(r.blast), main="Residuals against Regressor")
abline(h=0)
plsmoothline(r.smx, d.blast$dist, resid(r.blast), plargs=t.plargs)
## or, without using plsmoothlines:
## for (lg in 1:length(levels(r.smx$group))) {
##   li <- as.numeric(r.smx$group)==lg 
##   lines(r.smx$x[li],r.smx$y[li], col=lg+1, lwd=3)
## }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}

