\name{smoothM}
\alias{smoothM}
%% \alias{smoothMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smooth: wrapper function
}
\description{Generate fits of a smoothing function for multiple y's.
%%-  (\code{smoothM}) and multiple x's (\code{smoothMM}).
  Smooths can be calculated within given groups.
}
\usage{
  smoothM(x, y, weights = NULL, band = FALSE, group = NULL, power = 1,
          resid = "difference",
          par = 5 * length(x)^log10(1/2), parband = par*2^log10(2),
          iterations = 50, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of x values.
}
  \item{y}{vector or matrix of y values.
}
  \item{weights}{vector of weights to be used for fitting.
}
  \item{band}{logical: Should a band consisting of low and high smooth
    be calculated?
}
  \item{group}{NULL or a factor that defines the groups for which the
    smooth is calculated.
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
  \item{par, parband}{argument to be passed to the smoothing function,
    \code{parband} when calculating "band" smooths
}
  \item{iterations}{argument passed on to the smoothing function.
  }
  \item{...}{Further arguments, passed to the smoothing function.
  }
}
\details{These functions are useful for generating the smooths enhancing
  residual plots. (\code{smoothM}) generates a smooth for a single x
  variable and multiple y's. It is used to draw smooths from simulated
  residuals.
  (\code{smoothMM}) is the wrapper for multivariate regression.
  If argument \code{group} is specified, the smooths will be calculated
  within the specified groups.

  NA's in either \code{x} or any column of \code{y} cause dropping the
  observation (equivalent to \code{na.omit}).
}
\value{
  \code{smoothM}: A list with components:  
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
%%-   For \code{smoothMM}, the value is a list of such lists, with a component
%%-   for each column of \code{x}. 
}
\details{
  The smoothing function used to produce the smooth is
  \code{\link{smoothRegr}}, which relies \code{\link{loess}} by default.
  This may be changed via \code{options(smoothFunction = func)} where
  \code{func} is a smoothing function with the same arguments as
  \code{\link{smoothRegr}}.
}

\author{Werner A. Stahel, ETH Zurich}

\seealso{\code{\link{smoothRegr}}
}
\examples{
data(d.blast)
r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
r.smooth <- smoothM( fitted(r.blast), residuals(r.blast))
showd(r.smooth$y)
plot(fitted(r.blast), resid(r.blast), main="Tukey-Anscombe Plot")
abline(h=0)
lines(r.smooth$x,r.smooth$y, col="red")

## grouped data
r.smx <- smoothM( d.blast$dist, residuals(r.blast), group=d.blast$location)
plot(d.blast$dist, resid(r.blast), main="Residuals against Regressor")
abline(h=0)
for (lg in 1:length(levels(r.smx$group))) {
  li <- as.numeric(r.smx$group)==lg 
  lines(r.smx$x[li],r.smx$y[li], col=lg+1, lwd=3)

## example with data from stats
data(swiss)
  plot(Fertility~Agriculture, swiss)

  r.sm <- smoothM( swiss$Agriculture, swiss$Fertility )
  lines(r.sm)

  r.sm <- smoothM( swiss$Agriculture, swiss$Fertility, band=TRUE )
  li <- r.sm$ybandind
  lines(r.sm$x[li], r.sm$yband[li], lty=2)
  lines(r.sm$x[!li], r.sm$yband[!li], lty=2)

## example with groups and band
  t.group <- swiss$Catholic>70
  plot(Fertility~Agriculture, swiss, pch=2+t.group, col=2+t.group)
  r.sm <- smoothM( swiss$Agriculture, swiss$Fertility,
                   group=t.group, band=TRUE )

  for (lg in c(FALSE,TRUE)) {
    lig <- which(r.sm$group==lg)
    lines(r.sm$x[lig], r.sm$y[lig], lty=1, col=2+lg)
    li <- r.sm$ybandind[lig]
    ligh <- lig[li]
    ligl <- lig[!li]
    lines(r.sm$x[ligh], r.sm$yband[ligh], lty=2+2*lg, col=2+lg)
    lines(r.sm$x[ligl], r.sm$yband[ligl], lty=2+2*lg, col=2+lg)
  }

}

## use of simresiduals is missing

## multivariate regression
data(d.fossiles)
r.mregr <-
  regr(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+lChlorophyll+region+N,
                data=d.fossiles)
r.msmooth <- smoothMM( fitted(r.mregr), residuals(r.mregr), resid=TRUE )
showd(r.msmooth[[2]]$resid)
r.msmresid <- sapply(r.msmooth, function(x) x$resid)
  ## row.names are lost by sapply
showd(r.msmresid)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}

