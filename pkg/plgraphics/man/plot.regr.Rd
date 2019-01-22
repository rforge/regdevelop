\name{plot.regr}
\alias{plot.regr}
\alias{plresx}
\title{Diagnostic Plots for Regr Objects}
\description{
  Diagnostic plots for fitted regression models:
  Residuals versus fit (Tukey-Anscombe plot) and/or target variable
  versus fit;
  Absolute residuals versus fit to assess equality of error variances;
  Normal Q-Q plot (for ordinary regression models);
  Residuals versus leverages to identify influential observations;
  Residuals versus sequence (if requested);
  and residuals versus explanatory variables.
  These plots are adjusted to the type of regression model.
}
\usage{
\method{plot}{regr}(x, data = NULL, plotselect = NULL, xvar = TRUE,
  transformed = FALSE, sequence = FALSE, weights = NULL,
  addcomp = FALSE, smooth.legend = FALSE, markextremes = NA, ...)

plresx(x, data = NULL, xvar = TRUE, transformed=FALSE,
  sequence = FALSE, weights = NULL,
  addcomp = FALSE, smooth.legend = FALSE, markextremes = NA,
  plargs = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{regr %(or lm or glm)
    object (result of a call to \code{regr}). Thiw iw the only argument
    that is needed. All others have useful defaults.}
  \item{data}{data set where explanatory variables and the following
    possible arguments are found: \code{weights, plweights, pch, plabs}}
  \item{plotselect}{which plots should be shown? See Details}
  \item{xvar}{if TRUE, residuals will be plotted versus all 
    explanatory variables (or terms, according to argument 'transformed')
    in the model (\code{plot.regr} will call \code{plresx}).\cr
    If it is a character vector, it contains the variables to be used.\cr
    If it is a formula, its right hand side contains these variables.
    The model formula is updated by such a formula.
    Whence, the use of \code{\~{}.+} adds variables to those in the
    model.\cr    
    If any variables are not be contained in the model, the argument
    \code{data} is needed.
  }
%%  \item{formula}{same as \code{xvar}. If both arguments are given,
%%    \code{formula} is used.}
  \item{transformed}{logical: should residuals be shown against
    transformed explanatory variables? If \code{TRUE}, the variables are
    transformed as implied by the model.}
  \item{sequence}{if TRUE, residuals will be plotted versus the
    sequence as they appear in the data.
    If another explanatory variable is monotone increasing or
    decreasing, the plot is not shown, but a warning is given.}
  \item{weights}{if TRUE, residuals will be plotted versus
    \code{x$weights}. Alternatively, a vector of weights can be specified}
  \item{addcomp}{ logical: should component effects be added to
    residuals for residuals versus input variables plots? }
  \item{smooth.legend}{When a grouping factor is used
    (argument \code{smooth.group}, see below),
    this argument determines whether and where the legend
    for identifying the groups should be shown, see Details}
  \item{markextremes}{proportion of extreme residuals to be labeled.
    If all points should be labeled, let \code{markextremes=1}.}
  \item{plargs}{result of calling \code{pl.control}.
    If \code{NULL}, \code{pl.control} will be called to generate it.
  }
  \item{\dots}{Many further arguments are available to customize the plots,
    see below for some of the most useful ones, and
    \code{\link{plotregr.control}} for a complete list.}
}
\details{
  Argument \code{plotselect} is used to determine which plots will be
  shown. It should be a named vector of numbers indicating
  \describe{
    \item{0}{do not show}
    \item{1}{show without smooth}
    \item{2}{show with smooth (not for \code{qq} nor \code{leverage})}
    \item{3}{show with smooth and smooth band (only for \code{resfit}
      in \code{plot.regr} and in \code{plresx})}      
  }
  The default is
  \code{c( yfit=0, resfit=smdef, absresfit = NA, absresweights = NA, qq = NA,
    leverage = 2, resmatrix = 1, qqmult = 3)}, where
  \code{smdef} is 3 (actually argument \code{smooth} of
  \code{plotregr.control} plus 1) for normal random deviations and
  one less (no band) for others. 

  Modify this vector to change the selection and the sequence in
  which the plots appear.
  Alternatively, provide a named vector defining all plots that should
  be shown on a different level than the default indicates,
  like \code{plotselect = c(resfit = 2, leverage = 1)}.
  Adding an element \code{default = 0} suppresses all plots not
  mentioned. This is useful to select single plots, like
  \code{plotselect = c(resfit = 3, default = 0)}

  The names of \code{plotselect} refer to:
  \describe{
    \item{yfit}{response versus fitted values}
    \item{resfit}{residuals versus fitted values (Tukey-Anscombe plot)}
    \item{absresfit}{residuals versus fitted values, defaults to TRUE for
  	 ordinary regression, FALSE for glm and others}
    \item{absresweights}{residuals versus weights}
    \item{qq}{normal Q-Q plot, defaults to TRUE for
  	 ordinary regression, FALSE for glm and others}
    \item{leverage}{residuals versus leverage (hat diabgonal)}
    \item{resmatrix}{scatterplot matrix of residuals for
  	 multivariate regression}
    \item{qqmult}{qq plot for Mahlanobis lengths versus sqrt of chisquare
  	 quantiles.}
  }

  In the 'resfit' (Tukey-Anscombe) plot, the reference line indicates
  a "contour" line with constant values of the response variable,
  \eqn{Y=\widehat y+r=} constant. It has slope \code{-1}.
  It is useful to judge whether any curvature shown by the smooth
  might disappear after a nonlinear, monotone transformation of the
  response.

  If \code{smresid} is true, the 'absresfit' plot uses modified
  residuals: differences between the ordinary residuals and the smooth
  appearing in the 'resfit' plot.
  Analogously, the 'qq' plot is then based on yet another modification
  of these modified residuals: they are scaled by the smoothed scale
  shown in the 'absresfit' plot, after these scales have been
  standardized to have a median of 0.674 (=\code{qnorm(0.75)}).

  The smoothing function used by default is \code{\link{smoothRegr}},
  which calls \code{\link{loess}}. This can be changed by setting
  \code{ploptions(smooth.function=<func>)}, which must have the same
  arguments as \code{smoothRegr}.
  
  The arguments \code{lty, lwd, colors} characterize how the graphical
  elements in the plot are shown.
  They should be three vectors of length 9 each, defining the
  line types, line widths, and colors to be used for ...
  \describe{
    \item{[1] }{observations;}
    \item{[2] }{reference lines;}
    \item{[3] }{smooth;}
    \item{[4] }{simulated smooths;}
    \item{[5] }{component effects in plresx;}
    \item{[6] }{confidence bands of component effects.}
    In the case of \code{glm.restype="cond.quant"}
    \item{[7] }{(random) observations;}
    \item{[8] }{conditional medians;}
    \item{[9] }{bars showing conditional quantiles.}
  }

  If smooths are shown according to groups (given in
  \code{smooth.group}), then a legend can be required and positioned
  in the respecive plots by using the argument \code{smooth.legend}.
  If it is \code{TRUE}, then the legend will be placed in the
  \code{"bottomright"} corner.
  Alternatively, the corner can be specified as
  "bottomright", "bottomleft", "topleft", or "topright".
  A coordinate pair may also be given.
  These possibilities can be used individually for each plot by
  giving a named vector or a named list, where the names are
  one of "yfit", "resfit", "absresfit", "absresweight", "(xvar)" or
  names of x variables provided by the \code{xvar} argument.
  A component "(xvar)" selects the first x variable.
}

\value{
  The list of the evaluations of all arguments and some more useful
  items is returned invisibly.
}
%\references{}
\author{Werner A. Stahel, ETH Zurich}
\note{
  This is a function under development. Future versions may behave
  differently and may not be compatible with this version.
}
\seealso{plotregr.control, plot.lm}
\examples{
r.savings <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
plot.regr(r.savings)

data(d.blast)
r.blast <-
     lm(log10(tremor)~location+log10(distance)+log10(charge),
          data=d.blast)
plot.regr(r.blast, sequence=TRUE, transformed=TRUE) 

data(d.fossileSamples)
r.foss <-
  lm(cbind(sAngle,lLength,rWidth)~SST+Salinity+lChlorophyll+Region+N,
  data=d.fossileSamples)
plot.regr(r.foss, plotselect=c(resfit=3, resmatrix=1, qqmult=1))

data(d.babysurvival)
rr <- glm(Survival~Weight+Age+Apgar1, data=d.babysurvival, family=binomial)
plot.regr(rr, xvar= ~Weight, cex.plab=0.7, ylim=c(-5,5))
plot.regr(rr, condquant=FALSE)

data(housing, package="MASS")
rr <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
plot.regr(rr, factor.show="jitter")

}
\keyword{hplot}
\keyword{regression}
