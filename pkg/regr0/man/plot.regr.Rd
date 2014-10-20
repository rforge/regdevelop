\name{plot.regr}
\alias{plot.regr}
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
\method{plot}{regr}(x, data = NULL, markprop = NULL, lab = NULL, cex.lab = 0.7,
  mf = NULL, mfcol=FALSE, mar=c(3,3,2,1), mgp=c(2,0.7,0),
  oma = 2*(prod(mf)>1), cex=par("cex"), ask = NULL, 
  multnrows = 0, multncols = 0,
  lty = c(1, 2, 5, 3, 4, 6, 1, 1), lwd = c(1, 1, 2, 1, 1.5, 1, 1, 1),
  colors = getUserOption("colors.ra"), pch = NULL, col = NULL,
  main = NULL, cex.title = NULL,
  wsymbols = NULL, symbol.size = NULL,
  smooth = TRUE, smooth.par = NA, smooth.iter = NULL, smooth.sim = 19,
  nxsmooth = 51,
  plotselect = NULL, 
  weights = NULL, leverage.cooklim = 1:2,
  reslim = TRUE, reslimfac=4, reslimext=0.1,
  yaxp=NULL, resaxp=NULL, stresaxp=NULL,
  ylim = TRUE, ylimfac=3.0, ylimext=0.1, jitterbinary=TRUE, 
  glm.restype = "deviance", condprobrange=c(0.05,0.8),
  sequence = FALSE, xplot = TRUE, x.se = FALSE, x.smooth = smooth,
  addcomp=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{regr %(or lm or glm)
    object (result of a call to \code{regr}). Thiw iw the only argument
    that is needed. All others have useful defaults.}
  \item{data}{data set where explanatory variables are found}
  \item{markprop}{proportion of extreme residuals to be labeled}
  \item{lab}{labels of observations used to identify them}
  \item{cex.lab}{character size for lab. If the length is $>1$, the
    second element determines the size of non-labeled points,
    otherwise, their size is adjusted to the number of points
    }
  \item{mf}{vector of 2 elements, indicating the number of rows and
    columns of panels on each plot page.
    Defaults to \code{c(2,2)}, except for multivariate models, where
    it adjusts to the number of target variables}
  \item{mfcol}{if TRUE, the panel will be filled columnwise}
  \item{mar, mgp, oma, cex, ask}{see \code{?par}}
  \item{multnrows, multncols}{number of rows and columns of panels on
    one page, for residuals of multivariate regression only}
  \item{lty, lwd, colors}{three vectors of length 6 each, defining the
    line types, line widths, and colors to be used for ...
    \describe{
      \item{[1] }{observations;}
      \item{[2] }{reference lines;}
      \item{[3] }{smooth;}
      \item{[4] }{simulated smooths;}
      \item{[5] }{component effects in plresx;}
      \item{[6] }{for confidence bands of component effects.}
      \item{[7] }{(random) observations in the case of
	\code{glm.restype="cond.quant"};}
      \item{[8] }{bars showing conditional quantiles for
	\code{glm.restype="cond.quant"};}
    }
  }
  \item{pch}{plotting character to use for unlabeled points. Defaults to
    3 (corresponding to a plus sign) unless residuals are of class
    \code{condquant}}
  \item{col}{color(s) to use for points}
  \item{main}{main title to be used on each page,
    defaults to the formula of object. If \code{main} starts by ":",
    this colon is dropped and the formula is appended to the remainder.
    If \code{getUserOption("doc")} is $>1$, then the \code{tit} attribute of
    \code{x} (see \code{\link{tit}}) is printed on a separate line.}
  \item{cex.title}{character size for title}
  \item{wsymbols}{if TRUE, points are displayed by a bubble whose size
    indicates the weight as given in \code{object$weights}.
    If it is a vector of appropriate size, it will be used to determine
    the sizes of such bubbles.}
  \item{symbol.size}{determines the size of weight symbols}
  \item{smooth}{if TRUE, smooths are added to the plots where
    appropriate. If it is a function, it is assumed to define the smooths.
  }
%%  \item{smooth.rob}{if TRUE, smooth will be done robustly.
%%    Default is FALSE for binomial or multinomial regression.}
  \item{smooth.par}{argument to the smooth function. In case of default
    smooth function, it is the argument \code{span} of \code{loess}}
  \item{smooth.iter}{if the loess smooth is used,
    this argument controls robustness of smooth. If \code{==0}, least
    squares is used, otherwise, robust smoothing is,
    \code{smooth.iter} being the maximal number of iterations.
    Defaults to 0 if the family is \code{binomial} or
    \code{multinomail}, to 50 otherwise.}
  \item{smooth.sim}{number of simulated smooths added to each plot}
  \item{nxsmooth}{number of values on the horizontal axes where smooths
    are evaluated}
  \item{plotselect}{which plots should be shown? See Details}
  \item{weights}{if TRUE, residuals will be plotted versus
    \code{object$weights}. Alternatively, a vector of weights can be specified}
  \item{leverage.cooklim}{levels of Cook distance for which contours are
    plotted in the leverage plot}
  \item{reslim}{limits for plotting residuals directly. Beyond these
    limits, points are shown in the margins with squeezes scale}
  \item{reslimfac, reslimext}{arguments passed to \code{robrange} if
    selected}
  \item{resaxp, stresaxp}{axis parameters used for 
    residuals and standardized residuals, respectively, passed to
    argument \code{yaxp} of \code{plot}, see \code{?par}.}
  \item{ylim, ylimfac, ylimext, yaxp}{same, for plotting y values}
  \item{jitterbinary}{logical: should variables with only 2 distinct values be
    treated as factors, i.e., jittered?}
  \item{glm.restype}{type of residuals to be used for glm model.
    In addition to those allowed in \code{residuals()} for
    \code{glm} objects, type \code{cond.quant} is possible for
    (ungrouped) binary regression. See \code{?\link{residuals.polr}} for
    an explanation.}
  \item{condprobrange}{if conditional probabilities are shown by
    vertical bars:
    show only those that are contained in \code{condprobrange}
    (vector of length 2)}
  \item{sequence}{logical: should residuals  be plotted versus sequence?}
  \item{xplot}{if TRUE, residuals will be plotted versus all (raw)
    explanatory variables in the model (by calling \code{plresx})} 
  \item{x.se}{logical: Should confidence bands be drawn in plot of residuals
    versus explanatory variables?}
  \item{x.smooth}{like \code{smooth}, for the latter plots}
  \item{addcomp}{ logical: should component effects be added to
    residuals for residuals versus input variables plots? }
  \item{\dots}{further arguments to be passed to \code{plot}}
}
\details{
  Argument \code{plotselect} is used to determine which plots will be
  shown. It should be a named vector of numbers indicating
     \describe{
       \item{0}{do not show}
       \item{1}{show without smooth}
       \item{2}{show with smooth}
     }
  The default is
  c( yfit=0, ta=3, tascale = NA, weights = NA, qq = NA, leverage = 2,
  resmatrix = 1, qqmult = 3).
     Modify this vector to change the selection and the sequence in
     which the plots appear.
     Alternatively, provide a named vector defining all plot that should
     be shown on a different level than the default indicates,
     like \code{plotselect = c(ta = 2, leverage = 1)}.
     Adding an element \code{default = 0} suppresses all plots not
     mentioned. This is useful to select single plots, like
     \code{plotselect = c(ta = 2, default = 0)}

     The names of \code{plotselect} refer to:
     \describe{
       \item{yfit}{response versus fitted values}
       \item{ta}{residuals versus fitted values (Tukey-Anscombe plot)}
       \item{tascale}{residuals versus fitted values, defaults to TRUE for
	 ordinary regression, FALSE for glm and others}
       \item{weights}{residuals versus weights}
       \item{qq}{normal Q-Q plot, defaults to TRUE for
	 ordinary regression, FALSE for glm and others}
       \item{leverage}{residuals versus leverage (hat diabgonal)}
       \item{resmatrix}{scatterplot matrix of residuals for
	 multivariate regression}
       \item{qqmult}{qq plot for Mahlanobis lengths versus sqrt of chisquare
	 quantiles.}
     }

     In the Tukey-Anscombe plot, the reference line indicates a "contour"
     line with constant values of the response variable,
     \eqn{Y=\widehat y+r=} constant. It has slope \code{-1}.
     It is useful to judge whether any curvature shown by the smooth
     might disappear after a nonlinear, monotone transformation of the
     response. 
}
\value{
none
}
%\references{}
\author{Werner A. Stahel, ETH Zurich}
\note{
  This is a function under development. Future versions may behave
  differently and may not be compatible with this version.
}
\seealso{plot.lm}
\examples{
  r.savings <- regr(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
  plot(r.savings)

  data(d.blast)
  r.blast <-
       regr(log10(tremor)~location+log10(distance)+log10(charge),
            data=d.blast)
  plot(r.blast, seq=FALSE, xplot=~.+no) 

  data(d.fossiles)
  r.foss <-
    regr(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+lChlorophyll+region+N,
    data=d.fossiles)
  plot(r.foss, plotselect=c(ta=3, resmatrix=1, qqmult=1),xplot=FALSE)
}
\keyword{hplot}
\keyword{regression}
