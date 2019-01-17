\name{plotregr.control}
\alias{plotregr.control}
\title{
Further Arguments to \code{plot.regr}
}
\description{
  Specify some arguments of minor importance for the function
  \code{\link{plot.regr}}
}
\usage{
plotregr.control(x, formula = NULL, data = NULL, rawvars = TRUE,
  smooth.group = NULL, weights = NULL, plweights = NULL,
  plab = NULL, pch = NULL, pcol = NULL, cex.plab = -1, cex.pch = -1,
  markprop = NULL,
  glm.restype = "working", smresid = TRUE,
  jitter = TRUE, jitterbinary = TRUE, mbox = TRUE,
  partial.resid = TRUE, cookdistlines = 1:2, leveragelim = c(0.99, 0.5),
  smooth = 2, smooth.par = NA, smooth.iter = NULL, smooth.sim = NULL,
  reflines = TRUE, reflines.se = FALSE,
  ylim = TRUE, ylimfac = c(4, 3), ylimext = 0.1,
  plext = getUserOption("plext"), yaxp = NULL,
  lty = c(1, 2, 5, 3, 4, 6, 1, 1), lwd = c(1, 1, 2, 1, 1.5, 1, 1, 1),
  colors.ra = getUserOption("colors.ra"),
  colors.smgrp = getUserOption("colors"), smooth.pale = 0.2,
  main = NULL, cex.title = NULL, xlabs = NULL, ylabs = NULL,
  mf = TRUE, mfcol = FALSE, multnrow = 0, multncol = 0, oma = NULL,
%%  ask = NULL,
  fit = TRUE, hat = TRUE, stresiduals = TRUE)
}
\arguments{
  \item{x}{regr %(or lm or glm)
    object (result of a call to \code{regr}). This is the only argument
    that is needed. All others have useful defaults.}
  \item{formula, data, rawvars}{see \code{?\link{plot.regr}}}
  \item{smooth.group}{grouping factor. If given, smooth will be drawn
    for each group rather than for the whole dataset.}
  \item{weights, plweights}{if TRUE (or NULL or NA), \code{x$weights} will be
    used if available as weights for calculating smooths and optionally
    be used to plot weighted residuals against it, and as sizes of
    symbols that show the data (residuals), respectively.
    Alternatively, a vector of weights can be specified.
    If FALSE, no weights will be used.}
  \item{plab}{labels of observations used to identify them, depending on
    \code{markprop}.}
  \item{pch}{plotting character to use for unlabeled points. Defaults to
    3 (corresponding to a plus sign) for less that 100 observations, and
    to "." otherwise -- unless residuals are of class \code{condquant},
    when 15 (filled square) is used.}
  \item{pcol}{color(s) to use for points}
  \item{cex.plab}{character size for data labels}
  \item{cex.pch}{symbol size for data to be shown by \code{pch}}
  \item{markprop}{proportion of extreme residuals to be labeled.
    If all points should be labeled, let \code{markprop=1}.}
  \item{glm.restype}{type of residuals to be used for glm models.
    In addition to those allowed in \code{residuals()} for
    \code{glm} objects, type \code{cond.quant} is possible for
    (ungrouped) binary regression. See \code{?\link{residuals.polr}} for
    an explanation.
    Warning: type "deviance" will not work with simulated smooths since
    NAs will emerge.}
  \item{smresid}{logical: Should residuals from smooth be used for
    'tascale' and 'qq' plots?}
  \item{jitter}{amount of jitter to be used for plotting against factors}
  \item{jitterbinary}{
will disappear
}
  \item{mbox}{logical: For plotting residuals versus factors, should
    \code{\link{plmboxes}} be used?}
  \item{partial.resid}{
%%     ~~Describe \code{partial.resid} here~~
}
  \item{cookdistlines}{levels of Cook distance for which contours are
    plotted in the leverage plot}
  \item{leveragelim}{bound for leverages to be used in standardizing
    residuals and in calculation of standardized residuals from smooth
    (if \code{smresid} is \code{TRUE}).}
  \item{smooth}{if TRUE (or 1), smooths are added to the plots where
    appropriate. If \code{==2}, smmooths to positive and negative
    residuals-from-smooth are also shown.}
  \item{smooth.par}{argument to the smooth function. In case of default
    smooth function, it is the argument \code{span} of \code{loess}}
  \item{smooth.iter}{if the loess smooth is used,
    this argument controls robustness of smooth. If \code{==0}, least
    squares is used, otherwise, robust smoothing is,
    \code{smooth.iter} being the maximal number of iterations.
    Defaults to 0 if the family is \code{binomial} or
    \code{multinomail}, to 50 otherwise.}
  \item{smooth.sim}{number of simulated smooths added to each plot.
    If NULL (the default) 19 simulated smooths will be generated if
    possible and sensible (i.e., none if \code{smooth.group} is set).}
  \item{reflines}{logical: Should reference lines showing
    (approximately) constant response values be shown?
    If \code{addcomp} is \code{TRUE}, should component effects be
    shown?}
  \item{reflines.se}{logical: Should standard error bands be shown
    around reference lines?}
  \item{ylim}{Limits to be used for vertical axes for setting an
    "inner plotting range". Can be a vector of length 2 or a matrix
    with columns giving the ranges for each plot. [needs more
    explanation.] !!!
    If \code{TRUE}, they will be determined in
    the default way, using \code{\link{robrange}}.}
  \item{ylimfac}{factor to be used for determining the (default)
    inner plotting range, see \code{?\link{robrange}}}
  \item{ylimext}{extension of the inner plotting range that determines
    the outer plotting range}
  \item{plext}{extension factors for the four sides of the plots,
    producing a "padding" area around the plotted points.
    This corresponds to setting \code{par(xaxs="e", yaxs="e")}, but
    allows to set different values for the 4 sides.}
  \item{yaxp}{axis parameters for the vertical axis, as in
    \code{?\link{par}}.}
  \item{lty, lwd, colors.ra}{three vectors of length 9 each,
    defining the line types, line widths, and colors to be used for ...
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
    } }
  \item{colors.smgrp}{colors to be used for distinguishing groups as
    specified by \code{smooth.group}}
  \item{smooth.pale}{paling applied to the color(s) of the smooth lines
    to distinuish simulated lines from the line corresponding to the data}
  \item{main}{main title to be used on each page,
    defaults to the formula of object. If \code{main} starts by ":",
    this colon is dropped and the formula is appended to the remainder.
    If \code{getUserOption("doc")} is $>1$, then the \code{tit} attribute of
    \code{x} (see \code{\link{tit}}) is printed on a separate line.}
  \item{cex.title}{character size for title}
  \item{xlabs}{labels for x variables. Defaults to \code{vars}}
  \item{ylabs}{labels for vertical axes}
  \item{mf}{vector of 2 elements, indicating the number of rows and
    columns of panels on each plot page.
    Defaults to \code{c(2,2)}, except for multivariate models, where
    it adjusts to the number of target variables.
    \code{mf=c(1,1)} or \code{mf=1} asks for a single frame per page.
    \code{mf=NA} or \code{mf=0} leaves the framing (and \code{oma})
    unchanged.}
  \item{mfcol}{if TRUE, the panel will be filled columnwise}
  \item{multnrow, multncol}{number of rows and columns of panels on
    one page, for residuals of multivariate regression only}
  \item{oma}{vector of length 4 giving the number of lines in the outer
    margin. If it is of length 2, they refer to top an right margins.}
%%-   \item{ask}{logical, as in \code{?\link{par}}}
  \item{fit, hat, stresiduals}{logicals: should these items be generated?}
}
%%- \details{
%%- %%  ~~ If necessary, more details than the description above ~~
%%- }
\value{A list containing all the items needed to specify plotting
  in \code{\link{plot.regr}} and \code{\link{plresx}}
}
\note{This function is not explicitly called by the user, but by
  \code{\link{plot.regr}} and \code{\link{plresx}}.
  All the arguments specified here can and should be given as
  arguments to these functions.
}
  
\author{
WS
}
\seealso{
  \code{\link{plot.regr}} and \code{\link{plresx}}
}
\examples{
data(d.blast)
( r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast) )

plargs <- plotregr.control(r.blast, formula = ~.+distance, rawvars=FALSE,
smooth.group = location )
showd(plargs$pldata)
names(plargs)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

