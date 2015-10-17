\name{plresx}
\alias{plresx}
\title{ Plot Residuals Against Explanatory Variables }
\description{
  Diagnostic plots for fitted regression models:
  Residuals are plotted against explanatory variables,
  including smooth for actual residuals and simulated ones,
  as well as a reference reference line for constant values of the
  response variable.  
}
\usage{
plresx(x, data = NULL, resid=NULL, vars = NULL, sequence=FALSE,
            se = FALSE, partial.resid = TRUE, addcomp = FALSE, 
            glm.restype = "deviance", condprobrange=c(0.05,0.8),            
            weights = NULL, wsymbols=NULL, symbol.size=NULL,
            markprop=NULL, lab = NULL, cex.lab = 0.7,
            reflines = TRUE, mbox = FALSE, jitter=NULL, jitterbinary=TRUE,
            smooth = TRUE, smooth.par=NA, smooth.iter=NULL,
            smooth.sim=19, nxsmooth=51, smooth.group=NULL,
            smooth.col=NULL, smooth.pale = 0.2, smooth.legend=TRUE, 
            multnrows = 0, multncols = 0, 
            lty = c(1,2,5,3,6,4,1,1), lwd=c(1,1,2,1,1.5,1,1,1),
            colors = getUserOption("colors.ra"), pch=NULL, col=NULL,
            rug = FALSE, 
            main = NULL, cex.title = NULL,  xlabs = NULL, ylabs = NULL, 
            ylim=TRUE, ylimfac=4.0, ylimext=0.1, yaxp=NULL, 
            cex = par("cex"), ask = NULL,
            ...)
}
\arguments{
  \item{x}{fitted regression object, result of calling \code{link{regr}}}
  \item{data}{data set. defaults to the data set used when calling \code{regr}}
  \item{resid}{residuals. If not \code{NULL}, these will be used instead
    of \code{residuals(x)}}
  \item{vars}{variables against which the residuals shall be plotted.
    Defaults to the variables contained in the formula of \code{object}.
    \code{vars} can be a formula (of which all variables on the right
    hand side will be used) or a character vector of variable names.
    All variables must be contained in \code{data}}
  \item{sequence}{if TRUE, the residuals are also plotted against their
    sequence number}
  \item{se}{if TRUE, standard errors of the fit are plotted}
  \item{partial.resid}{if TRUE (default),
    residuals for observations will be displayed}
  \item{addcomp}{if TRUE, the component effect is added to the residuals
    as is common for a "partial residual plot"}
  \item{glm.restype}{type of residuals to be plotted in case of a glm
    model. In addition to those allowed in \code{residuals()} for
    \code{glm} objects, type \code{cond.quant} is possible for
    (ungrouped) binary regression. See \code{?\link{residuals.polr}} for
    an explanation.}
  \item{condprobrange}{if conditional probabilities are shown by bars:
    show only those that are contained in \code{condprobrange}
    (vector of length 2)}
  \item{weights}{weights to be used for smoothing as well as
    plotting of observations without label}
  \item{wsymbols}{if TRUE, weights are used as symbol size}
  \item{symbol.size}{maximal symbol size, relative to \code{0.02*par("pin")[1]}}
  \item{markprop}{proportion of extreme residuals to be labeled}
  \item{lab}{character string labels to be used for displaying observed
    residuals}
  \item{cex.lab}{character expansion for labels}
  \item{reflines}{logical: Should reference lines be added to the plots
    of residuals on x? Default is TRUE if available.}
  \item{mbox}{For plotting residuals versus factors, should
    \code{\link{plmboxes}} be used?}
  \item{jitter}{amount of jitter to be added for displaying factors}
  \item{jitterbinary}{logical: should variables with only 2 distinct values be
    treated as factors, i.e., jittered?}
  \item{smooth}{if TRUE, a smooth line will be fitted to each plot.
    By default, \code{\link{loess}} is used. Alternatively,
    another smoothing function can be given as argument \code{smooth}}
  \item{smooth.par}{parameter driving the smooth: argument \code{span}
    if \code{loess} is used, which defaults to
    \code{3*nobs^log10(1/2)*(1+glm)}, where \code{glm==1} for glm models}
  \item{smooth.iter}{if the loess smooth is used,
    this argument controls robustness of smooth. If \code{==0}, least
    squares is used, otherwise, robust smoothing is,
    \code{smooth.iter} being the maximal number of iterations.
    Defaults to 0 if the family is \code{binomial} or
    \code{multinomail}, to 50 otherwise.}
  \item{smooth.sim}{number of simulated smooth lines to be drawn}
  \item{nxsmooth}{number of x values for which the smooth is calculated}
  \item{smooth.group}{grouping factor. If given, smooth will be drawn
    for each group rather than for the whole dataset.}
  \item{smooth.col}{colors to be used for the different groups}
  \item{smooth.pale}{argument to \code{\link{colorpale}} for attenuating
    the color of simulated smooths with respect to 'original' smooths}
  \item{smooth.legend}{if a legend for the groups is desired, its
    position should be given here}
  \item{multnrows, multncols}{number of rows and columns for multiple
    frames in case of multivatiate residuals}
  \item{lty, lwd, colors}{vector of line types and colors to be used for \cr
    \describe{
      \item{[1] }{observations (only \code{colors[1]} matters);}
      \item{[2] }{horizontal reference line;}
      \item{[3] }{smooth;}
      \item{[4] }{simulated smooths;}
      \item{[5] }{reference lines related to component effect;}
      \item{[6] }{for confidence bands of component effects;}
      \item{[7] }{(rnadom) observations in the case of
	\code{glm.restype="cond.quant"};}
      \item{[8] }{bars showing conditional quantiles for
	\code{glm.restype="cond.quant"};}
      }
%%-     [1]         observations (only \code{colors[1]} matters),\cr
%%-     [2]         horizontal reference line ,\cr
%%-     [3]         smooth,\cr
%%-     [4]         simulated smooths,\cr
%%-     [5]         reference lines,\cr
%%-     [6]         confidence bands of terms.\cr
  }
  \item{pch}{plotting character to use for unlabeled points}
  \item{col}{color(s) to use for points}
  \item{rug}{if TRUE, the "rug" is displayed. It consists of tick marks
    for the observed values.}
  \item{main}{main title}
  \item{cex.title}{character expansion for main title}
  \item{xlabs}{labels of x variables to be used instead of those
    given by \code{vars}}
  \item{ylabs}{label for residuals}
  \item{ylim}{if TRUE, the range of the y axis will be restricted by
    calling \code{\link{robrange}}}
  \item{ylimfac}{argument passed to \code{robrange} if selected}
  \item{ylimext}{argument passed to \code{plcoord}}
  \item{yaxp}{axis parameters used for residuals, passed to \code{plot},
    see \code{?par}.}
  \item{cex}{character expansion}
  \item{ask}{Should R ask before erasing the screen?
    The default is 'getOption('ask'), followed by TRUE for certain
    graphical devices.}
  \item{\ldots}{further arguments passed to \code{plot}}
}
\details{
  The plots have several unusual features:
  \itemize{
    \item
    Showing outliers in the plot margins and marking them with their
    labels (if \code{markprop>0}).
    \item
    A reference line is shown. It follows a "contour" of an approximately
    constant value of the response variable -- more precisely,
    it shows the "component effect" with opposite sign.
    It helps decide if a curvature shown by the smooth can be removed
    by transforming the explanatory variable forming the horizontal
    axis. The component effect is defined as the curve of fitted values
    obtained by varying the explanatory variable, keeping all the other
    variables at their "central value" (the mean of continuous variables
    and the mode of factors).
}
}
\value{
  The names of the variables against which the residuals are plotted are
  returned invisibly
}
\references{There will be a website}
\author{Werner A. Stahel, ETH Zurich}
\note{This function is still under construction
}
\seealso{\code{\link{plot.regr}}, \code{\link{termplot}}}
\examples{
  r.savings <- regr(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings)
  plot(r.savings)

  data(d.blast)
  r.blast <-
       regr(log10(tremor)~location+log10(distance)+log10(charge),
            data=d.blast)
  plresx(r.blast, xplot=~.+no) }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{regression}% __ONLY ONE__ keyword per line
