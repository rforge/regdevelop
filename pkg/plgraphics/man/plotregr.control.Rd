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
plotregr.control(x, data = NULL, xvar = TRUE, transformed = FALSE,
  weights = NULL,
  glm.restype = "working", condquant = TRUE, smresid = TRUE,
  partial.resid = TRUE, cookdistlines = 1:2,
  leveragelim = c(0.99, 0.5), condprobRange = NULL,
  refline = TRUE, reflineband = FALSE, testlevel = 0.05,
  smooth = 2, %% smoothPar = NA, smoothIter = NULL,
  smooth.sim = NULL,
  xlabs = NULL, reslabs = NULL, markextremes = NULL,
  mf = TRUE, mfcol = FALSE, multnrow = 0, multncol = 0, oma = NULL, ...)
} %%    jitter = TRUE, jitterbinary = TRUE, mbox = TRUE,
\arguments{
  \item{x}{regr %(or lm or glm)
    object (result of a call to \code{regr}). This is the only argument
    that is needed. All others have useful defaults.}
  \item{data}{see \code{?\link{plot.regr}}}
  \item{xvar}{variables for which residuals shall be plotted.
    Either a formula like \code{~ x1 + x2} or a character vector of
    names. Defaults to all variables (or terms, see \code{transformed})
    in the model.
  }    
  \item{transformed}{see \code{?\link{plot.regr}}}
  \item{weights}{logical: should residuals be plotted against weights?
    Used in \code{plresx}.}
  \item{glm.restype}{type of residuals to be used for glm models.
    In addition to those allowed in \code{residuals()} for
    \code{glm} objects, type \code{condquant} is possible for
    (ungrouped) binary regression.
    See \code{?\link{residuals.regrpolr}} for an explanation.
    Warning: type "deviance" will not work with simulated smooths since
    NAs will emerge.}
  \item{condquant}{logical: should conditional quantiles be shown for
    censored observations, binary and ordered responses?}
  \item{smresid}{logical: Should residuals from smooth be used for
    'tascale' and 'qq' plots?}
  \item{partial.resid}{if TRUE (default),
    residuals for observations will be displayed.}
  \item{cookdistlines}{levels of Cook distance for which contours are
    plotted in the leverage plot}
  \item{leveragelim}{bound for leverages to be used in standardizing
    residuals and in calculation of standardized residuals from smooth
    (if \code{smresid} is \code{TRUE}).}
  \item{condprobRange}{numeric vector of length 2.
    In the case of residuals of class \code{condquant}, quartile bars
    are only drawn for residuals with probability between
    \code{condprobRange[1]} and \code{condprobRange[1]}.
    Default is \code{c(0.05,0.8)} for less than 50 observations, and
    \code{c(0,0)}, suppressing the bars, otherwise.}
  \item{refline}{logical: should reference line be shown?}
  \item{reflineband}{logical: should a confidence band be drawn for the
    component effects?}
  \item{testlevel}{level for statistical tests}
  \item{smooth}{if TRUE (or 1), smooths are added to the plots where
    appropriate. If \code{==2}, smmooths to positive and negative
    residuals-from-smooth are also shown.}
%%-   \item{smoothPar}{argument to the smooth function. In case of default
%%-     smooth function, it is the argument \code{span} of \code{loess}}
  \item{smooth.sim}{number of simulated smooths added to each plot.
    If NULL (the default) 19 simulated smooths will be generated if
    possible and sensible (i.e., none if \code{smooth.group} is set).}
  \item{xlabs}{labels for x variables. Defaults to \code{vars}}
  \item{reslabs}{labels for vertical axes}
  \item{markextremes}{proportion of extreme residuals to be labeled.
    If all points should be labeled, let \code{markextremes=1}.}
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
  \item{...}{further arguments in the call, to be ignored by
    'plotregr.control'}
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
  
\author{Werner A. Stahel, Seminar for Statistics, ETH Zurich}

\seealso{
  \code{\link{plot.regr}} and \code{\link{plresx}}
}
\examples{
data(d.blast)
( r.blast <-
  lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast) )

plargs <- plotregr.control(r.blast, formula = ~.+distance, transformed=TRUE,
smooth.group = location )
showd(plargs$pdata)
names(plargs)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

