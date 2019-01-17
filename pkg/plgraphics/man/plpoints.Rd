\name{plpoints}
\alias{plpoints}
\alias{plbars}
\alias{plframe}
\alias{pltitle}
\alias{plaxis}
\alias{plsmooth}
\alias{plsmoothline}
\alias{plrefline}
\alias{plmark}

\title{Low level plotting functions for the 'pl' system
}
\description{
  These functions do the plotting of points and enhancing lines
  based on the 'pl' paradigm
}
\usage{
plframe(x, y, xlab = NULL, ylab = NULL, plextext = NULL,
  axcol = rep(1, 4), plargs = NULL, ploptions = plargs$ploptions)

pltitle(main=NULL, sub=NULL, cex=NULL, cexmin=NULL, 
  side=3, line=NULL, adj=NULL, outer.margin=NULL, col="black",
  doc=NULL, show=TRUE, plargs=NULL, ploptions = plargs$ploptions, ...)

plaxis(axis, x, lab=TRUE, range=NULL, varlabel=NULL, col=1,
  tickintervals=NULL, plargs = NULL, ploptions = plargs$ploptions, ...) 

plpoints(x, y, type="p", plab=NULL, pch=NULL, col=NULL, psize=NULL, lty=NULL,
  condquant = TRUE, plargs = NULL, ploptions = plargs$ploptions)

plbars(x, y, plargs = NULL, ploptions = plargs$ploptions)

plsmooth(x, y, band=NULL, plargs = NULL, ploptions = plargs$ploptions)
plsmoothline(smoothline, x, y, plargs=NULL, ploptions = plargs$ploptions)

plrefline(refline, x=NULL, innerrange=NULL, y=NULL,
  cutrange = c(x = TRUE, y = FALSE), plargs=NULL, ploptions=plargs$ploptions)

plmark(x, y=NULL, markextremes=NULL, plabel=NULL,
  plargs=NULL, ploptions=plargs$ploptions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{coordinates for the horizontal axis
  }
  \item{y}{coordinates for the vertical axis
  }
  \item{plab}{labels for displaying points. Overrides labels provided by
    \code{plargs$pdata[["plab"]]}.
  }
  \item{axcol}{colors for drawing axes scales
  }
  \item{xlab, ylab}{axis labels}
  \item{main, sub}{main title and subtitle}
  \item{varlabel}{variable name}
  \item{plextext}{extension of the plotting area beyond the possibly
    already extended range. This helps to accomodate sybmbols.
  }
  \item{side}{For \code{pltitle}: in which margin should the text
    be shown?}
  \item{cex}{character size. May be vector of length 3, giving size for
    main title, subtitle, and \code{tit} attribute of title,
    respectively. The default is given by \code{ploptions("title.cex")}.
  }
  \item{cexmin}{minimal character size, to be used to adjust the
    character size to the length of the text (if \code{cex} is \code{NULL})}
  \item{line}{line in margin on which the main title is placed -- or the
    subtitle if \code{main} is \code{NULL}}
  \item{adj}{text adjustment, scalar between 0 and 1}
  \item{outer.margin}{logical: should title text be placed in outer margin?}
  \item{col}{color for the title text or axis line and tickmarks}
  \item{axis}{integer indicating which axis is to be drawn}
  \item{lab}{logical: should labels for tickmarks be displayed?}
  \item{range}{range in which tickmarks are set}
  \item{doc}{logical: should the \code{tit} attribute of \code{main}
    be displayed if available?}
  \item{show}{logical: if \code{FALSE}, nothing will be done if
    there are multiple frames and the current one is not the first.
    If it is negative, no title will be shown, but the value will be
    returned.
  }
  \item{tickintervals}{number of intervals used by
    \code{\link{pretty}} to determine the axis ticks.
  }
  \item{type}{type of displaying points. See \code{?\link{points}}.}
  \item{pch, psize, lty}{... and \code{col} in \code{plpoints}:
    plotting character(s), size and color of plotting points,
    and line type. %% to be used for \code{plpoints}.
    Overrides other settings, defined in \code{plargs}.
  }
  \item{condquant}{logical: should bars be drawn for censored data?
    If \code{FALSE}, censored observations will be set to the median of
    the conditional distribution and shown by a different plotting
    character, see argument \code{censored} of
    \code{\link{ploptions}}.
    If \code{NULL}, the standard plotting character will be used.
  }
  \item{smoothline}{for \code{plsmoothline}: result of a smooth fitting}
  \item{band}{logical: should a band (e.g., a confidence band)
    be drawn together with the smooth?
  }
  \item{refline}{for \code{plrefline}:
    A vector giving intercept and slope of a straight line, 
    or a function that returns this list, such as \code{lm}}
  \item{innerrange}{for \code{plrefline}: inner range in x direction -
    only needed of the refline should be clipped at a range different
    from the \code{innerrange} attribute of the horizontal variable}
  \item{cutrange}{for \code{plrefline}: logical vector of length 2:
    should the reference line(s) be cut at the inner plotting ranges in
  x- and y-direction?
    Otherwise, it will be continued outside it with the appropriate
    transformation.
  }
  \item{markextremes}{proportion of observations that should be labeled
    because they show extreme values for a variable
  }
  \item{plabel}{point labels to be used to mark extreme points
  }
  \item{plargs, ploptions}{result of \code{\link{pl.control}}, see Details
  }
  \item{\dots}{absorbs extra arguments}
}
\value{\code{pltitle} invisibly returns a character vector containing
  the main and sub title.\cr
  \code{plsmooth} invisibly returns the data.frame needed for
  drawing the smooth line.
}
\details{
  \code{plframe} draws axes according to argument \code{axes},
  by calling \code{plaxis}.
  It looks for attributes of \code{x} and \code{y}, such as
  \code{innerrange} and \code{ticksat}.
  Tick labels are shown at the values of the \code{ticklabelsat}
  attribute if available, otherwise at the values of \code{ticksat}.
  The labels can be given by the attribute \code{ticklabels}.
  This facilitates setting more tick marks than labels, see the
  example.\cr
  It also draws a grid.
  The positions of gridlines at \code{ticksat} by default.\cr
  Finally, it draws "zero" lines as determined by the pl option
  \code{zeroline}. The latter can be a numeric vector giving
  the positions of such threshold lines, or a list of two such vectors,
  the first for horizontal axis, the second for the vertical axis.
  
  \code{plargs} and \code{ploptions} may be specified explicitly,
  but are usually generated by calling \code{pl.control}.
}
  
\author{Werner A. Stahel
}

\seealso{\link{plsmooth}; \link{plrefline}; \link{gentimeaxis}; 
    \link{pl.control}
}
\examples{
plyx(Sepal.Width ~ Sepal.Length, data=iris)

## again, each step separately
t.dt <- .plargs$pldata
.plargs$grid <- TRUE
plframe(t.dt$Sepal.Length, t.dt$Sepal.Width, plargs=.plargs)
plsmooth(t.dt$Sepal.Length, t.dt$Sepal.Width, plargs=.plargs)
t.plab <- plmark(t.dt$Sepal.Length, t.dt$Sepal.Width, markextremes=0.03,
                 plargs=.plargs)
plpoints(t.dt$Sepal.Length, t.dt$Sepal.Width, plargs=.plargs, plab=t.plab) 
}

\keyword{aplot}

