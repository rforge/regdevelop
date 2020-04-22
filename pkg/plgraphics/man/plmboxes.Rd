\name{plmboxes}
\alias{plmboxes}
\alias{plmboxes.default}
\alias{plmboxes.formula}
\alias{plmbox}

\title{
Multibox plots
}
\description{
  Draw multibox plot(s) for given (grouped) values, possibly asymmetric.
  'plbox' draws a single multibox plot (low level graphical function).
  'plboxes' is a high level graphics function that draws multiboxes for
  grouped data. A secondary, binary grouping factor can be given
  to produce asymmetric multiboxes.
}
\usage{
plmboxes(x, ...)
\method{plmboxes}{formula}(x, data, ...)

\method{plmboxes}{default}(x, y, data, width=1, at=NULL, horizontal=FALSE,
    probs=NULL, outliers=TRUE, na=FALSE, backback=NULL, refline=NULL,
    add=FALSE, xlim=NULL, ylim=NULL, axes=TRUE, xlab=NULL, ylab=NULL,
    labelsperp=FALSE, xmar=NULL,
    widthfac=NULL, minheight=NULL, colors=NULL, lwd=NULL, 
    .subdefault=NULL, plargs = NULL, ploptions = NULL, ...)

plmbox(x, at=0, probs=NULL, outliers=TRUE, na.pos=NULL, horizontal=FALSE,
    width=1, wfac=NULL, minheight=NULL, adj=0.5, extquant=TRUE, 
    widthfac=c(max=2, med=1.3, medmin=0.3, outl=NA),
    colors=c(box="lightblue2",med="blue",na="gray90"),
    lwd=c(med=3, range=2), warn=options("warn"))
}

\arguments{
  \item{x}{For 'plmboxes.formula':
    a formula, such as 'y ~ grp' or 'y~grp+grp2',
    where 'y' is a numeric vector of data values to be split into groups
    according to the grouping variable 'grp' (usually a factor)
    and, if given, according to the binary variable 'grp2'.
    'y~1+grp2' produces a single asymmetric mbox.

    For 'plmboxes.default':
    factor to be used as the grouping variable or matrix or data.frame
    with 2 columns (for asymmetric mbox plot), where the second column
    is binary.
  }
  \item{y}{a numeric vector of data values}
  \item{data}{a data.frame from which the variables in 'formula'
          should be taken.
}
  \item{width}{a vector giving the widths of the multibox plot
  for each group 'grp1'. 
}
  \item{at}{horizontal position of the multiboxes. Must have length equal
    to the number of (present) levels of the factor 'grp'.
    if an element of 'at' is 'NA', the group will be skipped.
    Defaults to 1, 2, ...
  }
  \item{horizontal}{logical. If TRUE, boxes will be drawn horizontally.
    Note that 'y' is then the horizontal coordinate, i.e., still the
    quantitative variable defining the boxes, and 'x' is still the
    grouping.
  }
  \item{probs}{probability values for selecting the quantiles.
    If all 'probs' are <=0.5, they will be mirrored at 0.5 for
    'plmboxes'. The default is c(0.05,0.25,0.5) if the average number of
    data per group (for 'plmboxes', or the number of data, for 'plmbox')
    is less than 20, c(0.025,0.05,0.125,0.25,0.375,0.5), otherwise. 
  }
  \item{outliers}{logical: should outliers be marked?
  }
  \item{na, na.pos}{if 'na' is not NULL, NA values will be represented by
    a box.
    If 'na' is TRUE, the position of the box will be generated to be below
    the minimum of the data.
    If 'na' (for 'plmboxes') or 'na.pos' (for 'plmbox') is a scalar or
    a vector of length 2,
    the position of the box is at that value (with a generated width) or
    between the 2 values, respectively.
  }
  \item{backback}{logical: Should two back-to-back multiboxes be
    displayed if the (single) x factor is binary?
  }
  \item{refline}{vertical positions of any horizontal reference lines
  }
  \item{add}{logical. If TRUE, the mboxes will be added to an existing
    plot without calling 'plot'.
  }
  \item{xlim}{plotting limits for the horizontal axis.
  }
  \item{ylim}{plotting limits for the vertical axis.
  }
  \item{axes}{logical. If FALSE, no axes are drawn.
  }
  \item{xlab}{label for the x axis. Defaults to the "x factor" --
    the first name on the right hand side of 'formula'
  }
  \item{ylab}{label for the x axis. Defaults to the left hand side of
    'formula'.
  }
  \item{labelsperp}{logical: Should the labels for the levels of the
    "x factor" be shown in perpendicular to the axis?
    If it is numeric, it determines the maximum label length, with a
    maximum of 20.
  }
  \item{xmar}{plot margin for the "x factor" axis.
    Default tries to be suitable, i.e. expand the margin if
    \code{labelsperp} is \code{TRUE} according to the length of the
    levels' labels. If \code{xmar} has two more elements,
    they determine the margin lines where the variable label and the levels'
    labels are shown.
  }
  \item{widthfac}{named vector used to modify the following settings:
    \cr  
    max=2: determes the maximal width of the boxes.
    Boxes that should be wider are censored and marked as such.
    \cr
    med=1.3, medmin=0.3: determine the width of the mark for the data median.
    The width is 'med' times the maximal width of the boxes,
    but at least 'medmin'.
    \cr
    outl=NA: length of the marks for outliers.
    \cr
    sep=0.003: width of the gap between the "half" mbox plots in case of
    asymmetric mboxes (only needed for 'plmboxes').
    \cr
    For 'plmboxes', the argument needs to contain only the elements that
    should be different from the default values.
  }
  \item{colors}{named vector or list selecting the colors to be used,
    with named elements:
    
    box="lightblue2": color with which the central box(es)
    (those corresponding to probabilieties between 0.25 and 0.75)
    will be filled.
    
    med="blue": color of the mark for the median.
    
    na="grey90": color with which the box for NA values will be filled.
    
    For 'plmboxes', the argument needs to contain only the elements that
    should be different from the default values.
  }
  \item{lwd}{named vector or list selecting the line width to be used,
    with named elements:
    
    med=3: line width for the mark showing the median
    
    range=2: line width for the line along the range of the data
    
    For 'plmboxes', the argument needs to contain only the elements that
    should be different from the default values.
  }
  \item{plargs}{result of calling \code{pl.control}.
    If \code{NULL}, \code{pl.control} will be called to generate it.
    If not null, arguments given in \code{...} will be ignored.
  }
  \item{ploptions}{list of pl options.
  }
  \item{\dots}{additional arguments passed to 'plot'
  }
  \item{.subdefault}{text for the subtitle in case that it is not
    specified
  }
Specific arguments for 'plmbox':
\item{wfac}{factor by which the widths of the boxes must be multiplied.
  If given, it overrides 'width'
}
\item{minheight}{minimal class width ("height") for the boxes (in case
  two quantiled are [almost] identical). The default is 0.02 times
  the (median of the within group) IQR.
}
\item{adj}{adjustment of the boxes. 'adj=0' leads to boxes aligned on
  the left, 'adj=1', on the right, 'adj=0.5', centered. Other values of
  'adj' make little sense.
}
\item{extquant}{logical, passed to \code{\link{quinterpol}}:
  Should the quantiles be extrapolated beyond the range of the data?
}
\item{warn}{level of warning for the case when there is no non-missing data}
}
\details{A multibox plot is a generalization (and modification) of the
  ordinary box plot that draws more details of the distribution in the
  form of a histogram with variable class widths.
  The classes are selected such that preselected quantiles form the
  class breaks. By default, these quantiled include the median and the
  quartiles, thereby recovering the box of the traditional box plot.
}
\value{
  \code{plmboxes} invisibly returns the 'at' values that are finally
  used.\\
  \code{plmbox} returns a factor [to be defined...!!!]
  and, as attributes, the quantiles and widths used to draw the boxes
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{Werner A. Stahel
}
%% \note{
%% %%  ~~further notes~~
%% }

\seealso{
\code{\link{boxplot}}
}
\examples{
plmboxes(Sepal.Length~Species, data=iris)

plmboxes(Sepal.Length~Species, data=iris,
  widthfac=c(med=2), colors=c(med="red"), horizontal=TRUE)

plmboxes(Sepal.Length~factor(Species)+I(Sepal.Width<=3), data=iris[1:100,],
         labelsperp=TRUE, horizontal=TRUE)
}
\keyword{hplot}

