\name{plcond}
\alias{plcond}

\title{Plot Two Variables Conditional on Two Others
}
\description{A scatterplot matrix is generated that shows, in each
  panel, the relationship between two primary variables, with the
  dataset restricted by appropriate subranges of two 'conditioning'
  variables.
  This corresponds to \code{link{coplot}}.
  The points that are near to the the 'window' defining the panel's
  restriction are also shown, in a distinct style.
}
\usage{
plcond(x, y = NULL, condvar = NULL, data = NULL,
  panel = NULL, nrow = NULL, ncol = NULL,
  xaxmar = NULL, yaxmar = NULL, xlab = NULL, ylab = NULL,
  oma = NULL, plargs = NULL, ploptions = NULL, assign = TRUE, ...)
}

\arguments{
  \item{x, y}{the two variables used to generate each panel.
    They may be specified as vectors, as column names of \code{data} or
    by formulas as in \code{\link{plyx}}.
  }
  \item{condvar}{two (or one) variables that define the restrictions of
    the data for the different panels.
    A numerical variable is cut into intervals, see Details.
    A factor defines the 'ranges' as its levels. For each combination
    of intervals or levels of the two variables, a panel is generated.
  }
  \item{data}{data.frame in which the variables are found if needed
  }
  \item{panel}{function that generates each panel.
    If set by the user, it must accept the arguments
    \code{x, y, ckeyx, ckeyy, pcol, pale, cex, smooth, smooth.minobs,
      ploptions}. The default is \code{ploptions("plcond.panel")},
      which in turn is initiated as the function \code{plpanelCond}.
  }
  \item{nrow, ncol}{number of maximum rows and columns on a page
  }
  \item{xaxmar, yaxmar}{margin in which the axis (tick marks and
    corresponding labels) should be shown: either 1 or 3 for
    \code{xaxmar} and 2 or 4 for \code{yaxmar}.
  }
  \item{xlab, ylab}{labels of the variables \code{x} and \code{y}
  }
  \item{oma}{width of outer margins, see \code{\link{par}}.
    Note that a minimum of 2.1 is generally needed for showing
    tick and axis labels.}
%%  \item{cex}{general character expansion, will be multiplied by
%%    \code{par("cex")}
%%  }
  \item{plargs}{result of calling \code{pl.control}.
    If \code{NULL}, \code{pl.control} will be called to generate it.
    If not null, arguments given in \code{...} will be ignored.
  }
  \item{ploptions}{list of pl options.
  }
  \item{assign}{logical: Should the plargs be stored as \code{.plargs}
    in the \code{.GlobelEnv}?
  }
  \item{\dots}{further arguments passed to the \code{panel} function
    and possibly further to functions called by the panel function\cr
  }
}
\details{A numerical conditioning variable (\code{condvar}) will be
  split by default into classes by splitting its robust range
  (\code{\link{robrange}}) into \code{ploptions("plcond.nintervals")}
  equally long intervals. Alternatively, the variable may contain
  an attribute \code{cutpoints} which then defines the intervals.

  For numerical conditioning variables, each panel also shows
  neighboring points with a different color and diminished size.
  The size of the neighborhood is defined by the proportion of extension 
  \code{ploptions("plcond.ext")}.
  The point size of the respective 'exterior' points is given by
  \code{ploptions("plcond.cex")} The color are given by
  the 4 elements of \code{ploptions("plcond.col")}:
  The first element is used to paint the neighboring points
  to the left of the current range of the conditioning x variable, 
  the second element paints those to the right,
  and the third and fourth are used in the same way for the
  conditioning y variable. The neighboring points that are outside
  both ranges get a color mixing the two applicable colors according to
  this rule.
  Finally, paling is applied to these colors with a degree that is
  linear in the distance from the interval, determined by the range
  given by \code{ploptions("plcond.pale")}.
}
\value{None.
}
\author{Werner A. Stahel
}

\seealso{\code{\link{coplot}}
}
\examples{
plcond(Sepal.Width~Sepal.Length, data=iris, condvar=~Species+Petal.Length)
}
\keyword{plot}

