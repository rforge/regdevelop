\name{plpanel}
\alias{plpanel}
\alias{panelSmooth}
\alias{plpanelCond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Panel function for multiple plots}
\description{
  Draw a scatterplot or multibox plot, usuallly after \code{pl.control}
  and \code{plframe} have been called. 
  May also be used to augment an existing plot.
}
\usage{
plpanel(x, y, indx = NULL, indy = NULL, type = "p", frame = FALSE,
  title = FALSE, plargs = NULL, ploptions = NULL, ...)

panelSmooth(x, y, indx, indy, plargs = NULL, ...)

plpanelCond(x, y, ckeyx, ckeyy, pch = 1, pcol = 1, psize = 1,
  pale = c(0.2, 0.6), csize = 0.8,
  smooth = NULL, smooth.minobs = NULL, plargs = NULL, ...)
}

\arguments{
  \item{x}{values of the horizontal variable}
  \item{y}{values of the vertical variable}
  \item{indx}{index of the variable shown horizontally, among the
    \code{y} variables}
  \item{indy}{index of the variable shown horizontally, among the
    \code{y} variables}
  \item{type}{type of plot as usual in R: "p" for points, ...}
  \item{frame}{logical: should \code{plframe} be called?
  }
  \item{title}{logical: should \code{pltitle} be called?
  }
  \item{ckeyx, ckeyy}{vectors of 'keys' to calculate paling values and
    weights for smoothing. NA means that points should not be shown
    in this panel. 0 means no paling and weight 1.
    Other values are between -1 and 1,
    \code{cpl=(1-abs(ckeyx))*(1-abs(ckeyy))} is used for paling and weights.
  }
  \item{pch, pcol, psize}{vector of plotting symbols, colors and sizes
    for plotting points
  }
  \item{pale}{vector of length 2 indicating the range of paling values
    obtained from \code{cpl} values from 1 to 0.
  }
  \item{csize}{factor applied to the character expansion of the points
    with \code{cpl<1}
  }
  \item{smooth}{should a smooth line be drawn?
  }
  \item{smooth.minobs}{minimum number of points required for calculating
    and showing a smooth line
  }
  \item{plargs, ploptions}{result of calling \code{pl.control}.
    If \code{plargs} is \code{NULL}, \code{pl.control} will be called to
    generate it.
    The components are often needed to generate the panel.
  }
  \item{...}{further arguments passed to
    \code{plpoints, plmboxes, plsmooth}
  }
}
\details{
  The panel function \code{plpanel} draws a scatterplot if both
  \code{x} and \code{y} are
  numerical, and a multibox plot if one of them is a factor and
  \code{ploptions$factor.show == "mbox"}.\cr
  Grouping, reference and smooth lines and properties of the points
  are determined by the component of \code{plargs} in \code{plpanel}.

  This function is usually called by the high level pl functions
  \code{\link{plyx}} and \code{\link{plmatrix}}.
  A different suitable function can be used by setting their
  argument \code{panel}.

  The first arguments, \code{x} and \code{y},
  can be formulas, and an argument \code{data} can be given.
  These arguments then have the same meaning as in \code{\link{plyx}},
  with the restriction that only one variable should result for
  the \code{x} and \code{y} coordinates in the plot.
  When \code{frame} is true, \code{plpanel} can be used instead of
  \code{\link{plyx}} for generating a single plot.
  Note that \code{plpanel} does not set \code{.plargs} in the
  global environment, in contrast to \code{plyx}.

  \code{plpanelCond} shows selected points only and may show
  some of them with reduced size and paled color.
  It is appropriate for the high level function \code{\link{plcond}}.
}
\value{
  none
}
\note{These functions are rarely called by the user.
  The intention is to modify ond of them and then call the modified
  version when using \code{\link{plyx}, \link{plmatrix}} or
  \code{\link{plcond}} by setting \code{panel=mypanel}.
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{plyx}} is essentially a wrapper function of
  \code{plpanel} which calls \code{pl.control} and provides additional
  features.
  \code{\link{plmatrix}} also uses \code{plpanel}, whereas
  \code{\link{plcond}} uses \code{plpanelCond}.
}
\examples{
t.plargs <-
  pl.control(~Species+Petal.Length, ~Sepal.Width+Sepal.Length,
             data=iris, smooth.group=Species, pcol=Species)
t.plargs$ploptions$group.col <- c("magenta","orange","cyan")
plpanel(iris$Petal.Length, iris$Petal.Width, plargs=t.plargs,
        frame=TRUE)
}
\keyword{hplot}
