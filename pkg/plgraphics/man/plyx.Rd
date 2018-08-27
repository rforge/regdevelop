\name{plyx}
\alias{plyx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scatterplot, enhanced
}
\description{A scatterplot or a bunch of them is produced according to
  the concept of the  plplot  package
}
\usage{
plyx(x, y = NULL, group = NULL, data = NULL, type = "p",
  xlab = NULL, ylab = NULL, markextremes = 0, rescale = TRUE,
  mar = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{either a formula or the data to be used for the horizontal
    axis. If a formula of the type 'y~x', the variable 'y' in
    'data' will be plotted against the variable(s) 'x'.
    If a data.frame with more than one column is given,
    each column will be used in turn to produce a plot.
}
  \item{y}{data to be used as the y axis.
  }
  \item{group}{grouping factor: for each group, a plot will be shown
    for the respective subset of the data}
  \item{data}{data.frame containing the variables if 'x' is a formula
  }
  \item{xlab, ylab}{axis labels}
  \item{type}{type of plot, see \code{?\link{plot.default}}}
%%-  \item{linecol}{colors to be used for the y variables. This is useful in
%%-    case there is more than one.
%%-    If \code{NULL}, \code{plyx} will use the \code{linecol} attribute of
%%-    the variables in \code{data}. If this is not available,
%%-    \code{plargs$plappearance$default$linecol} will be used for
%%-    \code{type}s \code{"l", "b"} and
%%-    \code{plargs$plappearance$default$col} for other \code{type}s.
%%-  }
  \item{markextremes}{proportion of extreme residuals to be labeled.
    If all points should be labeled, let \code{markextremes=1}.}
  \item{rescale}{logical. Only applies if there are multiple y
    variables. If \code{TRUE}, the vertical axis will be adjusted
    for each of these variables.
  }
  \item{mar}{plot margins, see \code{\link{par}}
  }
  \item{\dots}{ more arguments, to be passed to \code{\link{pl.control}}
  }
}
%% \details{
%% 
%% }
\value{None.}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\note{There are many more arguments, obtained from \code{pl.control},
  see \code{?\link{pl.control}}. These can be passed to \code{plmatrix}
  by an argument \code{plargs} that is hidden in the ... argument list.
}

\author{Werner A. Stahel, ETH Zurich}

\seealso{
\code{\link{plmatrix}}
}
\examples{
plyx(Petal.Width ~ Sepal.Length, data=iris)
plyx(Petal.Width ~ Sepal.Length+Sepal.Width, data=iris, smooth=TRUE,
     smooth.group=Species)
plyx(Petal.Length + Petal.Width ~ Sepal.Length+Sepal.Width,
     group = Species, data=iris, smooth=TRUE)

}
\keyword{hplot}

