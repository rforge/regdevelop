\name{plyx}
\alias{plyx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scatterplot, enhanced
}
\description{A scatterplot or a bunch of them is produced according to
  the concept of the  plplot  package
}
\usage{
plyx(x = NULL, y = NULL, group = NULL, data = NULL, type = "p",
  panel = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
  markextremes = 0, rescale = TRUE,
  mar = NULL, plargs = NULL, ploptions = NULL, ...)
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
  \item{xlim, ylim}{plot ranges}
  \item{type}{type of plot, see \code{?\link{plot.default}}}
  \item{panel}{panel function to do the actual drawing. See Details.
  }
  \item{markextremes}{proportion of extreme residuals to be labeled.
    If all points should be labeled, let \code{markextremes=1}.}
  \item{rescale}{logical. Only applies if there are multiple y
    variables. If \code{TRUE}, the vertical axis will be adjusted
    for each of these variables.
  }
  \item{mar}{plot margins, see \code{\link{par}}
  }
  \item{plargs}{result of calling \code{pl.control}.
    If \code{NULL}, \code{pl.control} will be called to generate it.
    If not null, arguments given in \code{...} will be ignored.
  }
  \item{ploptions}{list of pl options.
  }
  \item{\dots}{ more arguments, to be passed to \code{\link{pl.control}}
  }
}
\details{
    Defaults to \code{plpanel}, which results essentially in
    \code{\link{points}} or \code{\link{text}}
    depending on the argument \code{pch} including a smooth line,
    to \code{\link{plmboxes}} if  'x'  is a factor and  'y'  is not or
    vice versa,
    or to a modification of \code{sunflowers} if both are factors.\cr
    The function must have the arguments \code{x} and \code{y}
    to take the coordinates of the points and may have the arguments
    \code{indx} and \code{indy} to transfer the variables\' index and
    \code{panelargs} for any additional objects to be passed on.
}
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

