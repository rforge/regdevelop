\name{plpanel}
\alias{plpanel}
\alias{panelSmooth}
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

panelSmooth(x, y, indx, indy, plargs = NULL, ploptions = NULL, ...)
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
  \item{plargs, ploptions}{result of calling \code{pl.control}.
    If \code{plargs} is \code{NULL}, \code{pl.control} will be called to
    generate it.
    The components are often needed to generate the panel.
  }
  \item{...}{further arguments passed to \code{plmboxes}
  }
}
\details{
  The panel function draws a scatterplot if both \code{x} and \code{y} are
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
}
\value{
  none
}
\note{This function is rarely called by the user.
  The intention is to modify it and then call the modified version
  when using \code{\link{plyx}} or \code{\link{plmatrix}} by setting
  \code{panel=mypanel}.
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{plyx}} is essentially a wrapper function of
  \code{plpanel} which calls \code{pl.control} and provides additional
  features.
}
\examples{
t.plargs <- pl.control(~Species+Petal.Length, ~Sepal.Width+Sepal.Length,
                       data=iris, smooth.group=Species, pcol=Species)
t.plargs$ploptions$group.col <- c("magenta","orange","cyan")
plpanel(iris$Petal.Length, iris$Petal.Width, plargs=t.plargs, frame=TRUE)
}
\keyword{hplot}
