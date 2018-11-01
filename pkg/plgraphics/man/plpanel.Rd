\name{plpanel}
\alias{plpanel}
\alias{panelSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Panel function for multiple plots}
\description{
  Draw a scatterplot or multibox plot after \code{pl.control} has been
  called. 
}
\usage{
plpanel(x, y, indx = NULL, indy = NULL, type = "p", plargs = NULL,
  frame = FALSE, title = FALSE, ...)

panelSmooth(x, y, indx, indy, plargs, ...)
}

\arguments{
  \item{x}{values of the horizontal variable}
  \item{y}{values of the vertical variable}
  \item{indx}{index of the variable shown horizontally, among the
    \code{y} variables}
  \item{indy}{index of the variable shown horizontally, among the
    \code{y} variables}
  \item{type}{type of plot as usual in R: "p" for points, ...}
  \item{plargs}{result of calling \code{pl.control}.
    If \code{NULL}, \code{pl.control} will be called to generate it.
    Its components are often needed to generate the panel.
  }
  \item{frame}{logical: should \code{plframe} be called?
  }
  \item{title}{logical: should \code{pltitle} be called?
  }
  \item{...}{further arguments passed to \code{plmboxes}
  }
}
\details{
  The panel function draws a scatterplot if both \code{x} and \code{y} are
  numerical, and a multibox plot if one of them is a factor and
  \code{plargs$ploptions$factor.show == "mbox"}.\cr
  Grouping, reference and smooth lines and properties of the points
  are determined by the component of \code{plargs} in \code{plpanel}.
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
