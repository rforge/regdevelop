\name{pl.control}
\alias{pl.control}
\title{Arguments for plotting functions}
\description{Arguments that can be specified calling \code{plyx} and
  other 'pl' functions are checked and data is prepared for plotting.  
}
\usage{
pl.control(x=NULL, y=NULL, data = NULL, transformed = TRUE,
           ## generate variables used for plotting
           psize = NULL, plab = FALSE, pch = NULL, pcol = NULL,
           cex = NULL, markextremes = NULL,
           ycol = NULL, ylty = NULL, ypch = NULL, smooth = NULL,
           main = NULL, sub = ":", .subdefault = NULL,
           mar = NULL, xlab = NULL, ylab = NULL,
           varlabels = NULL,
           ploptions = NULL, .environment. = parent.frame(), ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y, data}{as in \code{\link{plyx}}
  }
  \item{transformed}{logical: should transformed variables be used?
  }
  \item{psize, plab, pch, pcol}{Plotting characteristics of points,
    specified as a (unquoted) variable name found in \code{data}
    or as a vector.
    They set the size of the plotting symbols, labels (character
    strings), plotting character, and color, respectively.\cr
    \code{plabs = TRUE} asks for using the row names of \code{data}.
  }
  \item{cex}{character expansion, applied to both labels and plotting
    characters.
    By default, it adjusts to the number of observations.
  }
  \item{markextremes}{scalar: proportion of extreme points to be labelled
  }
  \item{ycol, ylty, ypch}{color, line type and plotting character
    to be used when multiple y-s are plotted (in the sense of
    \code{matplot})
  }
  \item{smooth}{logical: should a smooth line be added?}
  \item{main, sub}{string. Main title of the plot(s).
    If \code{sub} starts by \code{":"} (the default),
    \code{pl.control} tries to generate an informative subtitle,
    determined by the data or a model formula.
  }
  \item{.subdefault}{for internal use: default of subtitle}
  \item{mar}{plot margins}
  \item{xlab, ylab}{axis labels}
  \item{varlabels}{
    labels for variables replacing their names in the \code{x} and
    \code{y} arguments, either a simple vector of strings with an
    element for each variable,
    or a named vector, where names correspond to such variables.
  }
  \item{ploptions}{Plotting attributes, e.g., plotting character,
    line types, colors and the like, for different aspects of plots.
    Result of \code{\link{ploptions}}. Defaults to \code{.ploptions}.
  }
  \item{.environment.}{used by the calling function to provide the
    environment for evaluating \code{x} and \code{y}
  }
  \item{\dots}{further arguments, !!!
}
}
\details{
!!!
}
\value{A list containing all the arguments, possibly in modified form.
  Specifically, the evaluations of the variables contained in
  \code{x} and \code{y} along with 
  \code{psize, plab, pch, pcol, smoothGroup, smoothWeights}
  are collected in the component \code{pdata}.
  The component, \code{plappearances}, contains the result of calling
  the respective function or the corresponding component of the
  user options.
}
\author{Werner A. Stahel
}
\seealso{
  \code{\link{plyx}}, \code{\link{plmatrix}}, \code{\link{ploptions}}
}
\examples{
plyx(Sepal.Width~Sepal.Length, data=iris, axp=7, plab=TRUE, cex.plab=0.6)
## same as
plargs <- pl.control(Sepal.Width~Sepal.Length, data=iris)
plargs$pdata$plab <- row.names(iris)
plargs$cex.lab <- 0.6
plargs$axp <- 7
plyx(Sepal.Width~Sepal.Length, plargs=plargs)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
