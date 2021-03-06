\name{plmframes}
\alias{plmframes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple Frames for Plotting}
\description{
  This is a short-cut to set some graphical parameters
}
\usage{
plmframes(mfrow = NULL, mfcol = NULL, mft = NULL, byrow = TRUE, reduce = FALSE, 
  oma = NULL, mar = NULL, mgp = NULL, plargs = NULL, ploptions = NULL, \dots)
}
\arguments{
  \item{mfrow, mfcol}{number of rows and columns of panels. The default
    is 1 for both, which will reset the subdivision of the plotting page.}
  \item{mft}{total number of panels, to be split into \code{mfrow}
    and \code{mfcol} by the function. The result depends on the current
    aspect ratio (ratio of height to width) of the plotting area.}
  \item{byrow}{if TRUE, the panels will be generated by rows, otherwise, by
    columns}
  \item{reduce}{logical: If the number of rows or columns asked for by
    \code{mfrow} or \code{mfcol} exceeds the maximum numbers determined
    from \code{ploptions("mframesmax")}, suitable numbers for multiple
    pages are calculated. If \code{reduce} is \code{TRUE}, these
    suggested numbers are applied.
  }
  \item{mar}{plot margins.
    Any \code{NA}s in \code{mar} will be replaced by appropriate values.
  }
  \item{oma}{outer plot margins.
    Any \code{NA}s will be replaced by appropriate values.
  }
  \item{mgp}{margin-pars passed to \code{\link{par}(...)}.
    If \code{NULL}, it will be generated.
  }
  \item{plargs, ploptions}{result of calling \code{pl.control}, used
    for generating appropriate values for the margin parameters.
    If \code{NULL}, \code{pl.envir} will be used.
  }
  \item{\dots}{further graphical parameters passed to
    \code{\link{par}(...)}.
  }
}
\details{
  The function calls \code{par}. Its purpose is to simplify a call like
  \code{par(mfrow=c(3,4))} to \code{plmframes(3,4)} and to set some
  defaults differently from \code{par}.
}
\value{
  A named \code{list} containing the old values of the parameters,
  as for \code{par}.
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{par}}
}
\examples{
plmframes(2,3)
plmframes(mft=15)  ## will split the plotting area into >= 15 panels,
plmframes()  ## reset to 1 panel

t.plo <- ploptions(mframesmax=9, assign=FALSE)
t.mf <- plmframes(4,4, reduce=TRUE, ploptions=t.plo)
par("mfg")  
t.mf[c("mfigsug","npages")]
##  $mfigsug
##  [1] 2 4
##  $npages
##  [1] 2 1
##  if the device area was higher than wide,
##  the result is the other way 'round
t.mft <- plmframes(mft=12, reduce=TRUE, ploptions=t.plo)

}
\keyword{utilities}
