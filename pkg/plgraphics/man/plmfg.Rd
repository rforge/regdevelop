\name{plmfg}
\alias{plmfg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple Frames for Plotting}
\description{
  This is a short-cut to set some graphical parameters
}
\usage{
plmfg(mfrow = NULL, mfcol = NULL, mft = NULL, nrow = NULL, ncol = NULL,
  row = TRUE,
  oma = NULL, mar = NULL, mgp = NULL, \dots)
}
\arguments{
  \item{mfrow, mfcol}{number of rows and columns of panels. The default
    is 1 for both, which will reset the subdivision of the plotting page.}
  \item{mft}{total number of panels, to be split into \code{mfrow}
    and \code{mfcol} by the function. The result depends on the current
    aspect ratio (ratio of height to width) of the plotting area.}
  \item{nrow, ncol}{maximum number of rows and columns on one page}
  \item{row}{if TRUE, the panels will be used by rows, otherwise, by
    columns}
  \item{mar}{plot margins.
    Any \code{NA}s in \code{mar} will be replaced by the respective
    current values.
  }
  \item{oma}{outer plot margins. They will be widened if necessary to
    accomodate axis annotations.
  }
  \item{mgp, \dots}{further graphical parameters passed to
    \code{\link{par}(...)}.
  }
}
\details{
  The function calls \code{par}. Its purpose is to simplify a call like
  \code{par(mfrow=c(3,4))} to \code{plmfg(3,4)} and to set some
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
plmfg(2,3)
plmfg(mft=15)  ## will split the plotting area into >= 15 panels,
plmfg()  ## reset to 1 panel
}
\keyword{utilities}
