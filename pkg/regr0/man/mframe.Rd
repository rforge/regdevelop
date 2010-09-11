\name{mframe}
\alias{mframe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple Frames for Plotting}
\description{
  This is a short-cut to set some graphical parameters
}
\usage{
mframe(mfrow = NULL, mfcol = NULL, mft = NULL, row = TRUE,
  oma = c(0, 0, 2, 1), mar = getOption("mar"),
  mgp = getOption("mgp"), \dots)
}
\arguments{
  \item{mfrow, mfcol}{number of rows and columns of panels. The default
    is 1 for both, which will reset the subdivision of the plotting page.}
  \item{mft}{total number of panels, to be split into \code{mfrow}
    and \code{mfcol} by the function. The result depends on the current
    aspect ratio (ratio of height to width) of the plotting area.}
  \item{row}{if TRUE, the panels will be used by rows, otherwise, by columns}
  \item{oma, mar, mgp, \dots}{further graphical parameters passed to
    \code{\link{par}(...)}.}
}
\details{
  The function calls \code{par}. Its purpose is to simplify a call like
  \code{par(mfrow=c(3,4))} to \code{mframe(3,4)} and to set some
  defaults differently from \code{par}.
}
\value{
  a named \code{\link{list}} containing the old values of the
  parameters, as for \code{par}.
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{par}}; \code{\link[sfsmisc]{mult.fig}} from package
  \pkg{sfsmisc} with very similar goals.
}
\examples{
mframe(2,3)
mframe(mft=15)  ## will split the plotting area into >= 15 panels,
mframe()  ## reset to 1 panel
}
\keyword{utilities}
