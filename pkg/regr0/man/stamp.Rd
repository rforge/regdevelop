\name{stamp}
\alias{stamp}
\title{ Add a "Stamp", i.e., an Identification Line to a Plot }
\description{
  A line is added to the current plot in the lower right corner
  that contains project information and date.
}
\usage{
stamp(sure = TRUE, outer.margin = NULL,
      project = getUserOption("project"), step=getUserOption("step"),
      stamp = getUserOption("stamp"), ...)
}
\arguments{
  \item{sure}{ if FALSE, the stamp will only be added if
    \code{getUserOption("stamp")>0} }
  \item{outer.margin}{ if TRUE, the stamp is put to the outer margin of
    the plot.  This is the default if the plot is currently split into panels.}
  \item{project, step}{character string describing the project and the step of
    analysis.}
  \item{stamp}{controls default action, see details }
  \item{\dots}{arguments passed to \code{mtext}}
}
\details{
  The function is used to document plots produced during a data
  analysis. It is called by all plotting functions of this package.
  For getting final presentation versions of the plots, the stamp can be
  suppressed by changing the default by calling \code{options(stamp=0)}.

  In more detail:  If \code{stamp==0} (or \code{\link{getUserOption}("stamp")==0})
  the function will only do its thing if \code{sure==TRUE}.

  If \code{stamp==2}, it will certainly do it.

  If \code{stamp==1} and \code{sure==FALSE}, the stamp is added when a
  plot page is complete.
}
\value{
  invisibly returns the string that is added to the plot
}
\author{Werner A. Stahel, ETH Zurich}
\examples{
userOptions(project="Example A",  step="regression analysis")
plot(1:10)
stamp() #-> "stamp" at bottom of right border
}
\keyword{utilities}
\keyword{misc}
