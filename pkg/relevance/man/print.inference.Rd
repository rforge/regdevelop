\name{print.inference}
\alias{print.inference}
\alias{print.termeffects}
\alias{print.printList}

\title{Print Tables with Inference Measures}
\description{Print methods for objects of class
  \code{inference, termtable}, or \code{termeffects}
}
\usage{
\method{print}{inference}(x, show = getOption("show.inference"), print=TRUE,
  digits = getOption("digits.reduced"), transpose.ok = TRUE,
  legend = NULL, na.print = "  ", ...)

\method{print}{termeffects}(x, show = getOption("show.inference"),
  transpose.ok = TRUE, single = FALSE, print = TRUE, warn = TRUE, ...)

\method{print}{printList}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object to be printed}
  \item{show}{determines items (columns) to be shown}
  \item{digits}{number of significant digits to be printed}
  \item{transpose.ok}{logical: May a single column be shown as a row?}
  \item{single}{logical: Should components with a single coefficient
    be printed?}
  \item{legend}{logical: should the legend(s) for the symbols
    characterizing p-values and relevances be printed?
    Defaults to \code{regroptions("show.symbolLegend")}.
  }
  \item{na.print}{string by which \code{NA}s are shown}
  \item{print}{logical: if \code{FALSE}, no printing will occur,
    used to edit the result before printing it.}
  \item{warn}{logical: Should the warning be issued if
    \code{termeffects} has nothing to print since there are no
    terms with more than one degree of freedom}    
  \item{...}{further arguments, passed to \code{print.data.frame}}
}
\value{A kind of formatted version of \code{x}, with class
  \code{printList}.
  For \code{print.inference}, it will be
  a character vector or a data.frame with attributes
  \code{head} and \code{tail} if applicable.
  For \code{print.termeffects}, it will be a list of such elements,
  with its own \code{head} and \code{tail}.
  It is invisibly returned.
}
\details{The result, say \code{rr}, can be printed by using
  \code{print.printList} writing just \code{print(rr)}, which is just
  what happens internally unless \code{print=FALSE} is used.  
  This allows for editing the result before printing it, see Examples
  }
\author{Werner A. Stahel}

\seealso{
  \code{\link{twosamples}, \link{termtable}, \link{termeffects}, \link{inference}}
}
\examples{
data(d.blast)
r.blast <-
  lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
rt <- termtable(r.blast)
rt
print(rt, show="classical")
rr <- print(termeffects(r.blast), print=FALSE)
attr(rr, "head") <- sub("lm", "Linear Regression", attr(rr, "head"))
print(rr)
}
\keyword{regression}
\keyword{utility}
\concept{relevance}
