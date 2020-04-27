\name{plmarginpar}
\alias{plmarginpar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set Graphical Parameters according to those used in the Pl
  Function Called Last
}
\description{
  \code{plmarginpar} calls \code{par} to set the margin widths
  \code{mar} and \code{mgp} equal to those used in the last call of a
  high level pl function
}
\usage{
plmarginpar(plargs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plargs}{list from which the margin parameters are obtained.
    If \code{NULL}, the default, \code{.plargs} is used. 
  }
}
\value{The old settings of \code{par(c("mar","mgp"))} are returned
  invisibly.
}
\note{\code{plmarginpar} is used to complement a plot with
  low level ordinary R functions like \code{mtext} or
  \code{segments}, see Example.

  The same effect can be achieved by setting the pl option
  \code{keeppar} to \code{TRUE}, either by calling \code{ploptions}
  or by setting \code{keeppar=TRUE} in the call to the high level
  pl function.
}

\author{Werner A. Stahel
}

\examples{
par(mar=c(2,2,5,2))
plyx(Sepal.Width~Sepal.Length, data=iris) ## margins according to ploptions
par("mar") ## paramteres have been recovered
mtext("wrong place for text",3,1, col="red")  ## margins not appropriate for active plot
plmarginpar()
par("mar") ## margins used inside the call to  plyx . These are now active
mtext("here is the right place",3,1, col="blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
