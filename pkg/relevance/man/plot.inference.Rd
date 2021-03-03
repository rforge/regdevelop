\name{plot.inference}
\alias{plot.inference}
\alias{plot.termeffects}

\title{Plot Inference Results
}
\description{Plot the relevance interval(s) for one or several items
}
\usage{
\method{plot}{inference}(x, pos = NULL,
  plpars=list(lwd=c(2,1,2), endmarks=1, extend=NA, framecol="gray70"),
  xlab = "relevance", ...)
\method{plot}{termeffects}(x, pos = NULL, single=FALSE,
  plpars=list(lwd=c(2,1,2), endmarks=1, extend=NA, framecol="gray70",
              termeffects.gap = 1.1), xlab="relevance", ...)
}

\arguments{
  \item{x}{a vector or matrix of class \code{inference}
  }
  \item{pos}{positions of the bars in vertical direction
  }
  \item{plpars}{graphical options for this function, see Details
  }
  \item{xlab}{label for horizontal axis
  }
%-  \item{sub}{subtitle, to be shown on top of the figure}
  \item{single}{logical: should terms with a single degree of freedom
    be plotted?}
  \item{\dots}{further arguments to the call of \code{plot}
  }
}
\details{The graphical options are:
  \itemize{
    \item{lwd}{line width for the intervals}
    \item{endmarks}{determines the length of the end marks of the
      intervals as a multiplier of the default length}
    \item{extend}{extension of the vertical axis beyond the range}
    \item{framecol}{color to be used for the framing lines:
      axis and vertical lines at relevances 0 and 1}
  }
}
\value{none
}
\author{Werner A. Stahel
}

%- \seealso{}
\examples{
data(swiss)
rr <- lm(Fertility ~ . , data = swiss)
rt <- termtable(rr)
plot(rt)

data(d.blast)
rlm <- lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
rte <- termeffects(rlm)
plot(rte, single=TRUE)
}
\keyword{hplot}
\concept{relevance}
