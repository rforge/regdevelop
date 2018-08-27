\name{genvarattributes}
\alias{genvarattributes}
\alias{varattributes}
\title{Generate Variable Attributes for Plotting
}
\description{
  \code{genvarattributes} generates attributes of variables that are
  useful for the \code{plgraphics} functions.
  It is called by \code{\link{pl.control}}.\cr
  \code{varattributes} modifies or sets such attributes.
}
\usage{
genvarattributes(data, ynames = NULL, ycol = NULL, ylty = NULL, ypch = NULL,
  labels = NULL, innerrange.limits = NULL, ploptions = NULL, replace=FALSE)

varattributes(data, attributes = NULL)}

\arguments{
  \item{data}{data.frame consisting of the variables (columns) to be
    characterized by their attributes
  }
  \item{ynames}{names of variables to be treated as y variables}
  \item{ycol, ylty, ypch}{color, line type and plotting character
    to be used when multiple y-s are plotted (in the sense of
    \code{matplot})
  }
  \item{labels}{labels of the variables, in the case that the
    \code{names} of \code{data} are not appropriate.
  }
  \item{innerrange.limits}{limits of the inner plotting range, if
    predetermined} 
  \item{ploptions}{list containing the plotting elements
    needed to set the attributes}
  \item{replace}{logical: should existing attributes be replaced?}
  \item{attributes}{(for \code{varattributes}) is a list of lists.
    Its names identify the variables for which the attributes are
    set or modified. Each component is a list which is added to the
    existing attributes of the respective variable or replaces them if
    they already exist.}
}
\details{
  If the attribute \code{innerrange} is replaced, then \code{plcoord} is
  also replaced.
}
\value{
  Data.frame, returning the original values, but the variables are
  supplemented by the following \code{attributes}, where available:
  \item{nvalues}{number of distinct values
  }
  \item{innerrange}{inner plotting range
  }
  \item{plcoord}{plotting coordinates
  }
  \item{axisat}{tick marks for axis
  }
  \item{label}{label to be used as axis label
  }
}
\author{Werner A. Stahel
}

\seealso{
\code{\link{par}}
}
\examples{
data(d.blast)
dd <- genvarattributes(d.blast)
str(attributes(dd$tremor))

ddd <- varattributes(dd, list( tremor=list(axisat=seq(0,24,2),
  axislabels = seq(0,24,10)) ) )
str(attr(ddd$tremor, "axislabels"))
}

\keyword{utilities}

