\name{genvarattributes}
\alias{genvarattributes}
\alias{setvarattributes}
\title{Generate or Set Variable Attributes for Plotting
}
\description{
  \code{genvarattributes} generates attributes of variables that are
  useful for the \code{plgraphics} functions.
  It is called by \code{\link{pl.control}}.\cr
  \code{setvarattributes} modifies or sets such attributes.
}
\usage{
genvarattributes(data, vnames = NULL, vcol = NULL, vlty = NULL, vpch = NULL,
  varlabels = NULL, innerrange.limits = NULL, replace=FALSE, ploptions = NULL)

setvarattributes(data, attributes = NULL, list = NULL, ...)
}

\arguments{
  \item{data}{data.frame consisting of the variables (columns) to be
    characterized by their attributes
  }
  \item{vnames}{names of variables to be treated as y variables}
  \item{vcol, vlty, vpch}{color, line type and plotting character
    to be used when multiple y-s are plotted (in the sense of
    \code{matplot})
  }
  \item{varlabels}{labels of the variables, in the case that the
    \code{names} of \code{data} are not appropriate.
  }
  \item{innerrange.limits}{limits of the inner plotting range, if
    predetermined} 
  \item{ploptions}{list containing the plotting elements
    needed to set the attributes}
  \item{replace}{logical: should existing attributes be replaced?}
  \item{attributes}{(for \code{setvarattributes}) is a list of lists.
    Its names identify the variables for which the attributes are
    set or modified. Each component is a list which is added to the
    existing attributes of the respective variable or replaces them if
    they already exist.}
  \item{list}{a list of attributes to be set.
    Each component must have a name giving the name of the variable
    attribute to be set, and be itself a list (or a vector).
    This list must have names that identify the variables in 
    \code{data} for which the attributes are set.
    See examples to understand this.}
  \item{...}{further arguments, which will be collected and used as or
    added to\code{list}}
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
  \item{ticksat}{tick marks for axis
  }
  \item{varlabel}{label to be used as axis label
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

ddd <- setvarattributes(dd, list( tremor=list(ticksat=seq(0,24,2),
  ticklabelsat = seq(0,24,10), ticklabels=c("low","medium","high")) ) )
str(attributes(ddd$tremor))

data(d.river)
plyx(O2+H2CO3+T ~ date, data=d.river, subset=as.Date(date)<as.Date("2010-02-28"))
dd <- setvarattributes(d.river,
  list=list(vcol=c(O2="blue", T="red")), vpch=c(O2=1, T="T", H2CO3=5) )
attributes(dd$O2)
plyx(O2+H2CO3+T ~ date, data=d.river, subset=as.Date(date)<as.Date("2010-02-28"))
}

\keyword{utilities}

