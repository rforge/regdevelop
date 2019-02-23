\name{createNAvars}
\alias{createNAvars}
\alias{factorNA}
\alias{xNA}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make NAs available for regression
}
\description{
  For numeric variables, introduces a new indicator (binary) variable
  for missing values and assigns a fixed, non-missing value for them 
  in the (renamed) original variable\cr
  For factors, introduces a new level for missing values and assigns
  it to all missing values.
}
\usage{
  createNAvars(data, vars=NULL, na.prop=0.1, na.label=".NA.",
    na.values = NULL, name.suffix=c(".X",".NA"), append=TRUE, ...)
  factorNA(data, na.label = ".NA.",na.prop=0,  ...)
  xNA(data, na.values = NULL, na.prop=0.1, name.suffix=c(".X",".NA"))
}

\arguments{
  \item{data}{a data.frame of quantitative variables. A vector will be
  converted to a one column data.frame
}
\item{vars}{variables (column names or indexes) to be treated.
  It is usually not needed to specify \code{vars} since only those
  variables are treated that contain missing values.
}
\item{na.prop}{limit for the proportion of NAs. Only those variables
  are treated that contain a higher proportion of missings.
}
  \item{na.label}{for factors, label to be used for the new level of the factor
}
\item{na.values}{for numeric variables, value used in the .X variable
  for the NAs
}
  \item{name.suffix}{suffix to be used for the new quantitative
  variables and the new indicator variables
}
\item{append}{logical: should the modified any newly created variables
  be returned with the unmodified variables?
  If \code{TRUE}, factors will be replaced by the modified version
  and the created variables based on continuous variables will be appended.}
\item{\dots}{other arguments, to be passed to the function \code{factor}
}
}
\value{
  data.frame, containing,
  \itemize{
    \item{
      for each numeric variable \code{Vx}, a variable
      \code{Vx.X} and a variable \code{Vx.NA} (if the default
      \code{name.suffix} is used),
    }
    \item{
      for each factor, the factor consisting of the original values,
      but with NAs replaced by a level '.NA.' (unless specified
      otherwise).
    }
  }
  If \code{append} is true, the data.frame will be the same as the
  input \code{data}, except that factors containing NAs will have been
  modified and the variables created for numeric variables containing
  NAs will be appended. Note that no variables will be appended
  for those that did not contain any NAs.
}
\author{
Werner A. Stahel
}

\seealso{
\code{\link{factorNA}}
}
\examples{
data(d.blast, package="plgraphics")
dd <- d.blast
dd$location[c(1,3,5)] <- NA
rr <- createNAvars(dd, na.prop=0)

xNA(c(3,NA,NA,6,5,9,1,11))

df <- data.frame(a=c(3,NA,NA,6,5,9,1,11), b=c(1,2,1,0,3,4,2,0))
( dfn <- xNA(df) )
attr(dfn, "xNA.values")

factorNA(c("a","c",NA,"a",NA), na.label="missing")
}
\keyword{manip}

