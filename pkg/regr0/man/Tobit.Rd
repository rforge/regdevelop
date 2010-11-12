\name{Tobit}
\alias{Tobit}

\title{
Prepare a Response for a Tobit Model
}
\description{
Returns a \code{Surv} object that allows for setting up a Tobit
regression model by calling \code{survreg} 
}
\usage{
Tobit(data, limit = 0, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    the variable to be used as the response in the Tobit regression
}
  \item{limit}{
    Lower limit which censors the observations.
    If \code{log} is \code{TRUE}, then the default is the minimal value
    of \code{\link{logst}}(data), and if \code{limit>0}, it refers to the
    untransformed data.
}
  \item{log}{
    logical. If \code{TRUE}, \code{data} will be log transformed
    by calling \code{\link{logst}}.
  }
}
\details{
Tobit regression is a special case of regression with left censored
response data. The function \code{survreg} is suitable for fitting.
In \code{regr}, this is done automatically.
}
\value{
A \code{Surv} object.
}
\author{
Werner A. Stahel
}

\seealso{
  \code{Surv}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(data, limit=NULL, log=FALSE)
{
  ## Purpose:   create a Surv object for tobit regression
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  1 Jan 2010, 21:49
  if (is.null(limit)) limit <- min(data)
  if (log) {
    if (limit<=0) {
      data <- logst(data)
      limit <- min(data)
    } else {
      data <- log10(data)
      limit <- log10(limit)
    }
  }
  data <- pmax(data,limit)
  if (sum(data==limit)<=1)
    warning(":Tobit: <= 1 observation equal to 'limit'")
  rr <- Surv(data, data!=limit, type='left')
  attr(rr,'distribution') <- 'gaussian'
  class(rr) <- c(class(rr), "matrix")
  rr
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
