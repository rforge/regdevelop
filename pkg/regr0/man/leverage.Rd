\name{leverage}
\alias{leverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get leverage values}
\description{
  Extracts the leverage component of a fit object using the
  \code{na.action} component if available
}
\usage{
leverage(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
    an object containing a component \code{fit$leverage} and possibly
    a component \code{fit$na.action}
  }
}
\value{
  The vector \code{fit$leverage}, possibly expanded by missing values
  if \code{fit$na.action} has class \code{na.replace}
}
\author{Werner A. Stahel, ETH Zurich}

\seealso{\code{\link{influence}}}

\examples{
data(d.blast)
r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
t.lev <- leverage(r.blast)
length(t.lev)
dim(d.blast) ## number of observations is 388
length(r.blast$leverage) ## 362
  ## only observations without missings occur in $leverage
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
% \keyword{diagnostics}
