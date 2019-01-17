\name{regr.control}
\alias{regr.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Further Arguments to \code{regr}
}
\description{
  Specify some arguments of minor importance for the function
  \code{\link{regr}}
}
\usage{
regr.control(contrasts = getUserOption("regr.contrasts"),
  factorNA = TRUE, na.action = as.name("nainf.exclude"),
  calcdisp = NULL, suffmean = 3, dist = NULL,
  model = FALSE, x = TRUE, termtable = TRUE, vif = TRUE,
  testlevel = 0.05, leveragelim = c(0.99, 0.5), tit = NULL,
  control = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contrasts}{
%%     ~~Describe \code{contrasts} here~~
}
  \item{factorNA}{
%%     ~~Describe \code{factorNA} here~~
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
}
  \item{calcdisp}{
%%     ~~Describe \code{calcdisp} here~~
}
  \item{suffmean}{
%%     ~~Describe \code{suffmean} here~~
}
  \item{dist}{
%%     ~~Describe \code{model} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{termtable}{
%%     ~~Describe \code{termtable} here~~
}
  \item{vif}{
%%     ~~Describe \code{vif} here~~
}
  \item{testlevel}{
%%     ~~Describe \code{testlevel} here~~
}
  \item{leveragelim}{
%%     ~~Describe \code{leveragelim} here~~
}
  \item{tit}{
%%     ~~Describe \code{tit} here~~
}
  \item{control}{parameters to be passed to specific fitting functions.
  }
}
\value{A list, containing the additional arguments needed by \code{regr}
}
\author{
WSt
}

\seealso{
  \code{\link{regr}}
}
\examples{
## This function is not usually called by a user.

}
\keyword{utilities}% use one of  RShowDoc("KEYWORDS")

