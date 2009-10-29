\name{lassogrp}
\Rdversion{1.1}
\alias{lassogrp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lassogrp(x, ...)

lassogrp.formula(formula, data, subset, weights = rep(1, length(y)), na.action, model = "gaussian", nonpen = ~1, lambda = NULL, lfac = 2^seq(-1, -10), contrasts = NULL, ...)

lassogrp.default(x, y, index, subset, weights = rep(1, length(y)), model = "gaussian", offset = rep(0, length(y)), lambda = NULL, lstep = 21, coef.init = rep(0, ncol(x)), penscale = sqrt, center = NULL, standardize = TRUE, save.x = NULL, save.y = NULL, control = lassoControl(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, ...)
  UseMethod("lassogrp")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
