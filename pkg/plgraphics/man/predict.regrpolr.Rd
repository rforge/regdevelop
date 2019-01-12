\name{predict.regrpolr}
\alias{predict.regrpolr}
\alias{fitted.regrpolr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict and Fitted for polr Models}
\description{
Methods of \code{predict} and \code{fitted} 
}
\usage{
\method{predict}{regrpolr}(object, newdata = NULL,
  type =  c("class", "probs", "link"), ...)
\method{fitted}{regrpolr}(object, type = c("class", "probs", "link"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{result of \code{polr}}
  \item{newdata}{data frame in which to look for variables with
    which to predict.  If \code{NULL}, fitted values are produced.
  }
  \item{type}{type of prediction:\cr
    \code{"link"} asks for the linear predictor values.\cr
    Other \code{type}s are available according to the standard
    methods of the \code{predict} function.
  }
  \item{...}{arguments passed to standard methods of \code{predict}
    or \code{fitted}}
}
\value{Vector of predicted or linear predictor values
}
\author{Werner A. Stahel, ETH Zurich
}
\seealso{
  \code{\link{predict}, \link{fitted}, \link{residuals.regrpolr}} 
}
\examples{
require(MASS)
data(housing, package="MASS")
rr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
fitted(rr)
}
\keyword{regression}% use one of  RShowDoc("KEYWORDS")

