\name{termeffects}
\alias{termeffects}
\alias{print.termeffects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
all coefficients of a model fit
}
\description{
A list of all coefficients of a model fit, possibly with respective statistics
}
\usage{
termeffects(object, se = 2, df = df.residual(object), ...) %% use.na = TRUE, 
\method{print}{termeffects}(x, columns=NULL, transpose=FALSE,
  single = FALSE...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a model fit, produced, e.g., by a call to \code{lm} or
    \code{regr}
}
%%-   \item{use.na}{as in \code{\link{dummy.coef}}
%%- }
  \item{se}{logical: Should inference statistics be generated?
}
  \item{df}{degrees of freedom for t-test
  }
  
  Arguments for \code{print}:
  \item{x}{object of class \code{termeffects}, produced by \code{termeffects(...)}}
  \item{columns}{the columns of the tables to be printed. \code{NULL}
  selects all the columns.}
  \item{transpose}{logical. should the tables be
    transposed for printing?}
  \item{single}{logical. should components with a single coefficient be printed?}
  \item{\dots}{absorbs further arguments of call to generic function
}
}
\value{list with a component for each term in the model formula.
  The components are either vectors of coefficients
  (if \code{se} is false, else data.frames with columns
  \item{coef}{}
  
}
\author{
Werner A. Stahel
}
\seealso{
  \code{\link{dummy.coef}}
}
\examples{
  data(d.blast, package="plgraphics")
  r.blast <-
  regr(logst(tremor)~location+log10(distance)+log10(charge), data=d.blast)
  termeffects(r.blast)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
