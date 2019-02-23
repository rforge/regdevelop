\name{residuals.regr}
\alias{residuals.regr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residuals for \code{regr} objects
}
\description{Calculates Residuals for \code{regr} objects.
  This is a kind of generic function which calls further mehtods
  according to the secondary class of the object
}
\usage{
\method{residuals}{regr}(object, type = NULL, standardized = FALSE, ...)
}

\arguments{
  \item{object}{an object of class \code{regr}
  }
  \item{type}{type of residual, used by some methods that will be called
  }
  \item{standardized}{logical: Should standardized residuals be
    calculated and returned?
  }
  \item{\dots}{argiments passed to the methods that will be called
}
}
\value{vector or matrix of residuals
}
\author{Werner A. Stahel
}
\seealso{
  \code{\link{residuals.regrpolr}}, \code{\link{residuals.regrcoxph}},
  \code{\link{residuals.regrsurvreg}}
}
\examples{
data(d.blast, package="plgraphics")
r.blast <-
  regr(logst(tremor)~location+log10(distance)+log10(charge), data=d.blast)
str(residuals(r.blast))

data(housing, package="MASS")
t.r <- regr(Sat ~ Infl + Type + Cont, weights = housing$Freq,
            data = housing)
showd(residuals(t.r))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
