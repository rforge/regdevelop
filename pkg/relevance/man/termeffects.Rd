\name{termeffects}
\title{All Coefficients of a Model Fit}
\alias{termeffects}
\description{
A list of all coefficients of a model fit, possibly with respective statistics
}
\usage{
termeffects(object, se = 2, df = df.residual(object), rlv = TRUE)
}

\arguments{
  \item{object}{a model fit, produced, e.g., by a call to \code{lm} or
    \code{regr}.}
  \item{se}{logical: Should inference statistics be generated?}
  \item{df}{degrees of freedom for t-test}
  \item{rlv}{logical: Should relevances be calculated?}
}
\value{a \code{\link{list}} with a component for each term in the model formula.
  Each component is a \code{\link{termtable}} for the coefficients
  corresponding to the term.
}
\author{
Werner A. Stahel
}
\seealso{
  \code{\link{dummy.coef}, \link{termtable}}
}
\examples{
  data(d.blast)
  r.blast <-
    lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
  termeffects(r.blast)
}

\keyword{regression}
\concept{relevance}
