\name{regr0-internal}
\title{Internal regr0 objects}
\alias{i.glm}
\alias{i.lm}
\alias{i.multinomial}
\alias{i.survreg}
\alias{i.polr}
\alias{i.smooth}
%
\alias{is.formula}
\alias{\%nin\%}
\alias{summary.mreg}
\description{
  Internal functions of package \pkg{regr0} that are exported, but
  should normally not be called by the user directly.
}
\usage{
## i.glm(...)
## i.lm(...)
## i.multinomial(...)
## i.survreg(...)
## i.polr(...)
## i.smooth(...)

is.formula(object)
x \%nin\% y
summary.mreg(object, ...)
}
\arguments{
  \item{x,y, object}{\R objects}
  \item{...}{...}
}
\keyword{internal}
