\name{regr0-internal}
\title{Internal regr0 objects}
\alias{i.glm}
\alias{i.lm}
\alias{i.multinomial}
\alias{i.multinomfit}
\alias{i.survreg}
\alias{i.polr}
\alias{i.polrfit}
\alias{i.smooth}
%
\alias{is.formula}
\alias{\%nin\%}
\alias{summary.mregr}
\alias{print.mregr}
\alias{DB}
\alias{BR}
\alias{IR}
\alias{nna}
\alias{subset}
\alias{extract.AIC.regr}
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
## summary.mregr(object, ...)
}
\arguments{
  \item{x,y, object}{\R objects}
  \item{...}{...}
}
\keyword{internal}
