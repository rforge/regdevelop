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
\alias{print.mregr}
\alias{DB}
\alias{BR}
\alias{IR}
\alias{subset}
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
\method{print}{mregr}(x, na.print=getUserOption("na.print"), ...)
DB(on=TRUE)
BR()
IR(condition)
subset(x, ...)
}
\arguments{
  \item{x,y, object}{\R objects}
  \item{na.print}{see print.regr}
  \item{condition}{condition}
  \item{...}{...}
}
\keyword{internal}
