\name{regr-internal}
\title{Internal regr objects}
\alias{i.def}
\alias{i.glm}
\alias{i.lm}
\alias{i.multinomial}
\alias{i.multinomfit}
\alias{i.survreg}
\alias{i.polr}
\alias{i.polrfit}
\alias{i.smooth}
\alias{notice}
%
\alias{print.mregr}
%% \alias{subset}
\description{
  Internal functions of package \pkg{regr} that are exported, but
  should normally not be called by the user directly.
}
\usage{
## i.glm(...)
## i.lm(...)
## i.multinomial(...)
## i.survreg(...)
## i.polr(...)
## i.smooth(...)

\method{print}{mregr}(x, na.print=getOption("na.print"), ...)
notice(..., notices = NULL)
}
\arguments{
  \item{x,y, object}{\R objects}
  \item{na.print}{see print.regr}
  \item{condition}{condition}
  \item{...}{...}
}
\keyword{internal}
