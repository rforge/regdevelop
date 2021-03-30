\name{relevance-internal}
\alias{DB}
%% ---
\title{Internal objects of package relevance}
\description{
%%  Internal functions of package \pkg{relevance} that are exported, but
%%  should normally not be called by the user directly.
\code{DB} helps debugging functions by changing the \code{error} option.
}
\usage{
DB(on=TRUE)
}
\arguments{
  \item{on}{
    \describe{
      \item{0 or \code{FALSE}}{for switching debugging off,}
      \item{1 or \code{TRUE}}{for switching debugging on,}
      \item{2}{for turning warnings into errors}
    }
  }
}
\value{No return value, called for side effects}
\author{\packageAuthor{relevance}}
\keyword{internal}
