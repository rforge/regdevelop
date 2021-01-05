\name{weekday}
\alias{weekday}
\alias{ymd}
\title{Get Day of Week or Year, Month, Day
}
\description{From Dates, obtain the day of the week or
  the year, month and day
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
weekday(date, month = NULL, day = NULL, out = NULL, factor = FALSE)
ymd(date)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{date}{date(s), given as a \code{Date} object, a character object
    that can be converted into a \code{Date}, or as julian,
    or the year, in which case \code{month} and \code{day} must be given.
  }
  \item{month, day}{If the first argument is the year, these arguments
    must also be given.
  }
  \item{factor}{logical: Should the result be a (ordered) factor?
  }
  \item{out}{selection of output: either
    \code{"numeric"} for numeric output
    (0 for Sunday, ... 6 for Saturday),
    \code{"full"} or \code{"long"} for full weekday names, or
    \code{"short"} for abbreviated (3 character) names.
    Defaults to \code{"full"} if \code{factor} is \code{TRUE},
    to \code{"numeric"} otherwise.
  }
}
\value{For \code{weekdays},
  the output is as described above, depending on
  \code{factor} and \code{out}.
}
\author{Werner A. Stahel
}
\note{The functions call functions from the \code{chron}
  package
}
\seealso{
  \code{\link[chron]{day.of.week}, \link[chron]{month.day.year}}
}
\examples{
weekday(c("2020-05-01", "2020-05-02"), factor=TRUE)
## [1] Thursday Sunday  
## Levels: Sunday < Monday < Tuesday < Wednesday < Thursday < Friday < Saturday

dt <- ymd(18100+1:5)
weekday(dt)
## [1] 3 4 5 6 0
}
\keyword{manip}
\keyword{chron}
