\name{d.pollZH16}
\alias{d.pollZH16}
\alias{d.pollZH16d}
\docType{data}
\title{
pollution data, 1 year of hourly measurements
}
\description{
pollZH16.csv contains hourly air pollution measurements from a station
in the city center of Zurich, in a courtyard, for the year 2016.
pollZH16d.csv  is the subset of measurements for  hour=15 .
}
\usage{
data("d.pollZH16")
data("d.pollZH16d")
}
\format{
  A data frame with 8784 observations (366 for \code{d.pollZH16d})
  on the following 13 variables.
  \describe{
    \item{\code{date}}{Date}
    \item{\code{hour}}{hour of the day, 0-23 scale}
    \item{\code{O3}}{O3 concentration}
    \item{\code{NO2}}{NO2 concentration}
    \item{\code{temp}}{temperature, Celsius scale}
    \item{\code{rad}}{radiation}
    \item{\code{prec}}{precipitation}
    \item{\code{dateshort}}{character: Short identification of the day:
      A-L encodes the month;  1-9, a-x encodes the day of the month.
      For \code{d.pollZH16d}, these are used as \code{row.names}}
    \item{\code{weekday}}{... as factor, levels \code{Mon} to \code{Sun}}
    \item{\code{month}}{... as factor, levels \code{Jan} to \code{Dec}}
    \item{\code{sumhalf}}{logical: TRUE from Apr to Sep}
    \item{\code{sunday}}{logical: TRUE for Sundays}
    \item{\code{daytype}}{factor with levels \code{work}, for workdays,
    \code{Sat} and \code{Sun}}
}
}
\source{Bundesamt f\"ur Umwelt (BAFU), Schw. Eidgenossenschaft,
https://www.bafu.admin.ch/bafu/de/home/themen/luft/zustand/daten/datenabfrage-nabel.html
}
\examples{
data(d.pollZH16d)
plot(NO2~rad, data=d.pollZH16d)
}
\keyword{datasets}
