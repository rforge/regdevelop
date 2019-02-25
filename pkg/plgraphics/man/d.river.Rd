\name{d.river}
\alias{d.river}
\docType{data}
\title{Chemical Compounds in a Swiss River, Time Series}
\description{This time series of chemical concentrations can be used
  to research the activities of photosynthesis and respiration
  in a river.
}
\usage{data("d.river")}
\format{
  A time series with 9792 observations (10 minutes interval)
  on the following 12 variables.
  \describe{
    \item{\code{date}}{Date of the observation}
    \item{\code{hour}}{Hour}
    \item{\code{pH}}{pH}
    \item{\code{O2}}{concentration of Oxygen}
    \item{\code{O2S}}{Oxygen saturation value}
    \item{\code{T}}{Temperature [deg C]}
    \item{\code{H2CO3}}{Carbon dioxide concentration in the water}
    \item{\code{CO2atm}}{Carbon dioxide concentration in the atmosphere}
    \item{\code{Q}}{flow}
    \item{\code{su}}{sunshine}
    \item{\code{pr}}{precipitation}
    \item{\code{ra}}{radiation}
  }
}
%%- \details{
%%-   The measurements stem from a station in the river Glatt near Zurich.
%%- }
\source{
  The measurements have been collected in the river Glatt near Zurich.
}
%%- \references{
%%- %%  ~~ possibly secondary sources and usages ~~
%%- }
\examples{
data(d.river)
range(d.river$date)
t.i <- d.river$date < as.Date("2010-03-31")

plyx(~date, ~O2, data=d.river[t.i & d.river$hour==14,], smooth=FALSE)

d.river$time <- gentimeaxis(d.river$date, hour=d.river$hour)
plyx(O2~time, data=d.river[t.i,], type="l")

plyx(O2+T+ra~date, data=d.river[t.i & d.river$hour==14,], 
  smooth.par=0.5, smooth.xtrim=0.03, ycol=c(O2="blue",ra="red"))
}
\keyword{datasets}
