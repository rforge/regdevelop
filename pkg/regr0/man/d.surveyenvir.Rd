\name{d.surveyenvir}
\alias{d.surveyenvir}
\docType{data}
\title{Survey on Environment}
\description{
  A survey was conducted on attitudes towards environmental problems
}
\usage{data(d.surveyenvir)}
\format{
  A data frame with 2038 observations on the following 10 variables.
  \describe{
    \item{\code{age}}{numeric. age in years}
    \item{\code{sex}}{gender, either \code{m} for men or \code{f} for women}
    \item{\code{education}}{level of education, with levels
      \code{no.training}, \code{apprentiship}, \code{no.degree},
      \code{college}, \code{uni}}
    \item{\code{location}}{type of living environment, with levels
      \code{Zentrum}, \code{Stadt.dicht}, \code{Stadt.locker}
      \code{Rand.dicht}, \code{Rand.locker}, \code{Land}, \code{sonstig}}
    \item{\code{townsize}}{size of town, with categories
      \code{<2000}, \code{2000-4999}, \code{5000-19999}, \code{20000-49999}
      \code{50000-99999}, \code{100000-499999}, \code{>500000} inhabitants}
    \item{\code{party}}{political party, \code{CDU}, \code{SPD},
      \code{FDP}, \code{NPD}, \code{DKP}, \code{Gruene},
      \code{sonstige}, \code{keine}}
    \item{\code{disturbance}}{level of disturbance, \code{nicht}
      \code{etwas}, \code{ziemlich}, \code{sehr}} 
    \item{\code{gov}}{Is the government active enough?
      Levels \code{does.enough}, \code{not.enough}}
    \item{\code{responsability}}{Who has the main responsibility?
      Answers \code{individuals}, \code{government}, \code{both}} 
    \item{\code{weight}}{weight of the observation according to sampling
      design}
  }
}
%%- \details{
%%-   ~~ If necessary, more details than the __description__ above ~~
%%- }
\source{"Umweltschutz im Privatbereich.",
    Erhebung des EMNID, Zentralarchiv fuer empirische Sozialforschung
    der Universitaet Koeln 
}
%%- \references{
%%-   ~~ possibly secondary sources and usages ~~
%%- }
\examples{
data(d.surveyenvir)
t.r <- regr(disturbance~age+education+location, data=d.surveyenvir)
}
\keyword{datasets}
