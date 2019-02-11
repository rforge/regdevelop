\name{d.babysurvival}
\alias{d.babysurvival}
\alias{d.babysurvGr}
\docType{data}
\title{Survival of Premature Infants}
\description{Survival of Premature Infants to be modeled
  using 5 potential explanatory variables. 
}
\usage{
  data("d.babysurvival")
  data("d.babysurvGr")
}
\format{
  \code{d.babysurvival}: 
  A data frame with 246 observations on the following 6 variables.
  \describe{
    \item{\code{Survival}}{binary, 1 means the infant survived}
    \item{\code{Weight}}{birth weight [g]}
    \item{\code{Age}}{pregnancy in weeks}
    \item{\code{Apgar1}}{A score indication the fitness of the infant at
    birth, scores 0 to 9}
    \item{\code{Apgar5}}{alternative score}
    \item{\code{pH}}{blood pH}
  }
  \code{d.babysurvGr}: 
  Grouped data: Number of Infants that died and survived for each class
  of birth weight.
  \describe{
    \item{\code{n}}{Number of infants in the weight class}
    \item{\code{Survival.0, Survivl.1}}{Number of infants that died and
      survived, respectively}
    \item{\code{Weight}}{birth weight}
  }
}
%%- \details{
%%- %%  ~~ If necessary, more details than the __description__ above ~~
%%- }
\source{Hibbard (1986)
}
%% \references{
%% %%  ~~ possibly secondary sources and usages ~~
%% }
\examples{
data(d.babysurvival)
summary(d.babysurvival)
rr <- glm(Survival~Weight+Age+Apgar1, data=d.babysurvival, family="binomial")
plot.regr(rr, xvar= ~Age+Apgar1)
}
\keyword{datasets}
