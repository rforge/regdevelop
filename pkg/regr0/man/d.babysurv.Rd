\name{d.babysurv}
\alias{d.babysurv}
\alias{d.babysurvGr}
\docType{data}
\title{Survival of early born babies
}
\description{
  The survival of prematurely born babies is modeled based on several
  characteristics. 
}
\usage{data("d.babysurv")
data("d.babysurvGr")}
\format{
  \code{d.babysurv} is a data frame with 246 observations on the
  following 6 variables. To be checked!!!
  \describe{
    \item{\code{Survival}}{Binary variable to be modelled: The babies
      with \code{Survival==1} have survived}
    \item{\code{Weight}}{Weight}
    \item{\code{Age}}{Age in ...}
    \item{\code{Apgar1, Apgar2}}{Two scores that desribe the vitality
      of the baby}
    \item{\code{pH}}{Blood pH at birth}
  }
  \code{d.babysurvGr} contains summarized data: All babies in the same
    weight class of width 100 are grouped. The other variables are
    neglected.
  \describe{
    \item{n}{Number of babies in each weight class.}
    \item{Survival.0, Survival.1}{Number of babies that have not survived
      in each weight class.}
    \item{Weight}{Mean of weight class}
  }
}
\details{
  ...
  }
\source{
Hibbard (1986)
}
\references{
Hibbard (1986)
}
\examples{
data(d.babysurv)
data(d.babysurvGr)
plot(I(Survival.1/n)~Weight, data=d.babysurvGr)
}
\keyword{datasets}
