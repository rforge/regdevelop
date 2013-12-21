\name{fitted.polr}
\alias{fitted.polr}
\title{Fitted Values for polr objects}
\description{
  Calculates Fitted Values for ordered response regression.
  In addition to the types of the usual version, this function allows
  for extracting the \code{link} type of predictions.
}
\usage{
\method{fitted}{polr}(object, type = "link", na.action = object, ...)
}
\arguments{
  \item{object}{an object that inherits from class \code{polr}}
  \item{type}{type of fitted values (as with \code{predict}).
    \describe{
      \item{\code{link}}{asks for the linear predictor,}
      \item{\code{class}}{returns the most likely (integer) Y values}
      \item{\code{probs}}{returns probabilities for each value of Y}
    }

  }
  \item{na.action}{object that possibly contains an \code{na.action}
    component. If it is of the class \code{exclude}, the fitted values
    will be expanded to fit the number of rows in the data.frame used
    for fitting the regression.}
  \item{\dots}{not used}
}
\value{
  Vector of fitted values
}
\author{Werner A. Stahel, Seminar for Statistics, ETH Zurich}
\examples{
require(MASS) # -> dataset 'housing' and polr() :
house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
house.fit <- fitted(house.plr, type="link") # only possible with regr version
summary( house.fit )
}
\keyword{regression}
