\name{fitted.regr}
\alias{fitted.regr}
\alias{fitted.polr}
\title{Fitted Values for regr and polr objects}
\description{
  \code{fitted.regr} extracts fitted values from \code{regr} objects.
  The function is needed for results of survival models, including Tobit
  models. 
  \code{fitted.polr} calculates Fitted Values for ordered response regression.
  In addition to the types of the usual version, this function allows
  for extracting the \code{link} type of predictions.
}
\usage{
\method{fitted}{regr}(object, ...)
\method{fitted}{polr}(object, type = "link", na.action = object, ...)
}
\arguments{
  \item{object}{an object of class \code{regr} (for \code{fitted.regr}),
  or one that inherits from class \code{polr} (for \code{fitted.polr})}
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
  \item{\dots}{further arguments, passed to \code{fitted} methods called
  by \code{fitted.regr}; not used in \code{fitted.polr}}
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
