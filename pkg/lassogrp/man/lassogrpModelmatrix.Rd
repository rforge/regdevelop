\name{lassogrpModelmatrix}
\Rdversion{1.1}
\alias{lassogrpModelmatrix}
\title{
Model Matrix for Lasso
}
\description{
Generates the model matrix and adds information about variables to be
penalized or non-penalized by the L1 term in the lasso fitting.
The user will rarely call this function.
}
\usage{
lassogrpModelmatrix(m, formula, nonpen = ~1, data, weights, subset,
  na.action, contrasts, env)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{call to the lasso function  
}
  \item{formula}{model formula
}
  \item{nonpen}{formula defining the non-penalized variables
}
  \item{data}{data.frame in which the variables are found
}
  \item{weights}{as in lm
}
  \item{subset}{as in lm
}
  \item{na.action}{as in lm
}
  \item{contrasts}{as in lm
}
  \item{env}{environment in which the formula is evaluated
}
}
\details{
  This function generates the model.matrix and the 'index' vector which
  stores the information about penalization
}
\value{
  \item{x}{model.matrix}
  \item{y}{response}
  \item{w}{weights}
  \item{off}{offset}
  \item{mf}{model.frame}
  \item{index}{vector of length \code{ncol(x)} determining groups for
    penalization. Negative values correspond to non-penalized variables.}
}
\author{Lukas Meier, modified by Werner Stahel\email{stahel@stat.math.ethz.ch}
}
\examples{
## not to be called by the user. Is called from function lasso.
}
\keyword{regression}
