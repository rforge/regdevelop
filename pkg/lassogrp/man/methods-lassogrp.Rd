\name{methods-lassogrp}
\title{Prediction and Other Methods for 'lassogrp' Objects}
\alias{predict.lassogrp}
\alias{fitted.lassogrp}
\alias{formula.lassogrp}
\description{
  Compute fitted values or predictions for \code{"lassogrp"} objects
  (typically resulting from \code{\link{lasso}()}), or get
  (\dQuote{extract}) other information from such an object.
}
\usage{
\method{fitted}{lassogrp}(object, \dots)

\method{predict}{lassogrp}(object, newdata = NULL,
  type = c("link", "response"), na.action = na.pass, \dots)

\method{formula}{lassogrp}(x, \dots)
}
\arguments{
  \item{object,x}{an object of class \code{"lassogrp"}, containing a lasso fit}
  \item{newdata}{data.frame containing the new regressors for which the
    predictions are required.}
  \item{type}{type of prediction}
  \item{na.action}{action for missing values, see \code{\link{lm}}.}
  \item{\dots}{not used (needed for compatibility)}
}
\details{The fitted values are equal to response predictions for the
  observations used for fitting.
}
\value{matrix of fitted or predicted values. Columns correspond to lambda values,
  rows, to observations in \code{newdata}.}
\author{Werner Stahel, \email{stahel@stat.math.ethz.ch}}
\examples{
data(asphalt)
dd <- asphalt
rr <- lasso(log10(RUT)~log10(VISC)+ASPH+BASE+FINES+VOIDS+RUN, data=dd)
fitted(rr)
predict(rr, newdata=
  data.frame(VISC=2, ASPH=5, BASE=5, VOIDS=5, FINES=70, RUN=1))
}
\keyword{misc}
