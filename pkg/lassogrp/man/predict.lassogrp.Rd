\name{predict.lassogrp}
\Rdversion{1.1}
\alias{predict.lassogrp}
\alias{fitted.lassogrp}

\title{Fitted values and Prediction for 'lassogrp' Objects}
\description{Calculate fitted values or predictions for 'lassogrp' objects}
\usage{
\method{fitted}{lassogrp}(object, ...)

\method{predict}{lassogrp}(object, newdata = NULL,
  type = c("link", "response"), na.action = na.pass, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class 'lassogrp' containing a lasso fit}
  \item{newdata}{data.frame containing the new regressors for which the
    predictions are required}
  \item{type}{type of prediction}
  \item{na.action}{action for missing values, see \code{?lm}}
  \item{\dots}{not used (needed for compatibility)}
}
\details{The fitted values are equal to response predictions for the
  observations used for fitting}
\value{matrix of fitted or predicted values. Columns correspond to lambda values,
  rows, to observations in 'newdata'}
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
