\name{linear.predictors}
\alias{linear.predictors}
\alias{linpred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{linear predictors from a (generalized) linear model
}
\description{extracts the linear.predictor component of a model object,
  taking 'na.resid' into account, in analogy to 'residuals' or
  'fitted.values'
}
\usage{
linear.predictors(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{model fit
}
}
\value{
  vector (or, for models inheriting from 'multinom', matrix)
  of linear predictor values
}
\author{Werner A. Stahel
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fitted.values}}
}
\examples{
## example from 'glm'
clotting <- data.frame(
u = c(5,10,15,20,30,40,60,80,100),
lot1 = c(118,58,NA,35,27,25,21,19,18), ## NA inserted instead of 42
lot2 = c(69,35,26,21,18,16,13,12,12))

r.gam  <- glm(lot1 ~ log(u), data = clotting, family = Gamma)
linear.predictors(r.gam)
## 8 elements; 3rd missing.
r.gex  <- glm(lot1 ~ log(u), data = clotting, family = Gamma,
              na.action=na.exclude)
linear.predictors(r.gex)
## 9 elements, third is NA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
