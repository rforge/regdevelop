\name{predict.regr}
\alias{predict.regr}
%% \alias{predict.polr}
\alias{predict.mlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction Methods for 'regr' Objects}
\description{
  Calculates predicted values for \code{regr} objects.
  The possible results depend on the class of the fitted model. 
}
\usage{
\method{predict}{regr}(object, newdata = NULL, scale = object$sigma,
  df=object$df.residual, type = NULL, ...)

%% \method{predict}{polr}(object, newdata=NULL,
%%   type = c("class", "probs", "link"), ...)

\method{predict}{mlm}(object, newdata=NULL, se.fit = FALSE,
  scale = NULL, df = Inf,
    interval = c("none", "confidence", "prediction"), level = 0.95,
    type = c("response", "terms"), terms = NULL, na.action = na.pass,
    pred.var = NULL, weights = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{regr}}
  \item{newdata}{An optional data frame in which to look for variables with
          which to predict.  If omitted, the fitted values are used.}
  \item{se.fit}{if TRUE, standard errors will be calculated if possible}
  \item{scale}{Scale parameter for std.err. calculation}
  \item{df}{Degrees of Freedom to be used, %if \code{se.fit} is false,
    for quantiles defining intervals}
  \item{interval}{Type of interval requested}
  \item{level}{Confidence level}
  \item{type}{Type of prediction: response or model term\cr
  For glm or ordered regression, type \code{link} gives estimated values
  of the linear predictor.}
  \item{terms}{If \code{type="terms"}, which terms (default is all terms)}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}.  The default is to insert \code{NA}.}
  \item{pred.var, weights}{see \code{\link{lm}}}
  \item{\dots}{further arguments passed to specific methods}
}
\details{
  \code{regr} is a "super class" which includes many specific model
  classes such as "lm", "glm", "polr", ... .
  \code{predict.regr} is a wrapper function that calls the specific methods
  corresponding to the specific model class.
}
\value{
  vector of predictions, or matrix with columns \code{fit}, \code{lwr},
  and \code{upr} if \code{interval} is set.

  If \code{se.fit} is \code{TRUE}, a list with the
  following components is returned:  

  \item{fit}{vector or matrix as above}
  \item{se.fit}{standard error of predicted means}
  \item{residual.scale}{residual standard deviations}
  \item{df}{degrees of freedom for residual~Describe the value returned}
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{predict.lm}}, \code{\link{fitted.polr}}
}
\examples{
data(d.blast, package="plgraphics")
r.blast <- regr(log10(tremor)~location+log10(distance)+log10(charge),
            data=d.blast)
t.pr <- predict(r.blast)
plgraphics::showd(t.pr)

data(d.fossileSamples, package="plgraphics")
r.mregr <-
  regr(cbind(sAngle,lLength,rWidth)~SST+Salinity+lChlorophyll+Region+N,
                data=d.fossileSamples)
t.pr <- predict(r.mregr)
plgraphics::showd(t.pr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
