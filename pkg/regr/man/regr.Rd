\name{regr}
\alias{regr}
\alias{summary.regr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting Regression Models }
\description{
  'regr' fits regression models of various types: ordinary, robust,
  quantile, Tobit, logistic, 
  generalized linear, multinomial, ordered response, multivariate, ...

  It is a wrapper function which calls the respective S fitting
  functions and yields an extensive result, to be used by the
  respective print and plot methods.
}
\usage{
regr(formula, data=NULL, family=NULL, robust = FALSE, method=NULL, 
     nonlinear = FALSE, start=NULL, subset=NULL, weights=NULL, offset=NULL, ...)

\method{summary}{regr}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ data frame containing the variables in the model. }
  \item{family}{ character string describing the type of model and the
    fitting procedure. By default, \code{regr} will fit a model that is
    suitable for the response variable, if possible.\cr
    Can be either
    "gaussian" for ordinary regression, including multivariate response,
    "binomial" for logistic regression,
    "poisson"  for Poisson regression,
    "multinomial" for multinomial regression,
    "polr" for ordered response (cumulative logits) regression,
    "weibull", "lognormal", "loglogistic", "extreme"
    for survival regression with the respective distribution,
    "ph" for Cox regression,
    or a suitable argument for \code{glm}, i.e.: A
    description of the error distribution and link
    function to be used in the model. This can be a character string
    naming a family function, a family function or the result of a call
    to a family function.  (See 'family' for details of family
    functions.) }
  \item{nonlinear}{if TRUE, a non-linear model is expected, to be fitted
    by \code{nlm}
  }
  \item{start}{if \code{nonlin} is \code{TRUE}, starting values for
    \code{nlm}
  }
  \item{robust}{if TRUE, robust fitting is used if a function is
    available, see Details.
    }
    \item{method}{"rq" or "quantreg" requires fitting by quantile
  regression. Then, the argument \code{tau} is used for setting the
  quantile: \code{regr(..., method="quantreg", tau=0.9)}.
  Otherwise, if suitable, a method argument that is passed to the
  fitting procedure, see also Details.}
%%  \item{init.reg}{ function used to initialize robust methods }
  \item{subset}{ an optional vector specifying a subset of observations to be
    used in the fitting process. Names of variables in \code{data} can
    be used to generate the vector. }
  \item{weights}{ an optional vector of weights to be used for
    fitting. }
  \item{offset}{sometimes used to help specify the model.}
  \item{\dots}{ additional arguments, which are passed to %%the fitting function
    \code{regr.control}, see \code{?regr.control}}
  \item{object}{'summary.regr' simply returns its only argument
    \code{object}.}
}
% \details{}
\value{
  \code{regr} returns a list object of class \code{regr} and secondary class as
  produced by the fitting function. The components include generally
  those of the results of the fitting function and its summary.
  The important additional components are:
  \item{termtable}{data.frame,
    table for testing terms both for single and mulitple
    degree terms (continuous or binary explanatory variables and
    factors).
    The columns are
    The \code{print} method for the whole \code{regr} object will
    by default not show all of them.
}
  \item{allcoef}{All coefficients. This is a list with a component for
    each term. Each component characterizes the coefficients similar to
    the \code{termtable} just explained, but without the columns
    \code{df, stcoef, R2.x}.
    Note that the interpretation of the coefficients of dummy variables
    depends on the contrasts used for the factor.
    Tests need to be interpreted correspondingly, and with caution.}
  \item{stres}{Standardized residuals}
  \item{sigma}{estimated standard deviation for normal errors,\cr
    \code{sqrt(dispersion)} for \code{glm}}
  \item{h}{leverage values}
  \item{allvars}{the variables used for fitting the model}
  \item{binfac}{levels of binary factors}
  \item{fitfun, funcall}{R function and call that has been used for
    fitting the model}
}
\details{
  Contrasts: If (ordered) factors are generated in the \code{formula}
  (like \code{y~ordered(f1)}), the contrast functions
  \code{\link{contr.wsum}} and \code{\link{contr.wpoly}} will be
  equivalent to their unweighted versions.
  
  Robust fitting is available for ordinary linear models, by setting
  \code{robust = TRUE}.
  By default, the function \code{\link[robustbase]{lmrob}} from
  package \pkg{robustbase}
  is called with the argument \code{setting="KS2011"}.
  This default will be changed in future versions.
  Alternative methods can be obtained by setting the \code{method}
  argument. Choosing
  \code{method="rlm"}
  {calls function \code{\link[MASS]{rlm}} from package \pkg{MASS}}
  with default \code{method="MM"}.
  Any further elements of \code{method} will be passed on to the methods
  argument of the called function.

  The 'summary' generic function does nothing to 'regr' objects, since
  all useful information is already collected in the object.

  Logical explanatory variables are converted into numeric, in contrast
  to the official R fitting functions, which treat them as factors.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Werner A. Stahel, Seminar for Statistics, ETH Zurich}
%% \note{ ~~further notes~~
%%  ~Make other sections like Warning with \section{Warning }{....} ~
%% }
\seealso{
  \code{\link{lm}}, \code{\link{glm}}, \code{\link[MASS]{rlm}},
  \code{\link[nnet]{multinom}}, \code{\link[MASS]{polr}},... }
\examples{
data(d.blast)
( r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast) )

## Anorexia
data(anorexia, package="MASS")
r.anorexia <- regr(Postwt ~ Prewt + Treat + offset(Prewt),
                   data = anorexia)

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
     ## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
d.dob <- data.frame(group = gl(2,10,20, labels=c("Ctl","Trt")),
                    weight = c(ctl, trt))
(r.dob <- regr(weight ~ group, data=d.dob))

## multinomial regression
data(d.surveyenvir)
d.surveyenvir$dist.unordered <- factor(as.character(d.surveyenvir$disturbance))
t.r <- regr(dist.unordered~age+education+location, data=d.surveyenvir)

## ordered regression
data(housing, package="MASS")
t.r <- regr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
plot(t.r)

## multivariate regression
data(d.fossileSamples, package="plgraphics")
r.mregr <-
  regr(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+lChlorophyll+region+N,
                data=d.fossileSamples)
}

\keyword{ regression }

