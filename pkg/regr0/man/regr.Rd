\name{regr}
\alias{regr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting Regression Models }
\description{
  'regr' fits regression models of various types: ordinary, robust,
  generalized linear, multinomial, ordered response, multivariate, ...

  It is a wrapper function which calls the respective S fitting
  functions and yields an extensive result, to be used by the
  respective print and plot methods.
}
\usage{
regr(formula, data, tit = NULL, family = NULL,
  calcdisp = NULL, suffmean = 3, nonlinear = FALSE, start=NULL,
  robust = FALSE, method=NULL, init.reg = "f.ltsreg",
  subset = NULL, weights = NULL, na.action = nainf.exclude,
  model = TRUE, x = FALSE, termtable = TRUE, vif = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ data frame containing the variables in the model. }
  \item{tit}{ title (becomes tit attribute of result) }
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
  \item{calcdisp}{ a logical value indicating whether, for
    \code{family=binomial} and \code{family=poisson}, dispersion
    should be calculated \code{(TRUE)} or set to 1 \code{(FALSE)} }
  \item{suffmean}{!!!}
  \item{nonlinear}{if TRUE, a non-linear model is expected, to be fitted
    by \code{nlm}
  }
  \item{start}{if \code{nonlin} is \code{TRUE}, starting values for
    \code{nlm}
  }
  \item{robust}{if TRUE, robust fitting is used if a function is
    available, see Details.
    }
  \item{method}{if suitable, a method argument that is passed to the
  fitting procedure, see Details.}
  \item{init.reg}{ function used to initialize robust methods }
  \item{subset}{ an optional vector specifying a subset of observations to be
    used in the fitting process. Names of variables in \code{data} can
    be used to generate the vector. }
  \item{weights}{ an optional vector of weights to be used for fitting. }
  \item{na.action}{ a function which indicates what should happen when the data
          contain 'NA's.  The default is set by the 'na.action' setting
          of 'options', and is 'na.fail' if that is unset.  The
          "factory-fresh" default is 'na.omit'.  Another possible value
          is 'NULL', no action.  Value 'na.exclude' can be useful. }
  \item{model, x}{ logical values indicating whether the model frame or
  the model.matrix, respectively, should be
          included as a component of the returned value. }
  \item{termtable}{if \code{FALSE}, the \code{regr} type term table will
    not be generated}
  \item{vif}{ a logical value indicating whether the collinearity
    measure R2.j ahould be calculated }
  \item{\dots}{ other argument to be passed to the fitting function }
}
% \details{}
\value{
  \code{regr} returns a list object of class \code{regr} and secondary class as
  produced by the fitting function. The components include generally
  those of the results of the fitting function and its summary.
  Additional components are described here in a draft version.
  \item{stres }{Standardized residuals}
  \item{sigma }{estimated standard deviation for normal errors,\cr
  sqrt(dispersion) for \code{glm}}
  \item{testcoef}{Table for testing terms both for single and mulitple
    degree terms (continuous or binary explanatory variables and
    factors)}
  \item{allcoef}{All coefficients. Estimated effects of factors are
    given for all their levels.}
  \item{h}{leverage values}
  \item{vartype}{types of the input variables (c for continuous,
    b for binary, f for factor, o for ordered factor)}
  \item{binfac}{levels of binary factors}
  \item{fitfun}{R function that has been called for fitting the model}
}
\details{
  Robust fitting is available for ordinary linear models, by setting
  \code{robust = TRUE}.
  By default, the function \code{\link[robustbase]{lmrob}} from package \pkg{robustbase}
  is called with the argument \code{setting="KS2011"}.
  Alternative methods can be obtained by setting the \code{method}
  argument. Choosing
  \code{method="rlm"}
  {calls function \code{\link[MASS]{rlm}} from package \pkg{MASS}}
  with default \code{method="MM"}.
  Any further elements of \code{method} will be passed on to the methods
  argument of the called function.
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
t.r <- regr(disturbance~age+education+location, data=d.surveyenvir)

## ordered regression
data(housing, package="MASS")
t.r <- regr(Sat ~ Infl + Type + Cont, weights = housing$Freq,
            data = housing)
plot(t.r)

## multivariate regression
data(d.fossiles)
r.mregr <-
  regr(cbind(sAngle,lLength,rWidth)~SST.Mean+Salinity+lChlorophyll+region+N,
                data=d.fossiles)
}

\keyword{ regression }

