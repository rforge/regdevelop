\name{step}
%% \alias{step}
\alias{step.regr}
%% \alias{step.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Choose a model in a Stepwise Algorithm
}
\description{
  Select a formula-based model
}
\usage{
%%  step(object, ...)
%% \method{step}{default}(object, scope, scale = 0,
%%   direction = c("both", "backward", "forward"), trace = 1, keep = NULL,
%%   steps = 1000, k = 2, ...) 
%% \method{step}{regr}
step.regr(object, scope = NULL, expand=FALSE, scale = 0,
  direction = c("both", "backward", "forward"), trace = FALSE, keep = NULL,
  steps = 1000, k = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object representing a model of an appropriate class.
    This is used as the initial model in the stepwise search.
  }
  \item{scope}{defines the range of models examined in the stepwise search.
    This should be either a single formula, or a list containing
    components 'upper' and 'lower', both formulae.  See the
    details for how to specify the formulae and how they are
    used.
  }
  \item{expand}{logical. If TRUE and \code{scope} is not specified,
    all squared terms and interactions are used as upper scope
    (using \code{\link{terms2order}})
    and the formula of the \code{object} as the lower.}    
  \item{scale}{used in the definition of the AIC statistic for selecting the
    models, currently only for 'lm', 'aov' and 'glm' models.  The
    default value, '0', indicates the scale should be estimated:
    see 'extractAIC'.
  }
  \item{direction}{the mode of stepwise search, can be one of '"both"',
    '"backward"', or '"forward"', with a default of '"both"'.  If
    the 'scope' argument is missing the default for 'direction'
    is '"backward"'.
  }
  \item{trace}{if positive, information is printed during the running of
    'step'.  Larger values may give more detailed information.
  }
  \item{keep}{a filter function whose input is a fitted model object and
    the associated 'AIC' statistic, and whose output is
    arbitrary.  Typically 'keep' will select a subset of the
    components of the object and return them. The default is not
    to keep anything.
  }
  \item{steps}{the maximum number of steps to be considered.  The default is
    1000 (essentially as many as required).  It is typically used
    to stop the process early.
  }
  \item{k}{the multiple of the number of degrees of freedom used for the
    penalty.  Only \code{k = 2} gives the genuine AIC, \code{k = log(n)}
    is sometimes referred to as BIC or SBC.
    \code{k=4} usually leads to models with formally significant terms.
  }
  \item{\dots}{any additional arguments to 'extractAIC'.
  }
}
\details{
  The set of models searched is determined by the 'scope' argument.
  The right-hand-side of its 'lower' component must be included in
  the model, and right-hand-side of the model is included in the
  'upper' component.  If 'scope' is a single formula, it specifies
  the 'upper' component, and the 'lower' model is empty.  If 'scope'
  is missing, the initial model is used as the 'upper' model.

  Models specified by 'scope' can be templates to update 'object' as
  used by 'update.formula'.  So using '.' in a 'scope' formula means
  'what is already there', with '.^2' indicating all interactions of
  existing terms.

  Missing values lead to a reduced dataset:
  \code{step.regr} works on the dataset that includes all variables
  appearing in \code{scope} and then drops all lines containing
  missing values (by applying \code{na.omit}).
  The result is the model re-fitted to the dataset with only
  the variables used in the final model.
  This may lead to an increased number of rows.
  This differs from the behavior of \code{step} of package
  \code{stats}.

  [from \code{step} of package \code{stats}:]
  There is a potential problem in using 'glm' fits with a variable
  'scale', as in that case the deviance is not simply related to the
  maximized log-likelihood.  The '"glm"' method for function
  'extractAIC' makes the appropriate adjustment for a 'gaussian'
  family, but may need to be amended for other cases.  (The
  'binomial' and 'poisson' families have fixed 'scale' by default
  and do not correspond to a particular maximum-likelihood problem
  for variable 'scale'.)
}
\value{
     The selected model is returned, with up to two additional
     components.  There is an
     \item{anova}{steps taken in the search,}
     \item{keep}{if the \code{keep=} argument was supplied in the call.}\cr
       The '"Resid. Dev"'
     column of the analysis of deviance table refers to a constant
     minus twice the maximized log likelihood: it will be a deviance
     only in cases where a saturated model is well-defined (thus
     excluding 'lm', 'aov' and 'survreg' fits, for example).
}
\author{Werner A, Stahel}

\seealso{
  \code{step} in package \code{stats}, \code{stepAIC} in package
  \code{MASS}, \code{\link{add1.regr}}, \code{\link{drop1.regr}}
}
\examples{
r.fit <- regr(Fertility ~ ., data = swiss)
r.step <- step.regr(r.fit)
r.step$anova
r.st2 <- step.regr(r.fit, k=8, trace=FALSE)
r.st2
% the following causes an error in R CMD check, but runs when I try it.
% r.st3 <- step(r.fit, scope=terms2order(r.fit), k=5, trace=FALSE)
% r.st3 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}

