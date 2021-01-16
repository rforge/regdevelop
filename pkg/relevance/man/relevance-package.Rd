\name{relevance-package}
\alias{relevance-package}
\alias{relevance}
\docType{package}
\title{\packageTitle{relevance}}
\description{
  \packageDescription{relevance}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{relevance}
\packageIndices{relevance}

Relevance is a measure that expresses the (scientific) relevance of an effect.
The simplest case is a single sample of supposedly normally distributed
observations, where interest lies in the expectation,
estimated by the mean of the observations.
There is a threshold for the expectation, below which an effect is
judged too small to be of interest.

The estimated relevance \sQuote{\eqn{Rle}} is then simply the estimated effect divided by
the threshold. If it is larger than 1, the effect is thus judged
relevant. The two other values that characterize the relevance are the
limits of the confidence interval for the true value of the relevance,
called the secured relevance  \sQuote{\eqn{Rls}} and the potential relevance \sQuote{\eqn{Rlp}}.

If \eqn{Rle > 1}, then one might say that the effect is
\dQuote{significantly relevant}.

Another useful measure, meant to replace the p-value, is the
\dQuote{significance} \sQuote{Sg0}. In the simple case, it divides the
estimated effect by the critical value of the (t-) test statistic.
Thus, the statistical test of the null hypothesis of zero expectation
is significant if \sQuote{Sg0} is larger than one, \eqn{Sg0 > 1}.

These measures are also calculated for the comparison of two groups,
for proportions, and most importantly for regression models.
For models with linear predictors, relevances are obtained for
standardized coefficients as well as for the effect of dropping terms
and the effect on prediction.

The most important functions are
\describe{
  \item{\code{twosamples()}: }{
    calculate the measures for two paired or unpaired sampless or
    a simple mean. This function calls}
    \item{ \code{ciSgRl()}: }{
      calculates the confidence interval and siginificance based on an
      estimate and a standard error,
      and adds relevance for a standardized effect.}
    \item{\code{termtable()}: }{
      deals with fits of regression models with a linear predictor.
      It calculates confidence intervals and significances for
      the coefficients of terms with a single degree of freedom.
      It includes the effect of dropping each term
      (based on the \code{drop1} function)
      and the respective significance and relevance measures.}
    \item{\code{termeffects()}: }{
      calculates the relevances for the coefficients
      related to each term. These differ from the enties of \code{termtable}
      only for terms with more than one degree of freedom.}
  }
}
\author{
\packageAuthor{relevance}

Maintainer: \packageMaintainer{relevance}
}
\references{
  Stahel, Werner A. (2021).
  New relevance and significance measures to replace p-values.
  To appear in PLoS ONE
}
\seealso{Package \pkg{regr}, avaiable from \url{https://regdevelop.r-forge.r-project.org}}
\examples{
  data(swiss)
  rr <- lm(Fertility ~ . , data = swiss)
  termtable(rr)
}
\keyword{ package }
\concept{relevance}
