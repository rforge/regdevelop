\name{relevance-package}
\alias{relevance-package}
\alias{relevance}
\docType{package}
\title{
\packageTitle{relevance}
}
\description{
\packageDescription{relevance}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{relevance}
\packageIndices{relevance}
This package is being built.

Relevance is a measure that expresses the relevance of an effect.
The simplest case is a single sample supposedly normally distributed
observations, and we are interested in the expectation,
estimated by the mean of the observations.
There is a threshold for the expectation, below which an effect is
judged too small to be of interest.
The estimated relevance Rle is then simply the estimated effect divided by
the threshold. If it is larger than 1, the effect is thus judged
relevant. The two other values that characterize the relevance are the
limits of the confidence interval for the relevance, called
the secured relevance Rls and the potential relevance Rlp.
If Rle $>1$, then one might say that the effect is
"significantly relevant".

Another useful measure, meant to replace the p-value, is the
"significance" Sg0. In the simple case, it divides the
estimated effect by the critical value of the (t-) test statistic.
Thus, the statistical test of the null hypothesis of zero expectation
is significant if Sg0 is $>1$.

These measures are also calculated for the comparison of two groups,
for proportions, and most importantly for regression models.
For models with linear predictors, relevances are obtained for
standardized coefficients as well as for the effect of dropping terms.
%% The results are presented in tables

The most important functions are

rlstats

For models with a linear predictor, the results are given by tables.
A coeftable is produced by \code{rlstats}.
\code{termtable} contains the relevance of the terms as characterized by
the effect of dropping the term from the model formula.
It also includes the relevance of the coefficients as produced by
\code{rlstats} for terms with a single degree of freedom.
\code{termeffects} calculates the relevances for the coefficients
related to each term and is thus only of interest for terms with
more than one degree of freedom.
}
\author{
\packageAuthor{relevance}

Maintainer: \packageMaintainer{relevance}
}
\references{
  Stahel, Werner A. (2020). 
  Measuring Significance and Relevance instead of p-values.
  In preparation.
}
\keyword{ package }
%% \seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%% }
\examples{
  data(swiss)
  rr <- lm(Fertility ~ . , data = swiss)
  rt <- termtable(rr)
  rt
  names(rt)
  data.frame(rt)
}
