\name{leverage}
\alias{leverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get leverage values}
\description{
  Extracts the leverage component of a fit object using the
  \code{na.action} component if available
}
\usage{
leverage(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    an object containing a component \code{fit$leverage} and possibly
    a component \code{fit$na.action}
  }
}
\value{
  The vector \code{fit$leverage}, possibly expanded by missing values
  if \code{fit$na.action} has class \code{na.exclude}
}
\details{
  The difference to \code{hatvalues} is that \code{leverage} does not
  call \code{influence} and therefore does not require residuals.
  It is therefore simpler and more widely applicable.

  The function uses the \code{qr} decomposition of \code{object}.
  If necessary, it generate it.

  The leverage is the squared Mahalanobis distance of the observation
  from the center of the design X (\code{model.matrix}) with
  "covariance" X^T X. If there are weights (\code{object$weights}),
  the weighted center and "covariance" are used, and the distances are
  multiplied by the weights.
  To obtain the distances in the latter case, "de-weight" the leverages
  by dividing them by the weights.
}

\author{Werner A. Stahel, ETH Zurich}

\seealso{\code{\link{hat}; \link{hatvalues}; \link{influence}}}

\examples{
data(d.blast)
r.blast <-
  lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
t.lev <- leverage(r.blast)
length(t.lev)
dim(d.blast) ## number of observations is 388
length(r.blast$leverage) ## 362
  ## only observations without missings occur in $leverage
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
% \keyword{diagnostics}
