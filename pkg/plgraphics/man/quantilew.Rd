\name{quantilew}
\alias{quantilew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quantiles for weighted observations
%%  ~~function to do ... ~~
}
\description{
  Quantiles for weighted observations
  }
\usage{
quantilew(x, probs = c(0.25, 0.5, 0.75), weights = 1, na.rm=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric vector whose sample quantiles are wanted
    'NA' and 'NaN' values are not allowed unless 'na.rm' is 'TRUE'.

%%     ~~Describe \code{data} here~~
}
\item{probs}{
  numeric vector of probabilities with values in [0,1].
}
  \item{weights}{
  numeric vector of weights. They will be standardized to sum to 1.
}
  \item{na.rm}{
  remove NAs from 'x'? If FALSE and 'x' contains NAs, the value will be NA.
}
}
\value{
  Empirical quantiles corresponding to the given probabilities and
  weights. If a quantile is not unique since the cumulated weights
  hit the probability value exactly (the case of the median of a sample
  of even size), the mean of the corresponding values is returned.  
}
\author{
Werner A. Stahel
}

\seealso{
  \code{\link{quantile}}
}
\examples{
   x <- c(1,3,4,8,12,13,18,20)
   quantile(x, c(0.25, 0.5))
   quantilew(x, c(0.25, 0.5), weights=1:8)  ##  8  13
   ## relative weights  (1+2+3)/36  sum to  <0.25 , with the forth, they
   ##   are over 0.25, therefore, the quantile is the 4th value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}

