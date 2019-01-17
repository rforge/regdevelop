\name{pllimits}
\alias{pllimits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine Inner Plot Range 
}
\description{The inner plotting range is the range in which
  plotting functions of the regr0 package show unmodified coordinates.
  This function determines the range for one or more variables.
}
\usage{
pllimits(pllim, data, limfac = 4, FUNC=robrange)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pllim}{either a logical: shall an inner plotting range be
    determined? -- or a matrix with 2 rows and \code{NCOL(data)} rows,
    in which case the suitability will be checked.}
  \item{data}{vector or matrix or data.frame of data for which the inner
    plotting range is to be determined}
  \item{limfac}{scalar factor by which the range determined by
    \code{FUNC} is expanded}
  \item{FUNC}{function that determines the range of the data}
}
%%- \details{
%%- %%  ~~ If necessary, more details than the description above ~~
%%- }
\value{A matrix with 2 rows containing the minimum and the maximum
  of the inner plotting range. The columns correspond to those in
  \code{data}.
}
\author{WSt}
%%- \note{
%%- %%  ~~further notes~~
%%- }
\seealso{
\code{\link{plcoord}}
}
\examples{
  set.seed(0)
  xx <- rt(50, df=3)
  ( pll <- pllimits(TRUE, xx) )
  sum(xx<pll[1,] | xx>pll[2,])  ## 3
}
\keyword{ utitlities }

