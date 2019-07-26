\name{showd}
\alias{showd}
\title{Show a Part of a Data.frame}
\description{
  Shows a part of the data.frame which allows for grasping the nature of
  the data. The function is typically used to make sure that the data is
  what was desired and to grasp the nature of the variables in the phase
  of getting acquainted with the data.
}
\usage{
showd(data, first = 3, nrow. = 4, ncol. = NULL, digits=getOption("digits"))
}
\arguments{
  \item{data}{a data.frame, a matrix, or a vector}
  \item{first}{the first \code{first} rows will be shown and ...}
  \item{nrow.}{a selection of \code{nrow.} rows will be shown in
    addition. They will be selected with equal row number differences.
    The last row is always included.}
  \item{ncol.}{number of columns (variables) to be shown. The first and
    last columns will also be included. If \code{ncol.} has more than
    one element, it is used to identify the columns directly.}
  \item{digits}{number of significant digits used in formatting numbers}
}
\details{
  The \code{tit} attribute of \code{data} will be printed if available and
  \code{getUserOption("doc") > 0}, and any \code{doc} attribute,
  if \code{getUserOption("doc") >= 2} (see \code{\link{tit}}).
}

\value{
  returns invisibly the character vector containing the formatted data
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{
  \code{\link{head}} and \code{\link{tail}}.
}
\examples{
showd(iris)

data(d.birthrates)
names(d.birthrates)
## only show 7 columns, including the first and last
showd(d.birthrates, ncol=7)  

showd(cbind(1:100))
}
% R documentation directory.
\keyword{utilities}
\keyword{print}
