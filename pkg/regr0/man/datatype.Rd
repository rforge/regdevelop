\name{datatype}
\alias{datatype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine the Types of Variables}
\description{
  Determines the types of the variables (columns) of a data.frame
}
\usage{
datatype(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data.frame}
}
\details{
  A type is determined for every variable (column) of the data.frame.
  The types are:
  \describe{
    \item{c}{continuous variable}
    \item{b}{binary variable. There is no distinction between factors
      and logical variables or numerical variables having only 2 distint
      values.} 
    \item{f}{factor (not ordered)}
    \item{o}{ordered factor}
  }
  This information is used by \code{\link{regr}} to determine the
  adequate handling of the variables for modelling.
}
\value{
  the data.frame with two additional attributes: containing
  \item{vartype}{the variable types as described in Details}
  \item{binlevels}{a list containing the levels of the binary variables}
}
\author{Werner A. Stahel, ETH Zurich}
\examples{
data(d.blast)
dd <- datatype(d.blast[1:50,])
attr(dd,"vartype")
dr <- datatype(d.blast[as.numeric(d.blast$location)<=2,])
attr(dr,"vartype")
attr(dr,"binlevels")
}
\keyword{attribute}
