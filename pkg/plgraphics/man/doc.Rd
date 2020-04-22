\name{doc}
\alias{doc}
\alias{doc<-}
\alias{tit}
\alias{tit<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Define and obtain the doc or tit attribute }
\description{
  The attributes \code{doc} and \code{tit} describe an object, typically
  a data frame or a model. \code{tit} should be a short description (title),
  \code{doc} should contain all documentation useful to identify
  the origin and the changes made to the object.\cr
  The \code{doc} and \code{tit} functions set them and extract these
  attributes. 
}
\usage{
doc(x)
tit(x) 
doc(x) <- value
tit(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object to which the \code{doc} or \code{tit} attribute
    should be attached
    or from which it is obtained}
  \item{value}{character vector (\code{doc}) or string (\code{tit})
    to be stored}
}
\details{
  Plotting and printing functions may search for the \code{tit}
  attribute or even for the \code{doc} attribute, depending on
  \code{c.env$docout}.

  \code{doc(x) <- text} will append the existing \code{doc(x)} text to
  the new \code{text} unless its first element equals (the first element
  of) \code{text}.
  (This avoids piling up the same line by unintended multiple call to
  \code{doc(x) <- value} with the same \code{value}.)
  If the first element of \code{text} equals \code{"^"},
  the first element of \code{doc(x)} is dropped.
  \code{tit(x) <- string} replaces \code{tit(x)} with \code{string}. 
}
\value{
  \code{doc} and \code{tit} return the respective attributes of object
  \code{x} 
}
\author{Werner A. Stahel, ETH Zurich}
\examples{
data(d.blast)
doc(d.blast)
doc(d.blast) <- "I will use this dataset in class soon."
doc(d.blast)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{attribute}
