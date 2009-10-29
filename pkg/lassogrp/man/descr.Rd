\name{descr}
\alias{descr}
\alias{descr<-}
\alias{tit}
\alias{tit<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Define and obtain the descr or tit attribute }
\description{
  The attributes \code{descr} and \code{tit} describe an object, typically
  a data frame or a model. \code{tit} should be a short description (title),
  \code{descr} should contain all documentation useful to identify
  the origin and the changes made to the object.\cr
  The \code{descr} and \code{tit} functions set them and extract these
  attributes. 
}
\usage{
descr(x)
tit(x) 
descr(x) <- value
tit(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ object to which the \code{descr} or \code{tit} attribute
    should be attached
    or from which it is obtained}
  \item{value}{character vector (\code{descr}) or string (\code{tit})
    to be stored}
}
\details{
  Plotting and printing functions may search for the \code{tit}
  attribute or even for the \code{descr} attribute, depending on
  \code{c.env$docout}.

  \code{descr(x) <- text} will append the existing \code{descr(x)} text to
  the new one unless the first element of text equals \code{"^"},
  whereas \code{tit(x) <- string} replaces \code{tit(x)}. 
}
\value{
  \code{descr} and \code{tit} return the respective attributes of object
  \code{x} 
}
\author{Werner A. Stahel, ETH Zurich}
\examples{
data(asphalt)
tit(asphalt)
descr(asphalt)
descr(asphalt) <- "I will use this dataset in class soon."
descr(asphalt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{attribute}
