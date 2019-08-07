\name{plmark}
\alias{plmark}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Labels for Extreme Points
}
\description{
Determine extreme points and get labels for them.
}
\usage{
plmark(x, y = NULL, markextremes = NULL, plabel = NULL, plargs = NULL, ploptions = NULL)
}

\arguments{
  \item{x,y}{coordinates of points. If \code{x} is of length 0,
    it is retrieved from \code{plargs$pldata[,1]}.
  }
  \item{markextremes}{proportion of extreme points to be 'marked'.
    This may be a vector or list of proportions with names
    indicating the variables for which the proportion is to be applied.
    In the default case (\code{NULL}), it is obtained from
    \code{ploptions}, which in turn leads to calling the function
    \code{\link{markextremes}} with the argument equal to the number of
    (finite) observations.    
  }
  \item{plabel}{character vector of labels to be used for extreme
    points. If \code{NULL}, they are obtained from \code{plargs$plabel}.
  }
  \item{plargs, ploptions}{result of \code{\link{pl.control}}, cf
    \code{\link{plpoints}}
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A character vector in which the 'marked' observations contain
  the respective label and the others equal \code{""}.
}
\author{Werner A. Stahel
}

\seealso{\code{\link{plyx}}
}
\examples{
  plyx(Sepal.Width ~ Sepal.Length, data=iris)
  ( t.plab <-
    plmark(iris$Sepal.Length, iris$Sepal.Width, markextremes=0.03) )
}
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
