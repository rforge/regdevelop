\name{shortenstring}
\alias{shortenstring}

\title{Shorten Strings
}
\description{Strings are shortened if they are longer than
  \code{n} 

}
\usage{
shortenstring(x, n = 50, endstring = "..", endchars = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a string or a vector of strings
  }
  \item{n}{maximal character length
  }
  \item{endstring}{string(s) to be appended to the shortened strings
  }
  \item{endchars}{number of last characters to be shown at the end of
    the abbreviated string. By default, it adjusts to \code{n}.
  }
}
\value{Abbreviated string(s)
}
\author{Werner A. Stahel
}
\seealso{\code{\link{substring}}, \code{\link{abbreviate}}
}
\examples{
shortenstring("abcdefghiklmnop", 8)

shortenstring(c("aaaaaaaaaaaaaaaaaaaaaa","bbbbc",
  "This text is certainly too long, don't you think?"),c(8,3,20))

}
\keyword{utilities}

