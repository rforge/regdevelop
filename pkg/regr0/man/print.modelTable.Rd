\name{print.modelTable}
\Rdversion{1.1}
\alias{print.modelTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print a modelTable Object}
\description{Prints a modelTable object. These collect results 
for >1 regression models}
\usage{
print.modelTable(x, tex = FALSE, digits = 3,
  stars = c("***", "** ", "*  ", ":  ", ".  "), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a modelTable object}
  \item{tex}{logical. If TRUE, the output will be suitable for pasting
    into (la)tex
  }
  \item{digits}{number of digits used for printing
  }
  \item{stars}{symbols used for showing significance levels
  }
  \item{\dots}{other arguments, passed to \code{print}
  }
}
\value{formatted table that is printed (invisible)
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{
\code{\link{modelTable}}
}
\examples{
  example(modelTable)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
