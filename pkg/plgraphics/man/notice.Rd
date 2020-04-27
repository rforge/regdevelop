\name{notice}
\alias{notice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a Notice
}
\description{Generate a notice to be sent to output
}
\usage{
notice(..., printnotices = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{contents of the notice, will be pasted together
  }
  \item{printnotices}{logical: Should the notice be printed?
    Default is the respective pl option. 
}
}
\details{
  This function is very similar to 'message'
}
\value{None.
}
\author{Werner A. Stahel
}

\seealso{
\code{\link{message}}
}
\examples{
  ff <- function(x) {
    if (length(x)==0)  {
      notice("ff: argument 'x' is NULL. I return 0")
      return(0)
    }
    1/x
  }
  ff(3)
  ff(NULL)
  oo <- ploptions(printnotices=FALSE)
  ff(NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}% use one of  RShowDoc("KEYWORDS")

