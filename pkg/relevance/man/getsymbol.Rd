\name{getsymbol}
\alias{getsymbol}
\alias{rlvSymbols}
\alias{pSymbols}

\title{Get Symbol According to Categorization of Numerical Value
}
\description{Get the symbol that corresponds to a numerical value
  according to a translation list.
  Used to get the symbols characterizing p-values or relevance values
}
\usage{
getsymbol(value, symbols)

rlvSymbols
pSymbols
}

\arguments{
  \item{value}{tne numerical value (may be a vector)
  }
  \item{symbols}{list with the component \code{cutpoint} and \code{symbol}
  }
}
\value{character (vector) with the required symbols.
  Values in the $i$th interval will result in the $i$th symbol.
%%    symbols$symbol[as.numeric(cut(value, symbols$cutpoint))]
}
\details{
  \code{pSymbols} and \code{rlvSymbols} are lists that define the
  categories (by the component \code{cutpoint}) for
  p-values and relevance values, respectively, and the corresponding
  symbols (component \code{symbol}) to be printed. A third component
  is the legend that is printed if the argument \code{legend} of
  \code{print.termtable} is \code{TRUE}.
}
\author{Werner A. Stahel
}
\examples{
   pSymbols

   getsymbol(c(0.5, 0.02, 0.003, 0), pSymbols)

## [1] " "  "*"  "**" NA  
}
\keyword{manip}
\keyword{utility}
