\name{terms2order}
\alias{terms2order}
\title{Second Order Terms of a Formula}
\description{Create all square and interaction terms of a given formula
}
\usage{
terms2order(object, squared = TRUE, interactions = TRUE)
}
\arguments{
  \item{object}{An object for which \code{formula(object)} returns a
    formula. Usually a regression fit}
  \item{squared}{logical: Should square terms be added to the
    interactions?}
  \item{interactions}{logical: Should interaction terms be added?}  
}
\details{The 'squares' are the squares of the terms of \code{dataClass}
  'numeric' appearing in the formula (cf.
  \code{attr(object$terms, "dataClass")},
  except for terms that are already squares (contain \code{^}).
  Analogously for interactions: they combine all terms that do not
  contain \code{:}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{add1.regr}}}

\value{'term.labels' of the generated terms}

\author{Werner A. Stahel, ETH Zurich}
\examples{
  data(d.blast)
  t.r <- regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
  terms2order(t.r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
