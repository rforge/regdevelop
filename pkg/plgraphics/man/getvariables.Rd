\name{getvariables}
\alias{getvariables}
\alias{getvarnames}
\title{Extract Variables or Variable Names for Model Fitting
}
\description{
  Extracts the variables names (\code{getvarnames}
  occurring in a formula, in raw form
  (as \code{get\_all\_vars}) or in transformed form
  (as \code{model.frame} does it).
}
\usage{
getvariables(formula, data = NULL, transformed = TRUE,
  envir = parent.frame(), ...)
getvarnames(formula, data = NULL, transformed = FALSE)
}
\arguments{
  \item{formula}{a model 'formula' or 'terms' object or an R object,
  or a character vector of variable names}
  \item{data}{a data.frame, list or environment (or object coercible by
    'as.data.frame' to a data.frame), containing the variables in
    'formula'.  Neither a matrix nor an array will be accepted.}
  \item{transformed}{logical. If \code{TRUE}, variables will be extracted as
    transformed in \code{formula}, otherwise, untransformed variables
    are returned.}
  \item{envir}{environment in which the \code{formula} will be evaluated
    }
  \item{\dots}{further arguments such as
    \code{data, weight, subset, offset}
    used to create further columns in the resulting data.frame, with
    names between dots such as '".offset."'}
}
%%- \details{
%%- %%  ~~ If necessary, more details than the description above ~~
%%- }
\value{
  For \code{getvarnames}: list of character vectors containing
  \item{varnames}{names of all variables (\code{transformed=FALSE})
    or simple terms (\code{transformed=TRUE})}
  \item{xvar}{those from the right hand side of the formula}
  \item{yvar}{left hand side}

  For \code{getvariables}: 
  data.frame containing the extracted variables or simple terms.
}
\author{
Werner A. Stahel
}

\seealso{
  \code{\link{model.frame}}, \code{\link{get_all_vars}}
}
\examples{
data(d.blast)
getvariables(log10(tremor)~log10(distance)*log10(charge), data=d.blast)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{auxiliary}% use one of  RShowDoc("KEYWORDS")

