\name{getvars}
\alias{getvars}
\title{Extract Variables for Model Fitting from a Data Frame
}
\description{
  Extracts the variables occurring in a formula, in raw form
  (as \code{get\_all\_vars}) or in transformed form
  (as \code{model.frame} does it).
}
\usage{
getvars(formula, x = NULL, data = NULL, rawvars = TRUE, jitterbinary = TRUE, ...)
}
\arguments{
  \item{formula}{a model 'formula' or 'terms' object or an R object.}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{data}{a data.frame, list or environment (or object coercible by
          'as.data.frame' to a data.frame), containing the variables in
          'formula'.  Neither a matrix nor an array will be accepted.}
  \item{rawvars}{logical. If \code{TRUE}, untransformed variables will
    be extracted, otherwise, variables as transformed in \code{formula}.}
  \item{jitterbinary}{logical. Should binary variables be marked to be
    jittered in plots?}
  \item{\dots}{further arguments such as
    \code{data, weight, subset, offset}
    used to create further columns in the resulting data.frame, with
    parenthesised names such as '"(offset)"'}
}
%%- \details{
%%- %%  ~~ If necessary, more details than the description above ~~
%%- }
\value{
  data.frame containing the extracted variables.

  \code{attr(.,"is.fac")} characterizes the variables coming from
    \code{formula}: 0 for numeric, 1 for factor, 2 for binary
}
\author{
WSt
}

\seealso{
  \code{\link{model.frame}}
}
\examples{
data(d.blast)
getvars(log10(tremor)~log10(distance)*log10(charge), data=d.blast)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{auxiliary}% use one of  RShowDoc("KEYWORDS")

