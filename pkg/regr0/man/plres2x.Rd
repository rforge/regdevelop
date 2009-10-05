\name{plres2x}
\alias{plres2x}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Residuals vs. Two Explanatory Variables}
\description{
  Plot 2 variables, showing a third one with line symbols. Most suitable
  for showing residuals of a model as this third variable.
}
\usage{
plres2x(formula = NULL, reg = NULL, data = reg, restricted = NULL,
size = 0, slwd = 1, scol = 2, xlab = NULL, ylab = NULL,
xlim = NULL, ylim = NULL, main = NULL, cex.title = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{z~x+y}, where
    \code{x, y} are the 2 variables shown by the coordinates of points,
    and z is shown by line symbols: their orientation corresponds
    to the sign of \code{z}, and their length, to the absolute value.
  }
  \item{reg}{the result of the model fit}
  \item{data}{the data.frame where the variables are found}
  \item{restricted}{absolute value which truncates the size. Truncation
    is shown by stars at the end of the line symbols.
    }
  \item{size}{the symbols are scaled so that \code{size} is the length of
    the largest symbol, as a percentage of the length
    of the horizontal axis. Defaults to \code{5/log10(n)}, where
    \code{n} is the number of observations}
  \item{slwd}{line width of the line symbols}
  \item{scol}{color of the symbols}
  \item{xlab, ylab}{labels for horizontal and vertical axes.
    Default to the variable names}
  \item{xlim, ylim}{plot ranges for horizontal and vertical axes.
    They are expanded to accomodate the symbols}
  \item{main}{main title of plot.
    Defaults to \code{formula}}
  \item{cex.title}{character expansion for the main title}
  \item{\dots}{further arguments, passed to \code{plot}}
}
%%- \details{
%%-   ~~ If necessary, more details than the description above ~~
%%- }
\value{
  none.
}
\author{Werner A. Stahel and Andreas Ruckstuhl}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
  data(d.blast)
  t.r <- regr(log10(tremor)~location+log10(distance)+log10(charge),
            data=d.blast)
  plres2x(~distance+charge, t.r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{regression}% __ONLY ONE__ keyword per line
