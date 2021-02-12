\name{plres2x}
\alias{plres2x}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Residuals vs. Two Explanatory Variables}
\description{
  Plot 2 variables, showing a third one with line symbols. Most suitable
  for showing residuals of a model as this third variable.
}
\usage{
plres2x(formula = NULL, reg = NULL, data = NULL, restrict = NULL,
  size = 1, xlab = NULL, ylab = NULL, plextext = NULL, pale = 0.2,
  plargs = NULL, ploptions = NULL, assign = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula of the form \code{~x+y}, where
    \code{x, y} are the 2 variables shown by the coordinates of points,
    and residuals are shown by line symbols: their orientation
    corresponds to the sign of \code{the residual}, and their length,
    to the absolute value.
  }
  \item{reg}{the result of the model fit, from which the residuals are
    extracted
  }
  \item{data}{the data.frame where the variables are found. Only needed
    if the variable 'x' or 'y' is not available from the fitting results.}
  \item{restrict}{absolute value which truncates the size.
    if \code{TRUE}, the inner plotting limits of the residuals is used
    if available.
    Truncation is shown by stars at the end of the line symbols.
  }
  \item{size}{the symbols are scaled so that \code{size*par("cin")[1]}
    is the length of the largest symbol, as a percentage of the length
    of the horizontal axis.}
  \item{xlab, ylab}{labels for horizontal and vertical axes.
    Default to the variable names (or labels)}.
  \item{plextext}{scalar determining the proportion of the plotting area
    that is used as an extension to accomodate the symbols representing
    the residuals}.
  \item{pale}{scalar between 0 and 1: The points are shown in a more
    pale color than the segments as determined by
    \code{\link{colorpale}} with argument \code{pale=pale}.}
  \item{plargs}{result of calling \code{pl.control}.
    If \code{NULL}, \code{pl.control} will be called to generate it.
    If not null, arguments given in \code{...} will be ignored.
  }
  \item{ploptions}{list of pl options.
  }
  \item{assign}{logical: Should the plargs be stored 
    in the \code{pl.envir} environment?
  }
  \item{\dots}{further arguments, passed to \code{plotregr.control}}
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
  t.r <- lm(log10(tremor)~location+log10(distance)+log10(charge),
            data=d.blast)
  plres2x(~distance+charge, t.r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{regression}% __ONLY ONE__ keyword per line
