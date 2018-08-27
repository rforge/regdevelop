\name{fitcomp}
\alias{fitcomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Component Effects for a Model Fit}
\description{
  Determines effects of varying each of the given variables while all
  others are held constant. This function is mainly used to produce
  plots of residuals versus explanatory variables, also showing
  component effects.
}
\usage{
fitcomp(object, data = NULL, vars=NULL, transformed=FALSE, se = FALSE, 
  xm = NULL, xfromdata = FALSE, noexpand=NULL, nxcomp = 51)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a model fit, result of a fitting function}
  \item{data}{data frame in which the variables are found.
    If not provided, it is obtained from \code{object}.}
  \item{vars}{character vector of names of variables for which
    components are required. Only variables that appear in \code{data}
    will be used.}
  \item{transformed}{logical: should components be calculated for 
    transformed explanatory variables? If \code{TRUE}, the variables are
    transformed as implied by the model.}
  \item{se}{if TRUE, standard errors will be returned}
  \item{xm}{named vector of values of the fixed (central) point from
    which the individual variables are varied in turn. \cr
    Defaults to the componentwise median of quantitative variables and
    the modes of factors.}
  \item{xfromdata}{if TRUE, the components effects will be evaluated for
    the data values in \code{data}. Otherwise, the range of each
    numerical variable is filled with \code{nxcomp} equidistant points,
    whereas for factors, all levels are used.
    This is useful for residual plots with component effects.}
  \item{noexpand}{vector determining which variables should not be
    ``filled in'', probably because they are used like factors.
    Either a character vector of variable names or a vector of
    logical or numerical values with names, in which case the names
    corresponding to positive values will be identified.}
  \item{nxcomp}{number of points used for each (quantitative) variable
    if \code{xfromdata} is \code{FALSE}}
}
\details{
The component effect is defined as the curve of fitted values
obtained by varying the explanatory variable, keeping all the other
variables at their "central value" (the mean of continuous variables
and the mode of factors).
}
\value{A list consisting of
  \item{comp}{component effects}
  \item{x}{the values of the x variables for which the effects have been
    calculated}
  \item{xm}{the values of the x variables that are held fixed while one
    of the variables is varied}
  \item{se}{standard errors of the component effects, if required by the
    argument \code{se}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Werner A. Stahel, ETH Zurich}
%%- \note{ ~~further notes~~ 
%%- 
%%-  ~Make other sections like Warning with \section{Warning }{....} ~
%%- }
\seealso{\code{\link{predict}}}
\examples{
data(d.blast)
t.r <- lm(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
t.fc <- fitcomp(t.r,se=TRUE)
t.fc$comp[1:10,]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
