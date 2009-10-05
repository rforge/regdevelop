\name{modelTable}
\Rdversion{1.1}
\alias{modelTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Collect results for different regression models
}
\description{
Produces a table that shows which variables are contained in various
models. Typically, the different models are the results of various model
selection functions and model development strategies.
}
\usage{
modelTable(models, data = NULL, seq = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{a character vector containing the names of the models to
    be collected, or a list of model fits
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{seq}{The sequence in which the terms should appear in the table
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{a list of class 'modelTable', consisting of
  \item{coef}{the coefficients of the models, collected in a matrix.
    The columns of the matrix correspond to the models, the rows, to the
    terms. For terms with more than 1 degree of freedom, the matrix
    element is NA.}
  \item{p}{p values, in the same form}
  \item{sd.terms}{standard deviations for terms with 1 degree of
    freedom}
  \item{sigma}{estimated error standard deviations of the models}
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{Werner A. Stahel, ETH Zurich}

\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
