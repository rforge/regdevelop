\name{modelTable}
\Rdversion{1.1}
\title{Collect results for different regression models}
\alias{modelTable}
\alias{format.modelTable}
\alias{print.modelTable}
\alias{compareTerms}
\alias{[.modelTable}

\description{
  \code{modelTable}
  produces a table that shows which variables are contained in various
  models. Typically, the different models are the results of various model
  selection functions and model development strategies.

  \code{format.modelTable} and \code{print.modelTable} are used to turn
  the results into a nice form.

  \code{compareTerms} returns a table indicating which terms appear in
  which models.
}
\usage{
modelTable(models, data = NULL, seq = NULL)

\method{format}{modelTable}(x, digits=3,
       stars = c("***","** ","*  ",":  ",".  "), sep="", ...)
\method{print}{modelTable}(x, tex = FALSE, ...)

compareTerms(... , list=NULL, seq = NULL)
}

\arguments{
  \item{models}{a character vector containing the names of the models to
    be collected, or a list of model fits. }
  \item{data}{data on which all models are evaluated.}
  \item{seq}{The sequence in which the terms should appear in the table.
    If incomplete, the terms mentioned in \code{seq} will be the first
    rows in the result. }
  \item{x}{object of class \code{modelTable}}
  \item{digits}{to be used in formatting}
  \item{stars}{symbols used to indicate level of significance. They will
    be used for the intervals delimited by
    \code{c(-1, 0.001, 0.01, 0.05, 0.1, 1)}
  }
  \item{sep}{separator. "\&" is used to prepare LaTeX format.}
  \item{tex}{if TRUE, the output is can be used in (La-)TeX.}
  \item{\dots}{other arguments, passed to \code{print} or
    \code{format}.\cr
    In \code{compareTerms}, the models for which the terms should be
    compared can be listed as in \code{c} or \code{list} etc.}
  \item{list}{A list of models.}
}
%% \details{
%%   Note that \code{print.modelTable} allows you to prepare a (La-)TeX
%%   table.
%% }
\value{
  \code{modelTable} produces a list of class \code{"modelTable"},
  with components
  \item{coef}{the coefficients of the models, collected in a matrix.
    The columns of the matrix correspond to the models, the rows, to the
    terms.  For terms with more than 1 degree of freedom, the matrix
    element is NA.  If possible, standardized coefficients are given.}
  \item{p}{p values, in the same form}
  \item{sd.terms}{standard deviations for terms with 1 degree of
    freedom}
  \item{sigma}{estimated error standard deviations of the models}

  \code{format.modelTable} returns a character matrix and thereby
  prepares for \code{print.modelTable}.  If needed, the ouput of
  \code{format.modelTable} can be reordered or subsetted.
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\seealso{
  on formatting and printing of P-values, the standard \R functions
  \code{\link{symnum}}, \code{\link{format.pval}}.
}
\author{Werner A. Stahel, ETH Zurich}

\examples{
data(d.blast)
r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast)
r.bl2 <- update(r.blast, ~.+location:log10(distance)+I(log10(charge)^2) )
r.mt <- modelTable(list(large=r.bl2, original=r.blast))
r.mt

format(r.mt)
print(r.mt, tex=TRUE)
compareTerms(r.blast, r.bl2, seq=c("log10(distance)", "log10(charge)"))
}
\keyword{regression}
