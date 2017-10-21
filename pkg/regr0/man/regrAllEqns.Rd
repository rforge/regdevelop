\name{regrAllEqns}
\alias{regrAllEqns}
\alias{regrAllEqnsXtr}
\alias{print.regrAllEqns}
\alias{plot.regrAllEqns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{All Equations of a Linear Model
}
\description{Selection of best linear models by exhaustive search
}
\usage{
regrAllEqns(formula, data, weights = NULL, nbest = 50, nvmax = 20,
  force.in = NULL, force.out = NULL, codes = NULL, really.big = FALSE,
  ...)
regrAllEqnsXtr(object, nbest=1, criterion="cp")
\method{print}{regrAllEqns}(x, nbest = 20, criterion = "cp",
  printcriteria = FALSE, printcodes = TRUE, ...) 
\method{plot}{regrAllEqns}(x, criterion = "cp", critrange = 10,
  minnumber = 10, nbest=10, codes = x$codes, ncharhorizontal = 6,
  col="blue", legend = TRUE, mar=6, main="", cex=0.7*par("cex"),
  cex.lab=1.2*par("cex.lab"), ...) 
}

\arguments{
\item{formula}{formula for the full model
}
\item{data}{data.frame where the variables are found
}
\item{weights}{weights for weighted LS estimation
}
\item{nbest}{number of equations to record for each size.
  For \code{print.regrAllEqns}, number of equations printed in total.
}
\item{nvmax}{maximum size of subsets to examine
}
\item{force.in}{formula for terms to be kept in all models
}
\item{force.out}{formula for terms to be excluded from all models
}
\item{codes}{codes for the terms to be used for identifying the
  equations in the result 'outmat' and in the plot. Preferably one
  letter each.
}
\item{really.big}{Must be TRUE to perform exhaustive search on more
  than 50 variables.
}
\item{criterion}{Criterion to be used to rank the equations and for the
  vertical axis in plotting.
}
\item{object, x}{result of \code{regrAllEqns} to be used by
  \code{regrAllEqnsXrt} or \code{print} or \code{plot}
}
\item{printcriteria}{logicel: should the table of criteria be printed?
}
\item{printcodes}{logical: should the short codes for the terms
  (columns in the primary output 'outmat') be listed?
}
\item{critrange}{range of criterion that determines how many equations
  will be used in plotting: those that are less than \code{critrange}
  worse than the best.
}
\item{minnumber}{minimum number of equations shown; overrides
  \code{critrange} if needed.
}
\item{ncharhorizontal}{maximum length of string on margin 3 to be shown
  horizontally. Use \code{ncharhorizontal=0} to get all codes
  shown vertically. This avoids overwriting.
}
\item{col}{color used for the equation identifiers in the plot
} 
\item{legend}{logical: Should a list of short and long codes for terms
  be shown as a legend? Alternatively, a position for the legend:
  keyword like \code{"bottomleft"} or 2 coordinates.
}
\item{mar}{margins: either all 4 margins of the plot or 1 number for the
  top margin.
}
\item{main}{main title
}
\item{cex, cex.lab}{character expansion for codes in the plot and on the
  upper margin}
\item{\dots }{Other arguments, to be passed to \code{\link{regsubsets}},
  \code{print} or \code{plot}, respectively.
}
}
\details{This is a wrapper function for function
  \code{\link{regsubsets}} od package \code{leaps}.
  When factors (and other multicolumn terms) occur in the full model,
  it returns models that either include the factor or exclude it,
  and avoids those that contain some, but not all dummy variables
  related to it.

  \code{regrAllEqns} also provides convenient information for the 
  the plotting method.
  The results of \code{\link{regsubsets}} and its \code{summary}, which
  include the "avoided" eaquations, are also contained in the value. 
}
\value{\code{regrAllEqns} returns a list, containing
\item{which}{}
\item{criteria}{}
\item{codes}{}
\item{force.in}{}
\item{force.out}{}
\item{outmat}{}
\item{allsubsets}{}
\item{obj}{see \code{?\link{regsubsets}}}
\code{regrAllEqnsXtr} returns a formula or a list of formulas
(if \code{nbest>1}).
\code{plot.regrAllEqns} invisibly returns the coordinates
of the models used in the plot.
}
\author{Werner A. Stahel
}
\note{For further details and notes, see \code{?\link{regsubsets}}
}
\seealso{
\code{\link{regsubsets}}, \code{\link{leaps}}, \code{\link{step}}, 
}
\examples{
data(d.birthrates)
t.formula <- fertility ~ catholic + single24 + single49 + eAgric + eIndustry +
eCommerce + eTransport + eAdmin + gradeHigh + gradeLow + educHigh +
bornLocal + bornForeign + altitude + language
## altitude and language are factors

r.ae <- regrAllEqns(t.formula, data=d.birthrates, nbest=100, really.big=TRUE)
print(r.ae)
plot(r.ae, mar=9, ncharhorizontal=0, main="birthrates example",
     legend="bottomright", xlim=c(4,20))

## extract the formula of the best model
( t.formula <- regrAllEqnsXtr(r.ae) )
## ... and fit it
regr(t.formula, data=d.birthrates)

}
\keyword{regression}

