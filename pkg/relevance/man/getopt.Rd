\name{getopt}
\alias{getopt}
\title{Get Option that Influence the Output of Relevance Functions}
\description{
  Get options that set relevance thresholds and select columns
  of output tables to be shown as well as other items influencing output
}
\usage{
  getopt(x, opt = NULL)
}
\arguments{
  \item{x}{character (vector): name(s) of option(s) to query
  }
  \item{opt}{the list of options from which the option should be obtained.
    The effective default is a list of defaults obtained by
    \code{relevance:::optdefault}
  }
}
\value{The current value of the option (\code{NULL} if such an option is
  not defined). If \code{x} contains several (valid) names,
  the respective list is returned.
}
\note{options can be changed by setting them through the standard
  \code{options} function, see Examples.
}
\author{Werner A. Stahel}
\note{In contrast to the \code{\link{options}} of R base,
  the \code{.regroptions} list is stored on the \code{.GlobalEnv}
  (if any component is changed)
  and will be recovered if the latter is stored and loaded again.
}
\seealso{
  \code{\link{options}()}
}
\examples{
getopt("na.print")

## set options
options(digits=3, termcolumns=c("coef", "stcoef", "df", "signif", "p.symb"))
getopt("termcolumns") ## same as getOption("termcolumns")
}
\keyword{utilities}
