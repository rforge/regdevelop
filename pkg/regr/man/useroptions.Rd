\name{useroptions}
\alias{useroptions}
\alias{getUseroption}
\alias{.useroptions}
\alias{useroptionsDefault}
\title{Set and Get User "Session" Options that influence "regr0"s behavior}
\description{
  The user can set (and get) "Session" options which influence the
  behavior \pkg{regr0} functions, notably \code{\link{stamp}()}.
}
\usage{
useroptions  (x=NULL, default=NULL, list=NULL, useroptions = NULL,
    assign = TRUE, ...)
getUseroption(x, default = NULL)
.useroptions
useroptionsDefault
}
\arguments{
  \item{x}{option to query or set}
  \item{default}{default value}%% FIXME
  \item{list}{a list of options to be set}%% FIXME
  \item{...}{any Useroptions can be defined or modified,
    using 'name = value'.}
  \item{useroptions}{a list containing the actual user options.
    The effective default is \code{.useroptions}.}
  \item{assign}{logical: should the result be assigned to
    \code{.useroptions}?}
}
\value{For 'getUseroption', the current value set for Useroption 'x',
  or 'NULL' if the option is not set.

  For 'useroptions()', a list of all set Useroptions sorted by name.  
  For 'useroptions(name)', a list of length one containing the set
  value, or 'NULL' if it is not set.
  For uses setting one or more options, a list with the previous values
  of the Useroptions that are changed (returned invisibly).
}
\author{Werner A. Stahel}
\note{In contrast to the \code{\link{options}} of R base,
  the \code{Useroptions} are stored on the \code{.GlobalEnv}
  (if anyone is changed)
  and will be recovered if the latter is stored and loaded again.
}
\seealso{
  \code{\link{stamp}}; \R's own predefined \code{\link{options}()}.
}
\examples{
## see example(stamp)

## set options
useroptions(project="Example A",  step="regression analysis")
useroptions(show.termeffects=FALSE)  ## suppress dummy coefficients
                        ##  for factors when using  print.regr 
useroptions("regr.contrasts")  ##  c("contr.wsum","contr.wpoly")

useroptions(stamp=NULL)
getUseroption("stamp", default=2)
}
\keyword{utilities}
