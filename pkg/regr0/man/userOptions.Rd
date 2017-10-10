\name{userOptions}
\alias{userOptions}
\alias{getUserOption}
\alias{UserOptions}
\alias{UserDefault}
\title{Set and Get User "Session" Options that influence "regr0"s behavior}
\description{
  The user can set (and get) "Session" options which influence the
  behavior \pkg{regr0} functions, notably \code{\link{stamp}()}.
}
\usage{
userOptions  (x=NULL, default=NULL, list=NULL, ...)
getUserOption(x, default = NULL)
UserOptions
UserDefault
}
\arguments{
  \item{x}{option to query or set}
  \item{default}{default value}%% FIXME
  \item{list}{a list of options to be set}%% FIXME
  \item{...}{any UserOptions can be defined or modified,
    using 'name = value'.}
}
\value{For 'getUserOption', the current value set for UserOption 'x',
  or 'NULL' if the option is not set.

  For 'userOptions()', a list of all set UserOptions sorted by name.  
  For 'userOptions(name)', a list of length one containing the set
  value, or 'NULL' if it is not set.
  For uses setting one or more options, a list with the previous values
  of the UserOptions that are changed (returned invisibly).
}
\author{Werner A. Stahel}
\note{In contrast to the \code{\link{options}} of R base,
  the \code{UserOptions} are stored on the \code{.GlobalEnv}
  (if anyone is changed)
  and will be recovered if the latter is stored and loaded again.
}
\seealso{
  \code{\link{stamp}}; \R's own predefined \code{\link{options}()}.
}
\examples{
## see example(stamp)
## set options
userOptions(project="Example A",  step="regression analysis")
userOptions(show.dummycoef=FALSE)  ## suppress dummy coefficients
                        ##  for factors when using  print.regr 
userOptions("regr.contrasts")  ##  c("contr.wsum","contr.wpoly")

userOptions(stamp=NULL)
getUserOption("stamp", default=2)
}
\keyword{utilities}
