\name{regroptions}
\alias{regroptions}
\alias{getRegrOption}
\alias{.regroptions}
\alias{regroptionsDefault}
\title{Set and Get Regr Options that Influence \code{regr}'s behavior}
\description{
  The user can set (and get) options which influence the
  behavior of \pkg{regr} functions.
}
\usage{
regroptions(x=NULL, list=NULL, default=NULL, regroptions = NULL,
    assign = TRUE, ...)
getRegrOption(x, regroptions=NULL)
.regroptions
regroptionsDefault
}
\arguments{
  \item{x}{character (vector): name(s) of regroption(s) to query
      If \code{x} is set, all further arguments will be ignored.
    }
  \item{list}{a named list of options to be set, see Details}
  \item{default}{character vector of option names.
    These regroptions will be set according to \code{regroptionsDefault}.
    \code{default="all"} or \code{=TRUE} will reset all options.
    If \code{default} is set, all further arguments will be ignored.
  }
  \item{regroptions}{the list of options that should be modified.
    The effective default is \code{.regroptions}
  }
  \item{assign}{logical: should the list be assigned to
    \code{.regroptions}?}     
  \item{...}{any regroptions can be defined or modified,
    using \code{name = value}, as in \code{\link{options}} of basic R.}
}
\value{For \code{regroptions(x)}, where \code{x} is the name of a
  regr option, the current value of the option.
  or \code{NULL} if it is not such a name.
  If \code{x} contains several (valid) names, the respective list.
  
  For \code{regroptions()}, the list of all regr options sorted by name.  
  
  For uses setting one or more options, the important effect is a changed
  list \code{.regroptions} that is used by the packages functions
  (if \code{assign} is \code{TRUE}).
  The (invisibly) returned value is the same list, complemented by an
  attribute \code{"old"} containing the previous values of those options
  that have been changed.
  This list is useful for undoing the changes, usually at the end of
  graphical functions, to restore the previous status.

  \code{.regroptions} and \code{regroptionsDefault} are lists of options.
  The first one is generated whenever regroptions are set
  (and \code{assign} is \code{TRUE}), and is then
  stored in the \code{".GlobalEnv"} environment.
  The second is stored in the package \code{regr}.
}
\details{
  If the argument \code{list} is set, it must be a named list,
  and each component, with name \code{name} and value \code{value}
  is used as described for arguments in \code{...}, see above
  (in addition to such arguments).

  For a list of regr options, consult the vignette.
}
\author{Werner A. Stahel}
\note{In contrast to the \code{\link{options}} of R base,
  the \code{.regroptions} list is stored on the \code{.GlobalEnv}
  (if any component is changed)
  and will be recovered if the latter is stored and loaded again.
}
\seealso{
  \code{\link{stamp}}; \R's own predefined \code{\link{options}()}.
}
\examples{
## get options
regroptions("regr.contrasts")  ##  c("contr.wsum","contr.wpoly")

## set options
regroptions(digits=3,
  termcolumns=c("coef", "stcoef", "df", "signif", "p.symb"))
regroptions("digits")

regroptions(default="all") ## restore 
regroptions("digits")
t.opt <- regroptions(digits=6, assign=FALSE)
t.opt$digits
regroptions("digits") ## still the default

}
\keyword{utilities}
