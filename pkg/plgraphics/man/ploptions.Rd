\name{ploptions}
\alias{ploptions}
\alias{.ploptions}
\alias{ploptionsDefault}
\title{Set and Get User "Session" Options that influence "regr0"s behavior}
\description{
  The user can set (and get) "Session" options which influence the
  behavior \pkg{regr0} functions, notably \code{\link{stamp}()}.
}
\usage{
ploptions(x = NULL, default = NULL, list = NULL, ploptions = NULL,
          assign = TRUE,...)
.ploptions
ploptionsDefault
}
\arguments{
  \item{x}{character (vector) of name(s) of ploptions.
%%-     option to query (if character) or set (if a named object).
%%-     It may be a character vector of options if it is a query.\cr
%%-    For those ploptions that are lists, \code{x} must be a list, too.
%%-     It contains the components that should be modified, see the second
%%-     example below.
    If \code{x} is set, all further arguments will be ignored.
  }
  \item{list}{a named list of options to be set}
  \item{default}{character vector of option names.
    The named options will be set according to \code{ploptionsDefault}.
    \code{default="all"} or \code{=TRUE} will reset all options.
    If \code{default} is set and \code{x} is not, all further arguments
    will be ignored.
  }
  \item{ploptions}{the list of options that should be modified}
  \item{assign}{logical: should the list be assigned to
    \code{.ploptions}?}    
  \item{...}{any ploptions can be defined or modified,
    using \code{name = value}, as in \code{\link{options}} of basic R.}
}
\value{For \code{ploptions()}, a list of all set UserOptions sorted by name.  

  For \code{ploptions(x)}, where \code{x} is the name of a pl option,
  a \emph{list} of length one containing the option,
  or \code{NULL} if it is not such a name.

  For uses setting one or more options, the important effect is a changed
  list \code{.ploptions} that is used by the \code{pl} graphical
  functions (if \code{assign} is \code{TRUE}).
  The (invisibly) returned value is the same list, complemented by an
  attribute \code{"old"} containing the previous values of those ploptions
  that have been changed.
  This list is useful for undoing the changes, usually at the end of
  graphical functions, to restore the previous status.

  \code{.ploptions} and \code{ploptionsDefault} are lists of options.
  The first one is generated whenever ploptions are set
  (and \code{assign} is \code{TRUE}), and is then
  stored in the \code{".GlobalEnv"} environment.
  The second is stored in the package \code{plgraphics}.
}
\details{
  Here is a list of the components of \code{ploptionsDefault},
  describing the suitable alternative values to be set by calling
  \code{ploptions}:
  \cr very incomplete for the time being!!!
  %% ===========================================================
  \itemize{
    \item{tickintervals}{vector of length 2. The first element is the
      desired number of tick intervals for axes, to be used as argument
      \code{n} in \code{\link{pretty}}.
      The second determines how many tick labels are shown in the same
      way, and should therefore be smaller than (or equal to) the first.
    }      
    \item{gridlines}{can be \cr
      -- a logical indicating if gridlines should be drawn. If
      \code{TRUE}, gridlines will be drawn at the values given in
      \code{attr(.,"ticksat")};
      -- a vector of values at which the gridlines should appear;\cr
      -- a list of length 2 of such values;\cr
      -- a named list. If a name equals the attribute \code{varname}
      of either the x or y variable, the respective component will be
      used.
    }
    \item{\code{smoothline.lty, smoothline.col}}{line type and color.
      Note that if there is a \code{smooth.group} factor,
      \code{group.lty} and \code{group.col} are used.
    }
    \item{\code{smoothline.lwd}}{line width. If of length 2 (or more),
      the second element is the factor by which the line width is
      reduced for simulated smooths (that is, for the second to the last
      column of \code{smoothline$y}). It defaults to 0.7.
    }
    \item{smooth.xtrim}{proportion of fitted values to be trimmed off on
      both sides when drawing a smooth line, either a number or a function
      that takes the number of points as its argument.
      The default is the simple function \code{2^log10(n)/n}.
    }
    \item{condquant...}{Conditional quantiles for censored residuals.
      \itemize{
	\item{condquant}{logical: should bars be drawn for censored residuals?
	  If \code{FALSE}, censored observations will be set to the median of
	  the conditional distribution and shown by a different plotting
	  character, see argument \code{censored} of
	  \code{\link{ploptions}}.
	  If \code{NULL}, the standard plotting character will be used.
	}
	\item{condquant.probrange}{range for probabilities.
	  If the probability corresponding to the censored part of the
	  distribution is outside the range, bars will not be drawn.}
	\item{condquant.pale}{factor by which the \code{pcol} color
	  will be paled to show the points (\code{condquant.pale[1]})
	  and the bars (\code{...[2]}).
	}
      }
    }
  }
}
\author{Werner A. Stahel}
\note{In contrast to the \code{\link{options}} of R base,
  the \code{plptions} are stored in the \code{.GlobalEnv}
  (if anyone is changed)
  and will be recovered if the latter is stored and loaded again.
}
\seealso{
  \code{\link{stamp}}; \R's own predefined \code{\link{options}()}.
}
\examples{
## see example(stamp)
ploptions("stamp")
ploptions(c("jitterFactor", "gridlines"))

## set options
ploptions(stamp=FALSE, pch=0, col=c.colors[-1], anything="do what you want")
ploptions(c("stamp", "anything"))
ploptions(default=TRUE)
t.plopt <- ploptions(smoothline.col="green", assign=FALSE)
t.plopt$smoothline.col
ploptions("smoothline.col") ## unchanged
}
\keyword{utilities}
