\name{ploptions}
\alias{ploptions}
\alias{plpar}
\alias{.ploptions}
\alias{ploptionsDefault}
\title{Set and Get User "Session" Options that influence "regr0"s behavior}
\description{
  The user can set (and get) options -- mostly graphical "parameters" --
  which influence the behavior \pkg{plgraphics} functions. %%, notably \code{\link{stamp}()}.
}
\usage{
ploptions(x = NULL, list = NULL, default = NULL, ploptions = NULL,
          assign = TRUE, setpar = FALSE, ...)
plpar(x=NULL, list=NULL, default=NULL, ploptions = NULL,
      assign = FALSE, setpar = TRUE, ...)
.ploptions
ploptionsDefault
}
\arguments{
  \item{x}{character (vector) of name(s) of ploptions to query.
%%-     option to query (if character) or set (if a named object).
%%-     It may be a character vector of options if it is a query.\cr
%%-    For those ploptions that are lists, \code{x} must be a list, too.
%%-     It contains the components that should be modified, see the second
%%-     example below.
    If \code{x} is set, all further arguments will be ignored.
  }
  \item{list}{a named list of options to be set, see Details}
  \item{default}{character vector of option names.
    These ploptions will be set according to \code{ploptionsDefault}.
    \code{default="all"} or \code{=TRUE} will reset all options.
    If \code{default} is set, all further arguments will be ignored.
  }
  \item{ploptions}{the list of options that should be modified}
  \item{assign}{logical: should the list be assigned to
    \code{.ploptions}?}
  \item{setpar}{logical: should the 'margin parameters'
    \code{("mar","oma","mgp","cex")} be set permanently by calling
    \code{\link{par}}? If \code{FALSE}, these options will only be
    active while using pl functions.
  }
  \item{...}{any ploptions can be defined or modified,
    using \code{name = value}, as in \code{\link{options}} of basic R.}
}
\value{For \code{ploptions(x)}, where \code{x} is the name of a pl option,
  the current value of the option.
  or \code{NULL} if it is not such a name.
  If \code{x} contains several (valid) names, the respective list.
  
  For \code{ploptions()}, the list of all plptions sorted by name.  
  
  For uses setting one or more options, the important effect is a changed
  list \code{.ploptions} that is used by the packages functions
  (if \code{assign} is \code{TRUE}).
  The (invisibly) returned value is the same list, complemented by an
  attribute \code{"old"} containing the previous values of those options
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
  If the argument \code{list} is set, it must be a named list,
  and each component, with name \code{name} and value \code{value}
  is used as described for arguments in \code{...}, see above
  (in addition to such arguments).
  
  Here is a list of the components of \code{ploptionsDefault},
  describing the suitable alternative values to be set by calling
  \code{ploptions}:
  \cr incomplete for the time being!!!
  %% ===========================================================
  \itemize{
    \item{keeppar:}{ logical. If TRUE, the graphical parameter settings 
      "mar", "oma", "cex", "mgp", and "mfg" will be maintained when
      leaving high level pl functions, 
      otherwise, the old values will be restored (default).}
    \item{linewidth:}{ vector of \code{lwd} values to be used for
      the different line types (\code{lty}). The package
      sets \code{lwd} to a value
      \code{ploptions("linewidth")[lty]*lwd} intending to balance
      the visual impact of the different line types, e.g.,
      to allow a dotted line to make a similar impression as a solid
      line.}
    \item{ticklength:}{ vector of 4 scalars: tickmark length,
      corresponding to \code{par("tcl")}. The first 2 elements
      define the length of the regular tickmarks, the other two,
      of the ``small'' tichmarks given by
      \code{attr(ticksat, "small")} (\code{ticksat} is a possible
      attribute of each variable).
      There are two elements each in order to define tickmarks that
      cross the axis.
    }
    \item{tickintervals:}{ vector of length 2. The first element is the
      desired number of tick intervals for axes, to be used as argument
      \code{n} in \code{\link{pretty}}.
      The second determines how many tick labels are shown in the same
      way, and should therefore be smaller than (or equal to) the first.
    }
    \item{gridlines:}{  can be \cr
      -- a logical indicating if gridlines should be drawn. If
      \code{TRUE}, gridlines will be drawn at the values given in
      \code{attr(.,"ticksat")};
      -- a vector of values at which the gridlines should appear;\cr
      -- a list of length 2 of such values;\cr
      -- a named list. If a name equals the attribute \code{varname}
      of either the x or y variable, the respective component will be
      used.
    }
    \item{\code{smooth.lty, smooth.col}:}{ line type and color.
      Note that if there is a \code{smooth.group} factor,
      \code{group.lty} and \code{group.col} are used.
    }
    \item{\code{smooth.lwd}:}{ line width. If of length 2 (or more),
      the second element is the factor by which the line width is
      reduced for simulated smooths (that is, for the second to the last
      column of \code{smoothline$y}). It defaults to 0.7.
    }
    \item{smooth.xtrim:}{ proportion of fitted values to be trimmed off on
      both sides when drawing a smooth line, either a number or a function
      that takes the number of points as its argument.
      The default is the simple function \code{2^log10(n)/n}.
      The smoothing function may produce an attribute \code{xtrim}
      that is used as an additional factor to \code{smooth.xtrim}.
      This is applied, e.g., to suppress trimming if a straight line is
      fitted instead of a smooth by requiring \code{smoothLm} as the
      smoothing function.
    }
    \item{smooth.minobs:}{ minimal number of observations needed for
      calculating a smooth.}
    \item{smooth.band:}{ Indicator (logical) determining
      whether "low" and "high" smooth lines should be drawn.
      See above for their definition.}
    \item{condquant...:}{ Conditional quantiles for censored residuals.
      \itemize{
	\item{condquant:}{ logical: should bars be drawn for censored residuals?
	  If \code{FALSE}, censored observations will be set to the median of
	  the conditional distribution and shown by a different plotting
	  character, see argument \code{censored} of
	  \code{\link{ploptions}}.
	  If \code{NULL}, the standard plotting character will be used.
	}
	\item{condquant.probrange:}{ range for probabilities.
	  If the probability corresponding to the censored part of the
	  distribution is outside the range, bars will not be drawn.}
	\item{condquant.pale:}{ factor by which the \code{pcol} color
	  will be paled to show the points (\code{condquant.pale[1]})
	  and the bars (\code{...[2]}).
	}
      }
    }
    \item{plcond...:}{ features of \code{\link{plcond}}.
      \itemize{
	\item{plcond/panel:}{ panel function to be used}
	\item{plcond.nintervals:}{ number of intervals into which
	  numerical variables will be cut}
	\item{plcond.extend:}{ proportion of neighboring intervals for
	  which points are shown. 0 means no overlap.}
	\item{plcond.col:}{ 4 colors to be used to mark the points of
	  the neighboing intervals:
	  The first and second ones color the points lower or higher
	  than the interval of the horizontal conditioning variable,
	  and the other two regulating the same features for the
	  vertical variable.
	  The points which are outside the intervals of both
	  conditioning variables will get a mixed color.}
	\item{plond.pale:}{ minimum and maximum paling, to be applied
	  for distance 0 and maximal distance from the interval.}
	\item{plcond.cex:}{ proportion of \code{cex} used to show
	  the points outside the interval}
      }
    }
  }
}
\note{\code{plpar} is a version of \code{ploptions} with different
  defaults. Its purpose is to make it easy to restore the
  'margin parameters' according to the last call to a high level
  pl function in order to complement a plot with low level ordinary
  R functions like \code{mtext} or \code{segments}, see Examples.

  The same effect can be achieved by setting the pl option
  \code{keeppar} to \code{TRUE}, either by calling \code{ploptions}
  or by setting \code{keeppar=TRUE} in the call to the high level
  pl function.
}
  
\author{Werner A. Stahel}
\note{In contrast to the \code{\link{options}} of R base,
  the \code{.plptions} list is stored in the \code{.GlobalEnv}
  (if any component is changed)
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
t.plopt <- ploptions(smooth.col="green", assign=FALSE)
t.plopt$smooth.col
ploptions("smooth.col") ## unchanged

## setting 'margin parameters' mar, oma, mgp, cex : 
ploptions(mar=c(2,2,1,1))  ## only active in pl funcions
  ## (unless  ploptions("keeppar")  is  FALSE )
ploptions(mar=c(2,2,1,1), setpar=TRUE)  ## changes the margins permanently
ploptions(setpar=TRUE) ## does the same. All the margin pars are
  ## set permanently to the values in  .ploptions

## switching 'margin parameters' between those used
## outside and inside high level pl functions
par(mar=c(2,2,5,2))
plyx(Sepal.Width~Sepal.Length, data=iris)
par("mar")
mtext("wrong place for text",3,1, col="red")
t.plo <- plpar()
par("mar")
mtext("here is the right place",3,1)

par(attr(t.plo, "oldpar"))  ## resets the 'margin parameters'
par("mar")
plyx(Sepal.Width~Sepal.Length, data=iris, keeppar=TRUE)
par("mar")
}
\keyword{utilities}
