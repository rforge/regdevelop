\name{plmatrix}
\alias{plmatrix}
\alias{panelDefault}
\alias{panelSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scatterplot Matrix}
\description{
  Plots a scatterplot matrix, for which the variables shown horizontally
  do not necessarily coincide with those shown vertically. If desired,
  the matrix is divided into several blocks such that it fills more than
  1 plot page.

  code{panelDefault} is the function to be applied for each subplot,
  unless a different function is provided. 
}
\usage{
plmatrix(x, y = NULL, data = NULL, panel = panelDefault, 
         nrow = 0, ncol = nrow, reduce = TRUE, keeppar = FALSE,
         xaxmar=NULL, yaxmar=NULL, xlabmar=NULL, ylabmar=NULL,
         oma=NULL, mar=rep(0.2,4), cex.diaglabels = 1.5, ...) 

panelDefault(xx, yy, indx, indy, plargs, ...)

panelSmooth(xx, yy, indx, indy, plargs, ...)
%% args xx and yy instead of x and y in order to distinguish them here
%%  in the .Rd file
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data for columns (x axis), or formula defining column
    variables. If it is a formula containing a left hand side,
    the left side variables will be used last.}
  \item{y}{data or formula for rows (y axis). Defaults to \code{x}}
  \item{data}{data.frame containing the variables in case \code{x}
    or \code{y} is a formula}
  \item{panel}{a function that generates the marks of the individual
    panels, see Details.
    Defaults essentially to \code{\link{points}} or \code{\link{text}}
    depending on the argument \code{pch},
    to \code{\link{plmboxes}} if  'x'  is a factor and  'y'  is not,
    or to a modification of \code{sunflowers} if both are factors.
  }
%%-   \item{mbox}{logical: Should mixed plots of factors and continuous
%%-     variables be displayed as "mbox plots"? See \code{?\link{plmboxes}}.
%%-   }
  \item{nrow}{number of rows of panels on a page}
  \item{ncol}{number of columns of panels on a page}
  \item{reduce}{if y is not provided and \code{reduce==TRUE},
    the first row and the last column are suppressed.}
  \item{keeppar}{logical. If TRUE, the graphical settings
    "mfrow", "mar", "cex", and "mgp" will be maintained, otherwise,
    the old values will be restored (default).}
  \item{xaxmar, yaxmar}{margin in which the axis (tick marks and
    corresponding labels) should be shown: either 1 or 3 for
    \code{xaxmar} and 2 or 4 for \code{yaxmar}.
  }
%%   \item{markextremes}{proportion of extreme residuals to be labeled.
%%     If all points should be labeled, let \code{markextremes=1}.}
%%   \item{xaxmar, yaxmar}{in which margin should the x- [y-] axis be
%%     tick-marked?}
  \item{xlabmar, ylabmar}{in which margin should the x- [y-] axis be
    labelled?}
%%  \item{mar, mgp}{margin width and margin parameters, see
%%    \code{\link{par}}.
%%    Defaults are taken from \code{UserOptions}.
    %% mar is needed explicitly because otherwise
    %% it is confounded with  markextremes
%%  }
%%   \item{vnames}{labels for the variables}
%%   \item{main}{main title for the plot (to be repeated on each plot
%%     page)}
%%   \item{xlabs, reslabs}{vectors of labels for x- and y-axes, a string for
%%     each variable. This will override the \code{vnames} argument for the
%%     axis labels, but not for the diagonal panels.}
  \item{oma}{width of outer margins, ee \code{\link{par}}}
  \item{mar}{margins for the panels}
  \item{cex.diaglabels}{Character expansion for labels appearing in the
    "diagonal" of the scatterplot matrix (if present)}
  \item{\dots}{further arguments passed to the \code{panel} function
    and possibly further to functions called by the panel function\cr}

  Arguments for \code{panelDefault} and \code{panel.smooth}
  \item{xx}{values of the horizontal variable}
  \item{yy}{values of the vertical variable}
  \item{indx}{index of the variable shown horizontally, among the
    \code{y} variables}
  \item{indy}{index of the variable shown horizontally, among the
    \code{y} variables}
  \item{plargs}{result of calling \code{pl.control}.
    If \code{NULL}, \code{pl.control} will be called to generate it.
    Its components are often needed to generate the panel.
  }\cr
}
\details{
  The \code{panel} function can be user written. It needs \eqn{>=6}
  arguments which must correspond to the arguments of
  \code{panelDefault}, see above. If some arguments are not used,
  just introduce them as arguments to the function anyway in order to
  avoid (unnecessary) error messages and stops.
  The names \code{pch} and \code{col} must appear among the argument
  names.

  \code{panel.smooth} only calls \code{\link{panel.smooth}} of the
  \code{graphics} package, 
  absorbing the arguments that are needed but not used.

  Since large scatterplot matrices lead to tiny panels, \code{plmatrix}
  splits the matrix into blocks of at most \code{nrow} rows and
  \code{ncol} columns. If these numbers are missing, they default to
  \code{nrow=5} and \code{ncol=6} for landscape pages, and to
  \code{nrow=8} and \code{ncol=5} for portrait pages.
}
\value{
  none
}
\author{Werner A. Stahel, ETH Zurich}
\seealso{\code{\link{pairs}}, \code{\link{plyx}}
}
\note{There are many more arguments, obtained from \code{pl.control},
  see \code{?\link{pl.control}}. These can be passed to \code{plmatrix}
  by an argument \code{plargs} that is hidden in the ... argument list.
}
\examples{
plmatrix(iris, pch=as.numeric(Species))
plmatrix(~Sepal.Length+Sepal.Width, ~Petal.Length+Petal.Width,
    data=iris, smooth=TRUE, plab=substr(Species,1,2))
}
\keyword{hplot}
