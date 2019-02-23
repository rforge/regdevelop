\name{contr.wsum}
\alias{contr.wsumpoly}
\alias{contr.wsum}
\alias{contr.wpoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  contrasts for zeroing the weighted sum of coefficients or
  orthogonalizing polynomials}
\description{Whereas \code{contr.sum} generates contrasts for factors in
  regression models that lead to a (unweighted) sum of estimated
  coefficients of zero, this function does so for a weighted sum.
  Typically, the weights are the frequencies of the levels of the factor
  in the dataset serving to fit the model.

  \code{contr.wpoly} adjusts the generated variables for a polynomial
  to render them orthogonal for the given frequencies of an ordered
  factor instead of doing this only for equal weights.

  \code{contr.wsumpoly} will do either of these, depending on the
  arguments. 
}
\usage{
contr.wsumpoly(n, scores = NULL, y = NULL, w = NULL,
               contrasts = TRUE, sparse = FALSE, poly = NA)
contr.wsum(n, scores = NULL, y = NULL, w = NULL, contrasts = TRUE,
           sparse = FALSE)
contr.wpoly(n, scores = NULL, y = NULL, w = NULL, contrasts = TRUE,
            sparse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a factor, a vector of levels for a factor, the number of levels,
    or a data.frame.
  }
  \item{scores}{
    scores for the levels of the factor, see
    \code{?\link{contr.poly}}. Ignored in \code{contr.wsum}.
  }
  \item{y}{response variable, used to drop observations with NA response
  }
  \item{w}{weights for the levels. If \code{n} if this is a factor,
    the frequency vector is used and \code{w} is ignored.
  }
  \item{contrasts}{a logical indicating whether contrasts should be
    computed, see \code{?\link{contr.sum}}.
  }
  \item{sparse}{logical indicating if the result should be sparse (of class
          'dgCMatrix'), using package 'Matrix'.
  }
  \item{poly}{(only for \code{contr.wsumpoly})
    If TRUE, polynomial contrasts are generated;
    if FALSE, (weighted) sum contrasts are given;
    if NA (default), sum or polynomial contrasts will result for 
    unordered and ordered factors, respectively.
  }
}
\value{
  If \code{n} is a number, a vector of levels, or a factor,
  the value is the contrast matrix that can be assigned to
  \code{attr(factor, "contrasts")} and then will be used in setting up
  the model.matrix for a regression or anova model.
  If it is a data.frame, this attribute is attached to each unordered
  factor and the complemented data.frame is returned.
}
\author{Werner A. Stahel, Seminar f. Statistik, ETHZ
}
\note{Using these contrasts makes individual coefficients interpretable
  as the effect of the respective level compared to the weighted average
  of effects.
  It also makes the intercept -- for which the interpretation depends on
  the choice of contrasts -- as precisely estimable as possible.

  The function \code{regr} chooses these contrasts as a default
  (through the default of \code{useroptions}).
  This, however, only works for terms that are factors when before the
  formula of the model is translated into the model.frame.
  (This will probably change soon.)  
}
\seealso{
\code{link{contrasts}}, \code{link{contr.sum}}
}
\examples{
contr.wsum(factor(c(1,1,1,2,3,3)))
##    [,1] [,2]
##  1  1.0  0.0
##  2  0.0  1.0
##  3 -1.5 -0.5

contr.wpoly(factor(c(1,1,1,2,3,3)))
##       ^0     ^1     ^2
##  [1,]  1 -0.833  0.138
##  [2,]  1  0.167 -0.828
##  [3,]  1  1.167  0.207


data(d.blast, package="plgraphics")
contr.wsum(d.blast$location)
wg <- table(d.blast$location)
( r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast) )
sum(wg * r.blast$allcoef$location[,1]) ## approximately 0
}
\keyword{ utilities }
\keyword{ regression }% __ONLY ONE__ keyword per line
