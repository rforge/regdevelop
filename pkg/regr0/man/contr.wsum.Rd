\name{contr.wsum}
\alias{contr.wsum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  contrasts for zeroing the weighted sum of coefficients}
\description{Whereas \code{contr.sum} generates contrasts for factors in
  regression models that lead to a (unweighted) sum of estimated
  coefficients of zero, this function does so for a weighted sum.
  Typically, the weights are the frequencies of the levels of the factor
  in the dataset serving to fit the model.
}
\usage{
contr.wsum(n, y = NULL, w = NULL, contrasts = TRUE, sparse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a vector of levels for a factor, or the number of levels,
    or a factor, or a data.frame.
}
  \item{y}{response variable, used to drop observations with NA response
}
  \item{w}{weights for the levels
}
  \item{contrasts}{a logical indicating whether contrasts should be
    computed, see \code{?contr.sum}.
}
  \item{sparse}{logical indicating if the result should be sparse (of class
          'dgCMatrix'), using package 'Matrix'.
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
  (through the default of \code{userOptions}).
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

data(d.blast)
contr.wsum(d.blast$location)
wg <- table(d.blast$location)
( r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast) )
sum(wg * r.blast$allcoef$location[,1]) ## approximately 0
}
\keyword{ utilities }
\keyword{ regression }% __ONLY ONE__ keyword per line
