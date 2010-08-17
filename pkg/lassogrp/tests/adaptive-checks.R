library(lassogrp)

data(asphalt)
d_1 <- asphalt[ -1,]

set.seed(1)## adaptive=TRUE -- needs crossvalidation for starting lambda --> randomness !!
## This used to fail ....  'x' must be an array of at least two dimensions
## because " x[,-interc.which ] "  needed an 'drop = FALSE'
lassoadp <- lasso(RUT ~ ., data = d_1, adaptive = TRUE)
lassoadp

with(lassoadp,
 stopifnot(converged
	   ## 1.0-1:
	   ## lasso.terms == c(-1, 0,0,0, 1, 0,0, 1)
	   , names(lasso.terms) == names(terms)
	   , unname(lasso.terms) == c(-1, 0, 1, 1, 1, 1, 0)
	   , all.equal(unname(adaptcoef),
		       c(-23.8785241, 3.78992492, 2.91446966,
			 -9.85831061, 0.0388053332))
	   ## 1.0-1: c(10.863689855, -8.72937971090), check.attr=FALSE)
	   ))
