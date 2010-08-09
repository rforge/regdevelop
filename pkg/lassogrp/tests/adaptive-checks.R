library(lassogrp)

data(asphalt)
d_1 <- asphalt[ -1,]

set.seed(1)## adaptive=TRUE -- needs crossvalidation for starting lambda --> randomness !!
## This used to fail ....  'x' must be an array of at least two dimensions
## because " x[,-interc.which ] "  needed an 'drop = FALSE'
lassoadp <- lasso(RUT ~ ., data = d_1, adaptive = TRUE)

stopifnot(lassoadp$converged,
	  lassoadp$lasso.terms == c(-1, 0,0,0, 1, 0,0, 1)
	  ,
	  all.equal(lassoadp$adaptcoef,
		    c(10.863689855, -8.72937971090), check.attr=FALSE)
	  )
