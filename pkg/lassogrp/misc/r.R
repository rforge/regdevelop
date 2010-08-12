options(myPkgs='/u/stahel/R/regdevelop/pkg')
Mlibrary(lassogrp)
# source('../R/lassogrp.R')

##- load('../data/asphalt.rda')
##- load('../data/splice.rda')
data(asphalt)
dd <- asphalt
dd$Random <- rnorm(nrow(dd))
rr <- lasso(log10(RUT)~log10(VISC)+ASPH+BASE+FINES+VOIDS+RUN,
  data=dd)
rr[c(1,10,15)]
extract.lassogrp(rr, lambda=2.5)

rra <- lasso(rr)
plot(rra)
rcv <- cv.lasso(rra)
