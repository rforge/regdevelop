source('../R/lassogrp.R')

load('../data/asphalt.rda')
load('../data/splice.rda')
rr <- lasso(log10(RUT)~log10(VISC)+ASPH+BASE+FINES+VOIDS+RUN,
  data=asphalt)
rr[c(1,10,15)]
extract.lassogrp(rr, lambda=2.5)
