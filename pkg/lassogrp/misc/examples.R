source('/u/stahel/R/Pkgs/lassogrp/R/lassogrp.R')

d.asphalt <- read.table('/u/stahel/R/regdevelop/pkg/lassogrp/data/asphalt.dat',
                        header=T)
dd <- d.asphalt
rr <- lasso(log10(RUT)~log10(VISC)+ASPH+BASE+FINES+VOIDS+RUN, data=dd)
rr <- lasso(log10(RUT)~log10(VISC)+ASPH+BASE+FINES+VOIDS+RUN, nonpen=~1+RUN,
            data=dd)
rrf <- extract.lassogrp(rr, 13)

rra <- lasso(rr, adaptlambda=2.29)
rrr <- extract.lassogrp(rra, 11, data=dd, fitfun='lm')

d.blast <- read.table('/u/stahel/data/blast.dat',header=T)

dd <- d.blast
dd$xrand <- rnorm(nrow(dd))

rr <- lasso(log10(tremor)~charge+xrand+location, data=dd)
plot(rr)
rra <- lasso(rr, adaptlambda=-8)
plot(rra)

rrf <- extract.lassogrp(rr, 6)
rrr <- extract.lassogrp(rra, 3, data=dd, fitfun='regr')

dd1 <- dd
dd1[,c('distance','charge','xrand')] <-
  scale(dd1[,c('distance','charge','xrand')])


