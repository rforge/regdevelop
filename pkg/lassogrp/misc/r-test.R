source('../R/lassogrp.R')
load('../data/splice.rda')

library(grplasso)
##- example(grplasso)
data(splice)
##- contr <- rep(list("contr.sum"), ncol(splice) - 1)
##  names(contr) <- names(splice)[-1]
##- fit.splice <-
##-   grplasso(y ~ ., data = splice, model = LogReg(), lambda = seq(65,0,-8),
##-            center = TRUE, standardize = TRUE)
##- fit.splice$fn.val
##- t.fit <-
##-   grplasso(y ~ ., data = splice, model = LogReg(), lambda = seq(24,0,-8),
##-            nonpen = ~ Pos.2)

r.splice <-
  lasso(y ~ ., data = splice, model = LogReg(), lambda = seq(65,0,-8))
r.splice$fn.val

t.splice <-
  lasso(y ~ ., data = splice, model = LogReg(), lambda = seq(24,0,-8),
           nonpen = ~ Pos.2)
t.splice

## =========================================================
d.blast <- read.table('/u/stahel/data/blast.dat', header=TRUE)
t.r <- lasso(log10(tremor)~location+log10(distance)+log10(charge),
                   data=d.blast)
t.r <- lasso(log10(tremor)~location+log10(distance)+log10(charge),
                   data=d.blast, subset=-(1:3))
t.rr <- lasso(log10(tremor)~location+log10(distance)+log10(charge),
              data=d.blast, subset=-(1:3), adaptive=TRUE,
              control=lassoControl(trace=1))
t.cv <- cv.lasso(t.r)
plot(t.r, type='crit', cv=t.cv, se=T)

t.pr <- predict(t.r)


t.rr <- extract.lassogrp(t.r,8)
t.rr
