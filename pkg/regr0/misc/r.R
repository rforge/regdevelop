source('../R/regr.R')

r.blast <-
  regr(log10(tremor)~location+log10(distance)+log10(charge), data=d.blast) 
