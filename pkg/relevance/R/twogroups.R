twogroups <- #
  function(x, y, var.equal=TRUE, testlevel=0.05, RelTh=0.25)
{ ## effect (group difference) and relevance 
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  lpq <- 1-testlevel/2
  lrlvth <- RelTh
  lnx <- length(x)
  lny <- length(y)
  ln <- lnx+lny
  lmnx <- mean(x)
  lmny <- mean(y)
  lssx <- sum((x-lmnx)^2)
  lssy <- sum((y-lmny)^2)
  lv <-
    if (var.equal) ln*(1/lnx+1/lny)*(lssx+lssy)/(ln-2)
    else  ln*(lssx/(lnx-1)/lnx+lssy/(lny-1)/lny)
  lq <- qt(lpq, ln-2)
  lciwid <- lq*sqrt(lv/ln)
  leffci <- lmny-lmnx + c(0, -1, 1)*lciwid
  lrlvci <- leffci/sqrt(lv)/lrlvth
  list(estimate=leffci, rlv=lrlvci, V=lv,
       sig0=leffci[1]/lciwid, sigth=(lrlvci[1]-1)*2/diff(lrlvci[2:3]))
}
