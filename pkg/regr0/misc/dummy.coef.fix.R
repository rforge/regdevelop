dummy.coef.fix <- function (object, use.na = FALSE, ...) 
{
    xl <- object$xlevels
    if (!length(xl)) 
        return(as.list(coef(object)))
    Terms <- terms(object)
    tl <- attr(Terms, "term.labels")
    int <- attr(Terms, "intercept")
    facs <- attr(Terms, "factors")[-1, , drop = FALSE]
    Terms <- delete.response(Terms)
    mf <- object$model
    if (is.null(mf)) mf <- model.frame(object)
    xtnm <- dimnames(facs)[[1]]  ## names
    xtlv <- lapply(mf[,xtnm, drop=FALSE],function(x) levels(x)) ## levels
    xtnl <- pmax(sapply(xtlv,length),1)  ## number of levels
    lterms <- apply(facs, 2L, function(x) prod(xtnl[x > 0]))
    nl <- sum(lterms)
    ## df.dummy: data frame of vars
    args <- setNames(vector("list", length(xtnm)), xtnm)
    for (i in xtnm)
        args[[i]] <- if (xtnl[[i]] == 1)  rep.int(1, nl)    else
          factor(rep.int(xtlv[[i]][1L], nl), levels = xtlv[[i]])
    df.dummy <- as.data.frame(args) # do.call("data.frame", args)
    names(df.dummy) <- xtnm
    ## rnn: names of rows
    pos <- 0
    rn <- rep.int(tl, lterms)
    rnn <- rep.int("", nl)
    for (j in tl) {
        i <- unlist(xtnm[facs[, j] > 0])
        ifac <- i[xtnl[i] > 1]
        if (length(ifac) == 0L) {
            rnn[pos + 1] <- j
        }
        else if (length(ifac) == 1L) {
            df.dummy[pos + 1L:lterms[j], ifac] <- xtlv[[ifac]]
            rnn[pos + 1L:lterms[j]] <- as.character(xtlv[[ifac]])
        }
        else {
            tmp <- expand.grid(xtlv[ifac])
            df.dummy[pos + 1L:lterms[j], ifac] <- tmp
            rnn[pos + 1L:lterms[j]] <- apply(as.matrix(tmp), 
                1L, function(x) paste(x, collapse = ":"))
        }
        pos <- pos + lterms[j]
    }
    attr(df.dummy,"terms") <- attr(mf,"terms")
    lcontr <- object$contrasts
    lci <- sapply(df.dummy,is.factor)
    lcontr <- lcontr[names(lci)[lci]] ## factors with 1 level have disappeared (?) 
    mm <- model.matrix(Terms, df.dummy, lcontr, xl)
    if (any(is.na(mm))) {
        warning("some terms will have NAs due to the limits of the method")
        mm[is.na(mm)] <- NA
    }
    coef <- object$coefficients
    if (!use.na) 
        coef[is.na(coef)] <- 0
    asgn <- attr(mm, "assign")
    res <- setNames(vector("list", length(tl)), tl)
    for (j in seq_along(tl)) {
      keep <- which(asgn == j)
      cf <- coef[keep]
      na <- is.na(cf)
      ij <- rn == tl[j]
      if (any(na)) {
        rsj <- setNames(drop(mm[ij, keep[!na], drop = FALSE] %*% 
                               cf[!na]), rnn[ij])
        rsj[apply(mm[ij, keep[na], drop=FALSE], 1, function(x) any(x!=0))] <- NA
      } else
        rsj <- setNames(drop(mm[ij, keep, drop = FALSE] %*% cf), rnn[ij])
      res[[j]] <- rsj
    }
    if (int > 0) {
        res <- c(list(`(Intercept)` = coef[int]), res)
      }
    BR()
    class(res) <- "dummy_coef"
    res
}
