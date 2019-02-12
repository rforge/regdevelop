## ==================================================================
## Plotting functions, low and high-level
## -------------------------------------------------------------------------
pl.control <-
  function(x=NULL, y=NULL, data = NULL, subset = NULL, transformed = TRUE,
           gensequence = NULL,
           psize = NULL, plab = FALSE, pch = NULL, pcol = NULL,
           markextremes = NULL, smooth = NULL,
           cex = NULL, xlab = NULL, ylab = NULL, varlabels = NULL,
           ycol = NULL, ylty = NULL, ypch = NULL, 
           main = NULL, sub = ":", .subdefault = NULL, mar = NULL,
           ## needed because it hides  markextremes  otherwise
           ploptions = NULL, .environment. = parent.frame(), assign = TRUE, ... )
  ## get data for plotting, collect and check arguments
  ## do preparations that are common to all plots
  ## --------------------------------------------------------------
{
  ## ---
  lcall <- lcl <- match.call()
  ## ploptions
  ploptions <- i.def(ploptions, get(".ploptions", pos=1))
  lnmd <- setdiff(names(ploptionsDefault), names(ploptions) )
  if (length(lnmd)) ploptions <- c(ploptions, ploptionsDefault[lnmd])
  largsplo <- setdiff(names(lcl)[-1],
                      c("xlab","ylab", ## "cex","markextremes",
                       i.argPlcontr, i.argPldata))
  if (length(largsplo)) {
    lcl <- lcl[c("",largsplo)]
    lcl[1] <- list(quote(list))
    lcl <- as.call(lcl)
    lls <- eval(lcl, envir=parent.frame())
    ploptions <- ploptions(list=lls, ploptions=ploptions, assign=FALSE)
  }
  lcl <- lcall
  ## --- data
  ldtl <- NULL
  if (length(data)) {
    ldtl <- if (length(ltit <- tit(data))) ltit
    if (is.null(ldtl) & is.name(substitute(data))) ldtl <- substitute(data)
  }
  lftext <- lform <- NULL
  lformarg <- lynames <- lxnames <- NULL
  if (length(i.def(x, NULL, TRUE, NULL))) {
    if (is.atomic(x)&&is.character(x))  ## names of variables
      lxnames <- x
    else {  ## matrix or data.frame
      if (is.matrix(x)|is.atomic(x)) {
        x <- as.data.frame(x)
      }
      if (is.data.frame(x)) {
        lxnames <- names(x)
        if (is.null(ldtl))
          ldtl <- if (is.null(ltit <- tit(x)))
                    as.character(attr(x, "dname")) else ltit
        data <- if (is.null(data)) x else cbind(x, data)
      } else {  ## formula
        if (is.formula(x)) lform <- x
        else stop("!pl.control! Unsuitable argument 'x'")
      }
    }
    if (is.null(lxnames)) {
      lvars <- getvarnames(lform, transformed=transformed)
      lxnames <- lvars$xvar
      lynames <- lvars$yvar
    }
    if (is.null(lform))
      lform <-
        as.formula(paste(paste(lynames, collapse="+"), "~",
                         paste(lxnames, collapse="+")))
    environment(lform) <- environment() ## !!!???
    lftext <- format(lform)
  }
  ## y
  lyform <- NULL
  if (length(y)) {
    if (is.atomic(y)&&is.character(y)) 
      lynames <- y
    else {
      if (is.matrix(y)|is.atomic(y)) {
        y <- as.data.frame(y)
      }
      if (is.data.frame(y)) {
        lynames <- names(y)
        data <- if (is.null(data)) y else cbind(y, data)
      }
      else {
        if (is.formula(y)) lyform <- y
        else stop("!pl.control! Unsuitable argument 'y'")
      }
    }
    if (is.null(lynames))
      lynames <- getvarnames(lyform, transformed=transformed)$varnames
    if (is.null(lyform))
      lyform <- as.formula(paste("~", paste(lynames, collapse="+")))
    environment(lyform) <- environment() ## ?
    lftext <- paste(shorten(substring(format(lyform),2,100),50), lftext)
    ## fixes maximal length of formula text
  }
  ##
  lformarg <- list(lxnames,lynames)
  lvarnames <- c(lxnames, lynames)
  ## --- data
  ## ltransformed <- i.def(lcl$transformed, TRUE)
  if (!is.null(attr(data,"terms"))) { ## data is model.frame
    if (!transformed) {
      warning(":plotregr.control! Raw data not available.",
              "I can only use transformed data.")
      transformed <- TRUE
    }
  }
  ## get variables from  lform  and  data
  largs <- c(i.argPldata, "transformed")
  ## --- variables
  if (length(varlabels)) 
    if (length(names(varlabels))==0) {
      warning(":pl.control: 'varlabels' must have names")
      lcl$varlabels <- varlabels <- NULL
    }
  if (length(lvarnames)||length(varlabels)||
      any(names(lcl)%in%i.argPldata)) {
    lcl <- c(list(quote(getvariables), formula=lformarg, data=data),
               as.list(lcl[intersect(largs, names(lcl))]),
               envir=.environment.)
    lcl <- as.call(lcl)
    ##        ----
    lpldata <- eval(lcl, envir=environment(lform))
    ##        ----
    if (inherits(lpldata, "pl-error"))
      stop("!pl.control! ", attr(lpldata, "message"))
    ## set varlabels in a special case
    if (length(lxnames)==1 && is.formula(lcx <- lcall$x) && length(lcx)==2)
        attr(lpldata[[lxnames]], "varlabel") <- sub("~","",format(lcx))
    if (length(lynames)==1 && is.formula(lcy <- lcall$y) && length(lcy)==2)
        attr(lpldata[[lynames]], "varlabel") <- sub("~","",format(lcy))
    ##
    if (length(lgroup <- lpldata[["(group)"]]))
      if (length(llb <- as.character(lcall$group))<=20)
        attr(lpldata[["(group)"]], "varname") <- llb
    ## --- subset
    if (length(lcl$subset)) {
      lsub <- eval(lcl$subset, data)
      lpldata <- plsubset(lpldata, lsub)
    }
    ## --- attributes of variables
    ## plrange
    if(any(c("plrange","xlim","ylim")%in%names(lcall))) 
      lpldata <- i.setplrange(lpldata, ...)
    lvarnames <- attr(lpldata,"variables")
    if (length(lvarnames)) {
      lpldata[,lvarnames] <- 
        genvarattributes(lpldata[,lvarnames, drop=FALSE], ynames = lynames,
                         ycol = ycol, ylty = ylty, ypch = ypch,
                         varlabels = varlabels, ploptions=ploptions)
      lnr <- nrow(lpldata)
      lnobs <- lnr-median(sumna(lpldata[,lvarnames]))
    } else lnr <- 0
  } else {
    lpldata <- NULL
    lvarnames <- NULL
    lnr <- 0
  }
  if (lnr==0) {
    if (length(data)) {
      lnr <- NROW(data)
      lnobs <- lnr-median(sumna(data))
    } else {
      warning(":pl.control: data not found. I set 'nobs' to 100")
      lnobs <- 100
    }
  }
  ## single column
  ##  if (length(lpldata) && NCOL(lpldata)==1 && u.notfalse(x))
  ## { ## NCOL(NULL) is 1
  ##!!!
  if (length(lpldata) && u.notfalse(gensequence) &&
      (is.null(lxnames)|is.null(lynames))) { 
    lpldata <-
      cbind("(sequence)"= structure(1:NROW(lpldata), varlabel="sequence"),
            lpldata)
    attr(lpldata,"yvar") <- lvarnames
    attr(lpldata,"xvar") <- "(sequence)"
    lvarnames <- c("(sequence)", lvarnames)
  }
  ## -------------------------
  lvnm <-
    setdiff(intersect(names(data), c("(pch)","(plab)","(pcol)","(psize)")),
                  names(lpldata))
  if (length(lvnm))
    lpldata[,lvnm] <- data[,lvnm, drop=FALSE]
##-     lpldata <- transferAttributes(cbind(lpldata, data[,lvnm, drop=FALSE]),
##-                                   lpldata)
  ## labels anmd plotting character
  ## priorities:  plab , pch , row.names
  lpch <- lpldata$"(pch)"
  ## default plotting character
  if (length(lpch)) {
    if (is.factor(lpch)) lpch <- as.numeric(lpch)
    if (is.character(lpch) && any(nchar(lpch)>1)) {
      warning(":pl.control: 'pch' must be an integer or a single character.",
              " Use 'plab' to label points with strings")
      lpch <- substr(lpch,1,1)
    }
    if (is.logical(lpch)) lpch <- as.numeric(lpch)
    if (is.numeric(lpch)&&any(lpch<0|lpch>30)) { ## !!! 30? -> 
      warning(":plContol: unsuitable plotting character ('pch')")
      lpch <- NULL
    }
    lpldata$"(pch)" <- lpch
  }
  lplab <- lpldata$"(plab)"
  ## --- row.names
  lrown <- substring(row.names(data),1,3)
  if (length(lrown)==0) lrown <- as.character(1:lnr)
  lplabel <- lrown
  lIplab <- length(lplab)>0 && is.logical(lplab) && all(lplab) ## original arg T
  if (lIplab)  lpldata$"(plab)" <- lplab <- lplabel
  else
    if (length(lplab)) lplabel <- as.character(lplab)
  ## factors -> legend!!!
  ##     if (is.factor(lplab)) as.numeric(lplab) else lplab 
  ## now, lplabel always useful
  ## ----------------------------------------------------
  ## more ploptions
  ##
  if (length(ploptions$smoothline.col)==1)
    ploptions$smoothline.col[2] <-
      colorpale(ploptions$smoothline.col, i.getploption("smoothline.pale"))
  ## 
##-   ploptions$cex <- i.getplopt(cex)
##-   ploptions$markextremes <- i.getplopt(markextremes)
  ## ---
  lmardf <- i.getplopt(mar)
  if (length(lynames)>1) lmardf[4] <- lmardf[2] ## need space at the right
  ploptions$mar <- i.def(mar, lmardf)
  ## --- condprobRange
  ploptions$condprobRange <-
    if (length(ploptions$condprobRange)==0) {
      if (lnobs>50) c(0,0) else c(0.05,0.8) }
    else c(ploptions$condprobRange,1)[1:2]
  ## --- smooth
  lsmgrp <- lpldata$"(smooth.group)"
  lsmgrplab <- levels(lsmgrp)
  ## smooth
  lnsm <- lnobs
  lnsmgrp <- length(unique(lsmgrp))
  if (length(lsmgrp)) lnsm <- lnobs/lnsmgrp
  ploptions$smooth <- i.getplopt(smooth) ## i.def(ploptions$smooth, 0, 2, 0)
  ploptions$smooth.par <-
    i.def( ploptions$smooth.par, 5*lnsm^log10(1/2)*(1+inherits(x,"glm")) )
  ## --- main
  main <- i.def(main, "", "", "")
  sub <- i.def(sub, NULL, ":", NULL)
  ## --- more arguments
  ## refline <- i.def(refline, TRUE)
  ## ------------------------------------------------------------
  ## result of pl.control
  rr <- list(
    pldata = lpldata, formula = lform,
    ##  xvar = attr(lpldata,"xvar"), yvar = attr(lpldata,"yvar"),
    nobs = lnobs, transformed = transformed,
    pch = lpch, plabel = lplabel, plab = lIplab, ##plabna = lplabna, ???
    smooth.ngroups = lnsmgrp, smooth.grouplab = lsmgrplab,
    dataLabel = ldtl, main = main, sub = sub, .subdefault = .subdefault,
    ploptions = ploptions, datetime = date()
    )
  if (u.notfalse(assign)) assign(".plargs", rr, pos=1)
  rr
} ## end of  pl.control

## ===================================================================
getvarnames <-
  function(formula, transformed=FALSE)
{ ## get  varnames 
  if (is.character(formula))
    return(list(varnames=formula, xvar=formula, yvar=NULL))
  if (is.null(formula)) return(list(varnames=NULL, xvar=NULL, yvar=NULL))
  ##    formula <- as.formula(paste("~",paste(formula,collapse="+")))
  if (is.list(formula)) formula <- formula(formula)
  if (!is.formula(formula)) stop("!getvarnames! invalid argument 'formula'")
  lyv <- NULL
  lxv <- lvnm <-
    if (transformed)
      rownames(attr(terms(formula[c(1:2)]), "factors"))
    ## attr(..."variables") is language
    else all.vars(formula[1:2])
  if (length(formula)==3) {
    lyv <- lxv
    lxv <- if (transformed) rownames(attr(terms(formula[-2]), "factors"))
           else all.vars(formula[-2])
    lvnm <- c(lxv, lvnm)
  }
  list(varnames=lvnm, xvar=lxv, yvar=lyv)
}
## =====================================================================
getvariables <-
  function (formula, data = NULL, transformed = TRUE,
            envir = parent.frame(), ...)
{
  ## similar to get_all_vars , different error handling; generate is.fac
  if (is.list(formula)) {
    lvnm <- getvarnames(formula[[1]], transformed=transformed)
    lvnmy <- if (length(lfoy <- formula[[2]]))
               getvarnames(lfoy, transformed=transformed)  else NULL
    lvarnames <- unique(c(lvnm$varnames, lvnmy$varnames))
##    lenv <- environment(formula[[1]])
  } else {
    lvnm <- if (length(formula)) getvarnames(formula, transformed=transformed)
    lvnmy <- NULL
    lvarnames <- lvnm$varnames
##    lenv <- environment(formula)
  }
  ## data
  if (length(data)==0) 
    data <- environment(formula)
  else if (!is.data.frame(data) && !is.environment(data) && 
           !is.null(attr(data, "class"))) 
    data <- as.data.frame(data)
  else if (is.array(data)) 
    stop("!getvariables! 'data' must be a data.frame, not a matrix or an array")
  ##
  lenv <- envir ## parent.frame() # environment(formula)
  rownames <- .row_names_info(data, 0L)
  rr <- NULL
  lvn <- lvarnames
  if (transformed) {
    liv <- lvarnames%in%colnames(data)
    rr <- data[lvarnames[liv]]
    lvn <- lvarnames[!liv]
  }
  ##-   if (length(grep("cbind",yvar)))
  ##-     yvar <- all.vars(as.formula(paste("~", varnames[1]))) 
  if (length(lvn)) {
    inp <- parse(text = paste("list(", paste(lvn, collapse = ","),")"),
                 keep.source = FALSE)
    variables <- eval(inp, data, lenv)
    names(variables) <- lvn
    ## drop attributes for transformed variables
    ljtr <- lvn %in% names(data)
    if (any(!ljtr))
      for (lj in lvn[!ljtr]) {
        lia <- names(attributes(variables[[lj]])) %in%
          c("ticksat", "ticklabelsat", "ticklabels")
        attributes(variables[[lj]])[lia] <- NULL
      }
    lvmode <- sapply(variables, mode)
    if (any(li <- lvmode%nin%c("numeric","character","logical")))
      stop("!getvariables! variable(s)  ",paste(lvn[li],collapse=", "),
           "  has(have) wrong mode")  ## !!! convert into 'delayed error' as below
    variables <- data.frame(variables, check.names=FALSE)
    if (length(rr)) {
      if (nrow(rr)!=nrow(variables))
        stop("!getvariables! differing numbers of rows: ",
             nrow(rr), ", ", nrow(variables))
      rr <- cbind(rr, variables)
    } else rr <- variables
  }
  ## rownames
  if (is.null(rownames) && (length(resp <- attr(formula, "response")) > 0) ) {
    lhs <- rr[[resp]]
    rownames <- if (is.matrix(lhs)) rownames(lhs) else names(lhs)
  }
  ## --- extras
  extras <- substitute(list(...))
  extranames <- names(extras[-1L])
  if (length(extranames)) {
    extras <- eval(extras, data, lenv)
    names(extras) <-
      if (length(extranames)) paste("(",extranames,")",sep="") else NULL
    lexl <- sapply(extras, length)
    if (any(lexwrong <- lexl%nin%c(1,nrow(rr)))) {
      warning(":getvariables: differing numbers of rows:   ",
              nrow(rr), ", ", paste(lexl, collapse=", "),
              "\n  I drop   ", paste(extranames[lexwrong], collapse=", "))
      extras <- extras[!lexwrong]
    }
    ## --- bind together
    if (length(extras))
      rr <- cbind(rr,data.frame(extras, check.names=FALSE,
                                stringsAsFactors=FALSE))
  }
  len <- sapply(rr,NROW)
  rr <- rr[len>0]
  messg <- NULL
  len <- len[len>0]
  nobs <- if(is.data.frame(data)) nrow(data) else length(rr[[1]])
  if (any(li <- len%nin%c(1,nobs))) {
    messg <- paste(ifelse(sum(li)==1, "Variable  ","Variables  "),
                   paste(c(lvarnames, extranames)[li], collapse=", "),
                   ifelse(sum(li)==1, "  has inadequate length",
                          "  have inadequate lengths"),sep="")
    fatal <- any(li[1:length(rr)])
    if (!fatal) {
      warning(messg <- paste(messg,"\n   and will be ignored"))
      rr <- rr[!li]
      rr <- setNames(as.data.frame(rr), names(rr))
      if (!is.null(rownames)) 
        attr(rr, "row.names") <- rownames
    } else {
      class(rr) <- "pl-error" 
      attr(rr,"message") <- messg
      warning("!pl.control/getvariables! Error: ", messg) ## why not stop? 
      return(rr)
    }
  } 
  rr <- setNames(as.data.frame(rr), names(rr)) ## avoid modifying names
  attr(rr,"variables") <- lvarnames
  if (!is.null(rownames)) attr(rr, "row.names") <- rownames
  if (is.null(lvnmy)) {
    attr(rr, "xvar") <- lvnm$xvar
    attr(rr, "yvar") <- lvnm$yvar
  } else {
    attr(rr, "xvar") <- lvnm$varnames
    attr(rr, "yvar") <- lvnmy$varnames
  }
  if (!is.null(messg)) { ## warning
    attr(rr, "message") <- messg
    class(rr) <- c("pl-warning", class(rr))
  }
  rr
}
## ------------------------------------------------------------------------
i.getvarattribute <- function(attr, value=NULL, data, ploptionscomp, drop=0)
{ ## get plot property  attr  either from a direct argument  value ,
  ## or a stored property  ploptionscomp
  ## The function avoid a duplicated use of the property
  var <- names(data)
  pla <- setNames(rep(NA, length(var)), var)
  ## set values from argument  value
  if(length(value)) {
    lnm <- intersect(names(value), var)
    if(length(lnm)) {
      pla[lnm] <- value[lnm]
      var <- setdiff(var, lnm)
    } else {
      if(length(value)%in%c(1,length(pla))) {
        pla[] <- value
        return(pla)
      } else 
        warning(":getPlattribute: unsuitable length of (unnamed) attributes\n",
                paste(value, collapse=", "))
    }
  }
  ## set values from stored attribute in  data
  lpr <- sapply(data[,var], function(x) attr(x,attr))
  lpr <- unlist(lpr[sapply(lpr, length)>0]) ## drop empty 
  lnm <- intersect(names(lpr), var)
  if(length(lnm)) {
    pla[lnm] <- lpr[lnm]
    var <- setdiff(var, lnm)
  }
  ## get unset elements from default
  if(length(var)) {
    ldef <- ploptionscomp
    ## drop first element(s), like "black" for color
    if (drop>0) ldef <- ldef[-seq_len(drop)]
    ldef <- rep(setdiff(ldef, pla), length=length(var))
    pla[var] <- ldef
  }
  pla
}
## ------------------------------------------------------------------------
i.setvarattribute <- function(attr, value=NULL, data)
{ ## data <- i.setvarattributes("plrange", c("T","ra"), dd, list(T=c(15,20)))
  var <- names(data)
  if (!is.list(value))
    value <- setNames(rep(list(value), length(var)), var)
  lnm <- names(value)
  ljvar <- match(lnm, var, nomatch=0)
  if (any(ljvar==0)) {
    warning(":i.setvarattribute: name ",paste(lnm[ljvar==0], collapse=", "),
            "of 'value' not in names(data). Attribute '",attr,
            "' not set.")
  }
  lvar <- var[ljvar]
  if (length(lvar))
    for (lv in lvar)
      attr(data[,lv], attr) <- value[[lv]]
  data
}
## ----------------------------------------------------------------
i.setplrange <- function(data, plrange=NULL, xlim=NULL, ylim=NULL, ...)
{
  lf.checklim <- function(lim, varnames) {
    if (!is.list(lim)) lim <- list(lim)
    lxv <- attr(data, varnames)
    if (length(lnm <- names(lim))) {
      lim <- lim[intersect(lnm, lxv)]
    } else {
      if (length(lim)==length(lxv))
        names(lim) <- lxv
      else {
        warning(":pl.control/setplrange: 'xlim' or 'ylim' not suitable")
        lim <- NULL
      }
    }
    lim
  }
  ## -----
  if (!is.null(xlim)) xlim <- lf.checklim(xlim, "xvar")
  if (!is.null(ylim)) ylim <- lf.checklim(ylim, "yvar")
  lim <- c(xlim,ylim)
  if (length(plrange)) {
    if (!is.list(plrange)) plrange <- list(plrange)
    lim <- c(lim,plrange)
  }
  i.setvarattribute("plrange", lim, data)
}
## ===============================================================
plinnerrange <-
  function(innerrange, data, factor = 4.0, FUNC=robrange)
{ ## determine inner plot range
  ## if innerrange is a list or a matrix, leave it alone
  lIcq <- inherits(data, "condquant")
  ldt <- cbind( if (lIcq) c(data[,1:3]) else data )
  innerrange <- i.def(innerrange, TRUE, TRUE, FALSE)
  if (length(innerrange)>1)
    if (any(dim(cbind(innerrange))!=c(2,NCOL(ldt)))) {
      warning(":plot.regr/innerrange: unsuitable argument  innerrange ")
      innerrange <- TRUE
    }
  if (length(innerrange)==1 && is.logical(innerrange))
    innerrange <-
      if (innerrange) apply(ldt, 2, FUNC, fac=factor[1]) else
      matrix(FALSE, 2, NCOL(ldt))
  if ((!is.list(innerrange))&length(innerrange)==2)
    innerrange <- as.matrix(innerrange)
  innerrange
}
## ---------------------------------------------------------------
genvarattributes <-
  function(data, ynames = NULL, ycol = NULL, ylty = NULL, ypch = NULL,
           varlabels = NULL, innerrange.limits = NULL,
           ploptions = NULL, replace = FALSE)
{
  if (!is.data.frame(data))
    stop("!genvarattributes! 'data' must be a data frame")
  ## ---
  tickintervals <- i.getploption("tickintervals")
  ljt <- i.getploption("jitter")
  lljt <- length(ljt)
  if (lljt)
    ljt <- if (is.list(ljt)) ljt
            else setNames(rep(ljt, ncol(data)),
                          colnames(data))
  jitter.factor <- i.getploption("jitter.factor")
  ## innerrange
  lir <- i.def(i.getploption("innerrange"), FALSE, TRUE)
  lirg <- setNames(rep(list(lir), ncol(data)), colnames(data)) ## set T/F
  lirglim <- innerrange.limits  ## contains either a list of ranges or T/F
  if (length(lirglim)) {
    if (!is.list(lirglim))
      lirglim <- if (is.logical(lirglim)) as.list(lirglim) else list(lirglim) 
    if (length(lnm <- names(lirglim))) {
      lnm <- intersect(lnm, names(lirg))
      lirg[lnm] <- as.list(lirglim)
    } else lirg <- setNames(rep(lirglim, ncol(data)), colnames(data))
  }
  lirfactor <- i.getploption("innerrange.factor")
  ## varlabels
  llb <- NULL
  if (length(varlabels))
    if(!is.character(varlabels))
      warning("!genvarattributes! 'varlabels' must be of mode 'character'")
    else {
      if (length(lnm <- names(varlabels))==0) {
        if (length(varlabels)==ncol(data))
          llb <- setNames(varlabels, names(data))
        else warning("!genvarattributes! 'varlabels' must be of length  ",
                     ncol(data), "  or have names")
      }  else llb <- labels
    }
  ## variable names
  lnm <- colnames(data)
    if (anyNA(lnm)) colnames(data) <- lnm <-
      ifelse(is.na(lnm), paste("V",1:NCOL(data), sep=""), lnm)
  ##
  lrown <- row.names(data)
  ## Date
  for (lj in 1:ncol(data))
    if (inherits(ld <- data[,lj], "Date") && is.null(attr(ld, "numvalues")))
      data[,lj] <- gentimeaxis(setNames(ld,lrown))
  ## line color and type
  if (is.null(ynames))
    ynames <- union(union(names(ycol),names(ylty)),names(ypch))
  lny <- length(ynames)
  if (lImulty <- lny>0) {
    ldt <- data[,ynames, drop=FALSE]
    lypch <- 
      i.getvarattribute("pch", ypch, ldt, ploptions$variables.pch, drop=lny>1)
    lylty <- 
      i.getvarattribute("lty", ylty, ldt, ploptions$variables.lty, drop=lny>1)
    lycol <-
      i.getvarattribute("col", ycol, ldt, ploptions$variables.col, drop=lny>1)
  }
  ## --- loop
  for (lv in lnm) {
    lvv <- data[,lv]
    lcls <- class(lvv)[1]
    attr(lvv, "varname") <- lv
    lnv <- sum(!duplicated(lvv),na.rm=TRUE)
    ##
    if (replace || is.null(attr(lvv, "nvalues")) )
      attr(lvv, "nvalues") <- lnv
    ## turn character into factor
    if (lcls=="character") lvv <- factor(lvv)
    ##    if (lv %in% lfacgen)  class(lvv) <- c(class(lvv, "usedAsFactor"))
    if (lImulty) { ## line color and type, pch
      if(lv%in% names(lycol)) attr(lvv, "col") <- lycol[lv]
      if(lv%in% names(lylty)) attr(lvv, "lty") <- lylty[lv]
      if(lv%in% names(lypch)) attr(lvv, "pch") <- lypch[lv]
    }
    if (inherits(lvv, c("factor", "usedAsFactor"))) {
      ## factor
      lat <- seq_along(levels(lvv))
      if (replace || is.null(attr(lvv, "plrange")))
        attr(lvv, "plrange") <- c(0.5, max(lat)+0.5)
      if (replace || is.null(attr(lvv, "ticksat")))
        attr(lvv, "ticksat") <- lat
      if (replace || is.null(attr(lvv, "ticklabels")))
        attr(lvv, "ticklabels") <- levels(lvv)
      ## jitter
      if(replace || is.null(attr(lvv, "plcoord")) && (lij <- ljt[lv])) {
        attr(lvv, "plcoord") <-
          jitter(as.numeric(lvv), factor=jitter.factor,
                 amount=if(is.numeric(lij)) lij else NULL)
        attr(lvv, "plrange") <- c(0.5, length(levels(lvv))+0.5)
      }
    } else { ## continuous variable
      lvvv <- if(inherits(lvv,"Surv")) lvv[,1] else lvv
      names(lvvv) <- lrown
      lnmod <- c(0,0)
      if (is.numeric(lvvv)) {
        lrg <- NULL
        if (lir && u.notfalse(lirgv <- lirg[[lv]]) ) {
          ## inner range
          lrg <- attr(lvv, "innerrange")
          lplrg <- attr(lvv, "plrange")
          lIirg <- is.null(lrg) ## new range
          lIplrg <- is.null(lplrg) ## new range
          lrg <- if (length(lirgv)==2 && is.numeric(lirgv)) lirgv
                 else if (u.notfalse(lirgv))
                   plinnerrange(lirgv, lvvv, factor=lirfactor)
          if (replace | lIirg) {
            attr(lvv, "innerrange") <- lrg
            lpc <- plcoord(lvvv, lrg, ploptions=ploptions)
            ## attributes: avoid a level...
            lpca <- attributes(lpc)
            attributes(lvv)[names(lpca)] <- lpca
            if (!lIplrg) attr(lvv, "plrange") <- lplrg
            attributes(lpc) <- NULL
            attr(lvv, "plcoord") <- lpc
            lnmod <- lpca$nmod
            if (is.null(lnmod)) lnmod <- c(0,0)
            lrg <- attr(lvv, "innerrange") ## innerrange may have changed
            ## to avoid unnecessary inner bounds when plotting
          }
          if (!lIplrg) {
            lnmod <- c(attr(lvv, "nmod"), 0, 0)[1:2]
            lrg <- ifelse(lnmod>0, lrg, lplrg) 
          }
        } else { ## end inner range
          if (is.null(lrg))  
            lrg <- i.extendrange(range(lvvv, na.rm=TRUE), i.getploption("plext"))
        }
        ## set plrange
        if (replace || is.null(attr(lvv,"plrange")))
          attr(lvv, "plrange") <- lrg
        ## ticks
        ltint <- min(max(tickintervals[1],3)+sum(lnmod>0),20)
        lat <- clipat(pretty(lrg, n=ltint, min.n=ltint-2), lrg)
        llabat <-
          if (length(tickintervals)>1)
            clipat(pretty(lrg, n=tickintervals[2], min.n=1),lrg)  else lat
        ## not needed anymore:
##-         if (lnmod[1]) {
##-           lat <- lat[lat>=lrg[1]]
##-         }
##-         if (lnmod[2]) {
##-           lat <- lat[lat<=lrg[2]]
        ##-         }
        llabat <- clipat(llabat, lat)
        if (length(llabat)<2) {
          llabat <- clipat(pretty(lrg, n=4, min.n=1), lat)  
        } else lat
        if (replace  || is.null(attr(lvv, "ticksat"))) {
          attr(lvv, "ticksat") <- lat ## lat[lat>lrg[1]&lat<lrg[2]]
          ## do not replace  ticklabesat  if  ticksat  is set
          if (replace || is.null(attr(lvv, "ticklabelsat")))
            attr(lvv, "ticklabelsat") <- llabat ## lat[lat>lrg[1]&lat<lrg[2]]
        }
      }
    }
    attr(lvv,"varlabel") <- if (lv%in%names(llb)) llb[lv] else
      i.def(attr(lvv, "varlabel"), lv)
    attributes(data[[lv]]) <- attributes(lvv)
  }
  data
}
## =====================================================================
plframe <-
  function(x, y, xlab=NULL, ylab=NULL, ticklabels = TRUE, plextext=NULL,
           axcol=rep(1,4), plargs=NULL, ploptions=plargs$ploptions)
{
  ## -------------------
  lf.getcoord <- function(x, lext) {
    lx <- if (is.data.frame(x)) setNames(x[,1], row.names(x)) else x 
    if (length(lx)==0)
      stop("!plframe! unsuitable argument 'x' or 'y'")
    lxx <- i.def(i.def(attr(lx,"numvalues"), attr(lx,"plcoord")), lx)
    lrg <- attr(lx,"plrange")
    if (length(lrg)==0)
      lrg <- i.extendrange(range(lxx, na.rm=TRUE), lext[1:2])
    ltat <- attr(lx,"ticksat")
    if (length(ltat)==0) {
      if (is.factor(lx)) {
        llv <- levels(lx)
        attr(lx, "ticksat") <- seq_len(length(llv))
        attr(lx, "ticklabels") <- llv
      } else
        attr(lx, "ticksat") <- pretty(lrg, i.getploption("tickintervals")[1])
    }
    list(x = lx, xx = lxx, range = lrg)
  }
  ## ---
  lIlab <- as.logical(rep(ticklabels, length=4))
  lext <- rep(i.getploption("plext"),length=4)
  plextext <- rep(i.def(plextext, 0, i.getploption("plextext"), 0), length=4)
  axcol <- rep(i.def(axcol,1), length=4)
  ## plotting coordiates
  lcoord <- lf.getcoord(x, lext)
  lx <- lcoord$x
  lxx <- lcoord$xx
  lrgx <- lcoord$range
  lcoord <- lf.getcoord(y, lext[3:4])
  ly <- lcoord$x
  lyy <- lcoord$xx
  lrgy <- lcoord$range
##-   
##-   ly <- if (is.data.frame(y)) y[,1] else y ## x[,2] ## lpd[,y]
##-   lrgy <- attr(ly,"plrange")
##-   lyy <- i.def(attr(ly, "numvalues"), ly)
##-   if (length(lrgy)==0)
##-     lrgy <- i.extendrange(range(as.numeric(lyy), na.rm=TRUE), lext[3:4])
##-   if (length(ly)==0)
##-     stop("!plframe! unsuitable argument 'y'")
  ## 
  lmar <- rep(i.getploption("mar"), length=4)
  if (anyNA(lmar)) lmar <- ifelse(is.na(lmar), par("mar"), lmar)
  lmgp <- i.getploption("mgp")
  lop <- par(mar=lmar, mgp=lmgp)
  ## on.exit(par(lop))  ## produces artifact!
  plot(i.extendrange(lrgx, plextext[1:2]), i.extendrange(lrgy, plextext[3:4]),
       xlab = "", ylab = "", type="n", axes=FALSE, xaxs="i", yaxs="i")
  lmfg <- par("mfg")
  ## inner range
  lnmodx <- c(attr(lx, "nmod"),0,0)[1:2]
  lirgx <- c(attr(lx, "innerrange"), lrgx)[1:2]
  abline(v=unique(c(lirgx[lnmodx>0], lrgx)),lty=3)
  lnmody <- c(attr(ly, "nmod"),0,0)[1:2]
  lirgy <- c(attr(ly, "innerrange"), lrgy)[1:2]
  ## grid
  lgrl <- i.getploption("grid")
  ## !!! need attr("ticksat") -> otherwise, generate it!
  if (u.notfalse(lgrl)) {
    if (is.atomic(lgrl)) lgrlx <- lgrly <- lgrl
    else {
      if (length(lgnm <- names(lgrl))) {
        lgrlx <-
          if (length(lxnm <- attr(lx, "varname")) && lxnm%in%lgnm)
            lgrl[[lxnm]]  else  ploptions("grid")
        lgrly <-
          if (length(lynm <- attr(ly, "varname")) && lynm%in%lgnm)
            lgrl[[lynm]]  else  ploptions("grid")
      } else {
        if(length(lgrl)!=2) {
          warning(":plframe: unsuitable argument 'grid'")
          lgrlx <- lgrly <- ploptions("grid")
        } else {
          lgrlx <- lgrl[[1]]
          lgrly <- lgrl[[2]]
        }
      }
    }
    ## grid:x
    if (u.true(lgrlx)) {
      lgrlx <- attr(x,"ticksat")
      if (length(lgrlx)<2)
        lgrlx <- attr(x, "ticksat") <-
          pretty(lxx, i.getploption("tickintervals")[1]) 
    }
    lgrlx <- clipat(lgrlx, lirgx)
    if (length(lgrlx)) 
      abline(v=lgrlx,
             lty=i.getploption("grid.lty"),
             lwd=i.getploption("grid.lwd"),
             col=i.getploption("grid.col"))
    ## grid:y
    if (u.true(lgrly)) {
      lgrly <- attr(y,"ticksat")
      if (length(lgrly)<2)
        lgrly <- attr(ly, "ticksat") <-
          pretty(lyy, i.getploption("tickintervals")[1]) 
    }
    grl <- clipat(lgrly, lirgy)
    if (length(lgrly)) 
      abline(h=lgrly,
             lty=i.getploption("grid.lty"),
             lwd=i.getploption("grid.lwd"),
             col=i.getploption("grid.col"))
  }
  ## zero line
  if (length(lzl <- i.getploption("zeroline"))) {
    lzl <- i.def(lzl, 0, TRUE, NA)
    if (is.atomic(lzl)) lzlx <- lzly <- lzl
    else {
      lzlx <- lzl[[1]]
      lzly <- lzl[[2]]
    }
    if (is.logical(lzlx)) lzlx <- 0
    if (is.logical(lzly)) lzly <- 0
    abline(v=lzlx[lzlx>lirgx[1]&lzlx<=lirgx[2]],
           h=lzly[lzly>lirgy[1]&lzly<=lirgy[2]],
           lty=i.getploption("zeroline.lty"), lwd=i.getploption("zeroline.lwd"),
           col=i.getploption("zeroline.col"))
  }
  ## bounding boxes
  abline(h=unique(c(lirgy, lrgy)), lty=3)
  lxrg <- ifelse(lnmodx>0, lirgx, lrgx)
  lyrg <- ifelse(lnmody>0, lirgy, lrgy)
  lines(lxrg[c(1,2,2,1,1)], lyrg[c(1,1,2,2,1)])
  ## axes
  laxes <- i.getploption("axes")
  if (length(xlab)) attr(lx, "varlabel") <- xlab
  if (length(ylab)) attr(ly, "varlabel") <- ylab
  if (u.notfalse(laxes)) {
    ltint <- i.getploption("tickintervals")[1]
    if (1%in%laxes)
      plaxis(1, lx, lab = lIlab[1] && (lmar[1]>=lmgp[2]+1 | lmfg[1]==lmfg[3]),
             range=lxrg, col=axcol[1], tickintervals=ltint)
    if (2%in%laxes)
      plaxis(2, ly, lab = lIlab[2] && (lmar[2]>=lmgp[2]+1 | lmfg[2]==1), 
             range=lyrg, col=axcol[2], tickintervals=ltint)
    if (3%in%laxes)
      plaxis(3, lx, lab = lIlab[3] && (lmar[3]>=lmgp[2]+1 | lmfg[1]==1), 
             range=lxrg, col=axcol[3], tickintervals=ltint)
    if (4%in%laxes)
      plaxis(4, ly, lab = lIlab[4] && (lmar[4]>=lmgp[2]+1 | lmfg[2]==lmfg[4]),
             range=lyrg, col=axcol[4], tickintervals=ltint)
  }
  invisible(lop)
}
## --------------------------------------------------------------------
plaxis <-
  function(axis, x, lab=TRUE, range=NULL, varlabel=NULL, col=1,
           tickintervals=NULL, plargs = NULL, ploptions = plargs$ploptions,
           ...)
{
  lx <- i.def(i.def(attr(x,"numvalues"), attr(x,"plcoord")), x)
  range <- i.def(range, i.def(attr(x, "innerrange"), attr(x, "plrange")),
              valuefalse = range(lx,na.rm=TRUE))
  lat <- attr(x,"ticksat")
  latsmall <- attr(lat, "small")
  llabat <- attr(x,"ticklabelsat")
  llab <- attr(x,"ticklabels")
  varlabel <- i.def(i.def(varlabel, attr(x,"varlabel"), valuefalse=""),
                    attr(x,"varname") )
  lmar <- par("mar")
  lmgp <- par("mgp")
  lmfg <- par("mfg")
  col <- i.def(col, 1)
  lcex <- par("cex")
  lIouter <- switch(axis, lmfg[1]==lmfg[3], lmfg[2]==1,
                    lmfg[1]==1, lmfg[2]==lmfg[4])
  if(lab && (lmar[axis]>lmgp[1]+1 | lIouter) )
    mtext(varlabel, side=axis, line=lmgp[1], xpd=TRUE, col=col,
          cex=lcex*par("cex.lab"))
  ## tick labels
  if (length(lat)>1) {
    lat <- lat[lat>=range[1] & lat<=range[2]]
    if (length(lat)>1) {
      if (length(llabat)) {
        if (!is.numeric(llabat)) {
          warning(":plaxis: 'ticklabelsat' must be numeric")
          llabat <- NULL
        }
      } else llabat <- lat
    }
    if (length(llab)) {
      if (length(llab)!=length(llabat)) {
        if (length(llab)==1) llab <- rep(llab,length(llabat))
        else {
          warning(":plaxis: 'ticklabel' has wrong length")
          llab <- NULL
        }
      }
    }
    if (length(llab)==0) llab <- if (lab) format(llabat)
  } else {
    lat <- pretty(range, i.getplopt(tickintervals))
    llabat <- lat <- lat[lat>=range[1] & lat<=range[2]]
    llab <- if (lab) format(lat) else FALSE
  }
  if (length(llab)==0 | !lab) llab <- rep("",length(lat))
  ##
  axis(axis, at=lat, labels=rep("",length(lat)), col=col, ...)
  if (length(llabat))
    mtext(llab, axis, line=i.getploption("mgp")[2], at=llabat, col=col,
          cex=lcex*par("cex.axis"), ...)
  if (length(latsmall))
    axis(axis, at=latsmall, labels=rep("",length(latsmall)),
         tcl=0.5*par("tcl"), col=col, ...)
}
## ----------------------------------------------------------------------
pltitle <-
  function(main=NULL, sub=NULL, cex=NULL, cexmin=NULL, 
           side=3, line=NULL, adj=NULL, outer.margin=NULL, col="black",
           doc=NULL, show=TRUE, plargs=NULL, ploptions=plargs$ploptions, ...)
{
  ## Purpose:   title
  ## ----------------------------------------------
  main <- i.def(main, plargs$main)
  sub <- i.def(sub, plargs$sub, valuefalse=NULL) ## paste(":",plargs$sub,sep="")
  lsub <- if (length(sub) && substr(sub,1,1)==":")
            paste(substring(sub,2,30), plargs$.subdefault, sep="")
          else as.character(sub)
  if (length(main)==0 & length(lsub)==0) return()
  doc <- doc ## options("doc") !!!
  if (length(main) && all(lsub==main)) lsub <- NULL
  rr <- c(main=main, sub=sub)
  ##
  lf.text <- function(text, cex, cexdef, ...) {
    lcex <-
      i.def(cex, max(cexmin, min(cexdef, lfac/max(nchar(text)))),
            valuefalse = 0 )      
    lmaxchar <- lfac/lcex
    if (any(li <- nchar(text)>lmaxchar))
      text[li] <- paste(substr(text[li], 1, lmaxchar-3),"...")
    ladj <- i.def(adj, max(minadj,0.5*(lcex>cexmin)), 0.5, minadj)
    mtext(text, side, line-lcex, cex = lcex, adj=ladj, outer = outer.margin,
          col=col, ...)
    lcex
  }
  show <- i.def(show, 0.5)
  outer.margin <- i.def(outer.margin, par("oma")[3]>0,
                        valuefalse = FALSE)
  if (show<=0 || (show<1 && outer.margin && 1!=prod(par("mfg")[1:2])) )
    return(rr)
  lcexdef <- rep(i.getploption("title.cex"), length=3)
  title.cexmin <- cexmin
  cexmin <- i.getplopt(title.cexmin)
  scale <- 2.5
  minadj <- 0.2
  ##
  lwid <- if (outer.margin) par("mfg")[3] else 1
  lfac <- lwid * scale*par("pin")[1]/par("cin")[1]
  line <-
    i.def(line, max(if (outer.margin) par("oma")[side] else par("mar")[side]) )
  ##
  if (length(main) && main!="") {
    lcex <- lf.text(main, cex=cex[1], cexdef=lcexdef[1], ...)
    line <- line-lcex
  }
  if (length(lsub) && as.character(lsub)!=":") {
    lcex <- lf.text(lsub, cex=cex[2], cexdef=lcexdef[2], ...)
    line <- line-lcex
  }
  if (line>=0 && (!is.null(doc)) && doc && length(tit(main)))
    lf.text(tit(main), cex=cex[3], cexdef=lcexdef[3], ...)
  invisible(rr)
}
## -----------------------------------------------------------------
plpoints <-
  function(x = NULL, y = NULL, type = "p",
           plab = NULL, pch = NULL,
           col = NULL, lty = NULL, lwd = NULL, psize = NULL,
           condquant = TRUE,
           plargs = NULL, ploptions = plargs$ploptions, ...)
{
  if (is.null(plargs)) 
    plargs <- get(".plargs", ".GlobalEnv") ## list in calling fn
  pldata <- plargs$pldata
  ## intro, needed if formulas are used or data is given or ...
  lcl <- match.call()
  lcall <- sys.call() ## match.call modifies argument names
  lcnm <- i.def(names(lcall), names(lcl))
  names(lcall) <- ifelse(lcnm=="", names(lcl), lcnm)
  if (is.formula(x)|is.formula(y)|any(c("data","pcol")%in%lcnm)) {
    ## !!! this does not (always) work if plpoints is called by a call
    lcall$assign <- FALSE
    lcall$ploptions <- ploptions
    lplargs <- do.call(pl.control, as.list(lcall[-1]), envir=parent.frame())
    ploptions <- lplargs$ploptions
    plargs$pldata <- pldata <- lplargs$pldata
    x <- pldata[,2]
    y <- pldata[,1]
  } else {
  ## ---
    if (length(x)==0) x <- pldata[,2]
    if (length(y)==0) y <- pldata[,1]
  }
  ## ---
  if (is.data.frame(x))  x <- x[,1]
  lattrx <- attributes(x)
  lx <- i.def(lattrx$plcoord, i.def(lattrx$numvalues, x))
  if (is.data.frame(y))  y <- y[,1]
  lattry <- attributes(y)
  ly <- i.def(lattry$plcoord, i.def(lattry$numvalues, y))
  ## condquant
  condquant <- i.def(condquant, 1, 1, 0) 
  lIcq <- length(condquant)>0  ## condquant representation by bars
  lIcqx <- length(lcqx <- lattrx$condquant) >0
  lIcqy <- length(lcqy <- lattry$condquant) >0
  lIcensx <- inherits(x, "Surv")
  lIcensy <- inherits(y, "Surv")
  ##
  lnr <- length(x)
  pldata <- as.data.frame(plargs$pldata) ## 'as.data.frame' only needed for
  ## (not) finding the () columns if plargs is NULL
  psize <- i.def(psize, pldata[["(psize)"]])
  lpsize <-
    if (is.null(psize)) 1  else  sqrt(psize/median(psize, na.rm=TRUE))
  if (is.null(plab)) plab <- pldata[["(plab)"]]
  if (is.null(pch))
    pch <- i.def(i.def(pldata[["(pch)"]], lattry$ypch),
                 i.getploption("pch"), valuefalse="")
  pch <- rep(pch, length=lnr)
  lpcol <- i.def(col, pldata[["(pcol)"]])
  if (length(lpcol)) {
    lgrpcol <- i.getploption("group.col")
    if (is.factor(lpcol)) lpcol <- as.numeric(lpcol)
    if (is.numeric(lpcol)) lpcol <- rep(lgrpcol, length=max(lpcol))[lpcol]
  } else lpcol <- i.getploption("col")[1]
  pcol <- rep( lpcol, length=lnr)
  lty <- i.def(lty, i.getploption("lty"))
  lwd <- i.def(lwd, i.getploption("lwd")) * i.getploption("linewidth")[lty]
  lcol <- if (type=="l") i.getploption("lcol")[1] else lpcol
  lnobs <- sum(is.finite(lx) & is.finite(ly))
  cex <- i.getploption("cex")
  cex <- if (is.function(cex)) cex(lnobs) else i.def(cex, cexSize(lnobs))
##  cex.pch <- cex*i.getploption("default.cex")
  cex.plab <- cex*i.getploption("cex.plab")
  lsplab <- rep(abs(cex.plab*lpsize), length=lnr) 
  lspch <- rep(cex*lpsize, length=lnr)
  ## censored
  if (lIcensx|lIcensy) {
    lstx <- if (lIcensx) (1+(i.def(lattrx$type, "right")=="left"))*(1-x[,2]) else 0
    lsty <- if (lIcensy) (1+(i.def(lattry$type, "right")=="left"))*(1-y[,2]) else 0
    lipch <- lstx+3*lsty
    li <- lipch>0
    if (any(li)) {
      pch[li] <- rep(i.getploption("censored.pch"), length=8)[lipch[li]]
      pcol[li] <- colorpale(pcol[li], i.getploption("censored.pale"))
      lspch[li] <- lspch[li]*i.getploption("censored.size")
    }
  }
  ## condquant
  if (lIcq) {
    lpale <- rep(c(i.getploption("condquant.pale"), 0.5), length=2)
    lix <- if(lIcqx) lcqx[,"index"]
    liy <- if(lIcqy) lcqy[,"index"]
    lixy <- union(lix,liy)
    if (length(lixy)==lnr) lpale <- c(1,lpale[1]) ## all observations are cq
    pcol[lixy] <- colorpale(pcol[lixy], lpale[1])
    if (lIcqx) {
      li <- lix %nin% liy
      if (any(li)) {
        lsg <- if(length(lrgx <- lattrx$innerrange))
                 plcoord(lcqx[li,2:3], range=lrgx, ploptions=ploptions)
               else lcqx[li,2:3]
        lyy <- ly[lix[li]]
        segments(lsg[,1], lyy, lsg[,2], lyy,
                 col=colorpale(pcol[lix[li]], lpale[2]) )
      }
    }
    if (lIcqy) {
      li <- liy %nin% lix
      if (any(li)) {
        lsg <- if(length(lrgy <- lattry$innerrange))
                 plcoord(lcqy[li,2:3], range=lrgy, ploptions=ploptions)
               else lcqy[li,2:3]
        lxx <- lx[liy[li]]
        segments(lxx, lsg[,1], lxx, lsg[,2],
                 col=colorpale(pcol[liy[li]], lpale[2]) )
      }
    }
  }
##-     if (condquant>0) {
##-       if (lIcqx) 
##-         x[lcqix <- lcqx[,"index"]] <- NA
##-       if (lIcqy) 
##-         y[lcqiy <- lcqy[,"index"]] <- NA
##-       if (lIcqx) plbars(lcqx[,1:3], y[lcqix], plargs=plargs, ploptions=ploptions)
##-       if (lIcqy) {
##-         if(length(lrgy <- lattry$innerrange))
##-           lcqy <- plcoord(lcqy[,1:3], range=lrgy, ploptions=ploptions)
##-         plbars(x[lcqiy], lcqy, plargs=plargs, ploptions=ploptions)
##-       }
##-     } else {
##-       pch[lcqx[,"index"]] <- lpchcq[1]
##-       pch[lcqy[,"index"]] <- lpchcq[2]
##-       pch[intersect(lcqx[,"index"],lcqy[,"index"])] <- lpchcq[3]
##-     }
  lIpl <- (length(plab)>0) && any(!is.na(plab))
  plab <- as.character(plab)
  ## --- plot!
  if (lIpl) {
    text(lx, ly, plab, cex=lsplab*lpsize, col=pcol)
    lipch <- ifelse(is.na(plab), TRUE, plab=="")
  } else lipch <- rep(TRUE,length(x))
  if (any(lipch)) {
    if (type=="l") lines(lx, ly, col=lcol, lty=lty, lwd=lwd)
    else {
      if (type=="b") 
        points(lx, ly, pch=NA, col=lcol, cex=lspch, type="b", lty=lty, lwd=lwd)
      points(lx[lipch], ly[lipch], pch=pch[lipch], cex=lspch[lipch],
             col=pcol[lipch])
    }
  }
}
pllines <- function(x=NULL, y=NULL, type="l", ...)
  plpoints(x, y, type, ...)

## -----------------------------------------------------------------
plbars <-
  function(x, y, plargs = NULL, ploptions = plargs$ploptions)
{
  lcol <- i.getploption("bar.col")
  llty <- rep(i.getploption("bar.lty"), length=2)
  llwd <- rep(i.getploption("bar.lwd"), length=2)
  lmpw <- i.getploption("bar.midpointwidth")*
    diff(par("usr"))[c(1,3)]/100 
  if (is.matrix(x)) {
    segments(x[,1], y-lmpw[2], x[,1], y+lmpw[2], lty=llty, lwd=llwd[1],
             col=lcol)
    segments(x[,2], y, x[,3], y, lty=llty, lwd=llwd[2], col=lcol)
  }
  if (is.matrix(y)) {
    segments(x-lmpw[1], y[,1], x+lmpw[1], y[,1], lty=llty, lwd=llwd[1],
             col=lcol)
    segments(x, y[,2], x, y[,3], lty=llty, lwd=llwd[2], col=lcol)
  }
}
## -----------------------------------------------------------------
plmark <-
  function(x, y=NULL, markextremes=NULL, plabel=NULL,
           plargs=NULL, ploptions=plargs$ploptions)
{
  lf.pmark <- function(mprop, x) {
    lrk <- (rank(x, na.last="keep")-0.5)/sum(is.finite(x))
    lrk<mprop[1] | lrk>1-mprop[2]
  }
  ##
  if (is.null(plargs)) 
    plargs <- get(".plargs", ".GlobalEnv") ## list in calling fn
  plabel <- i.def(i.def(plabel, plargs$plabel), .plargs$plabel)
  if (is.null(plabel)) {
    warning(":plmark: no plabels found")
    return(NULL)
  }
  lx <- if (is.data.frame(x)) x[,1] else x 
  ly <- if (is.data.frame(y)) y[,1] else y 
  lnobs <-
    if (is.null(ly)) sum(is.finite(lx)) else sum(is.finite(lx)&is.finite(ly))
  lmxdef <- ceiling(sqrt(lnobs)/2)/lnobs
  lmx <- markextremes
  if (is.atomic(lmx)) 
    lmxx <- lmxy <- i.def(lmx, lmxdef)
  else {
    if (length(lmxnm <- names(lmx))) {
      lmxx <- i.def(unlist(lmx[attr(x,"varname")]), lmxdef)
      lmxy <- i.def(unlist(lmx[attr(y,"varname")]), lmxdef)
    } else if(length(lmx)==2) {
      lmxx <- lmx[[1]]
      lmxy <- lmx[[2]]
    } else {
      warning(":plmark: unsuitable argument 'markextremes'")
      lmxx <- lmxy <- 0
    }
  }
  lmxx[is.na(lmxx)] <- lmxdef
  lmxy[is.na(lmxy)] <- lmxdef
  rr <- if (any(c(lmxx, lmxy)>0)) {
    li <- lf.pmark(rep(lmxx, length=2), lx)
    if (length(y)) li <- li | lf.pmark(rep(lmxy, length=2), ly)
    ifelse(li, plabel, "")
  } else rep("", length(lx)) ## plabel
  invisible(rr)
}
## -----------------------------------------------------------------
plsmooth <-
  function(x, y, band=NULL, power = NULL, plargs=NULL, ploptions=plargs$ploptions)
{
  lIsm <- i.getploption("smooth")
  lsm <- NULL
  if (lIsm) {
    power <- i.def(power, 1,1,1)
    band <- i.def(i.def(band, lIsm>=2), FALSE, TRUE, FALSE)
    lsm <- gensmooth(x, y, band=band, power=power,
                     plargs=plargs, ploptions=ploptions)
    plsmoothline(lsm, x, y, plargs=plargs, ploptions=ploptions)
  }
  invisible(lsm)
}
## --------------------------------------------------------------
plsmoothline <-
  function(smoothline, x, y, plargs = NULL, ploptions = plargs$ploptions)
{
  if (!is.list(smoothline)) {
    warning(":plsmoothline: 'smoothline' is not suitable. No smooth lines")
    return()
  }
  lxtrim <- i.getploption("smooth.xtrim")
  lrgx <- attr(x,"innerrange")
  lx <- smoothline$x
  ly <- as.matrix(smoothline$y)
  lgrp <- as.numeric(smoothline$group)
  if (lInogrp <- length(lgrp)==0 || length(unique(notna(lgrp)))<=1) {
    lgrp <- rep(1, NROW(lx))
    lngrp <- 1
    lcol <- i.getploption("smoothline.col")
    llty <- i.getploption("smoothline.lty")
  } else {
    lngrp <- max(lgrp)
    lcol <- i.getploption("group.col")
    llty <- rep(i.getploption("group.lty"), length=lngrp)
  }
##  llty <- rep(llty, each=2) 
  lbd <- smoothline$yband
  lIband <- length(lbd)>0
  ## check if ordered
  lio <- order(lgrp, lx) 
  if (any(lio!=1:length(lio))) {
    lx <- lx[lio]
    ly <- ly[lio,, drop=FALSE]
    lgrp <- lgrp[lio]
    lngrp <- max(lgrp)
    if (lIband) {
      lbd <- lbd[lio]
      smoothline$ybandindex <- smoothline$ybandindex[lio]
    }
  }
  lx[lx<lrgx[1]|lx>lrgx[2]] <- NA
  lny <- ncol(ly)
  ## may be a 2-vector or  a matrix of 2 rows
  llwd <- rep(c(i.getploption("smoothline.lwd"), 0.7), length=2)
  llwid <- i.getploption("linewidth")
  lpale <- i.getploption("smoothline.pale")
  for (lgr in seq_len(lngrp)) {
    lig <- which(lgrp==lgr)
    lxg <- lx[lig]
    if (1< (lng <- length(lig))) {
      lndr <- round(lng * if(is.function(lxtrim))
                            lxtrim(lng) else i.def(lxtrim, 0,
                                                   smoothxtrim(lng), 0) )
      if (lndr) lxg[- ((lndr+1):(lng-lndr)) ] <- NA
      lcl <- lcol[min(lgr,length(lcol))]
      if (lny>1) 
        matlines(lxg, ly[lig,-1], lty=llty[lgr],
                 lwd=llwd[1]*llwd[2]*llwid[llty[lgr]],
                 col = colorpale(lcl, lpale))
      llw <- llwd[1]*llwid[llty[lgr]]
      lines(lxg, ly[lig,1], lty=llty[lgr], lwd=llw, col=lcl)  ## xxx
      if (lIband) {
      li <- smoothline$ybandindex[lig] ## separate upper and lower smooth
      if (any(li))
        lines(lxg[li], lbd[lig[li]], lty=llty[lgr], lwd=llw/2, col = lcl) 
      if (any(!li))
        lines(lxg[!li], lbd[lig[!li]], lty=llty, lwd=llw/2, col = lcl) 
    }
  }}
}
## -----------------------------------------------------------------
plrefline <-
  function(refline, x=NULL, innerrange=NULL, y=NULL,
           cutrange = c(x=TRUE, y=FALSE),
           plargs=NULL, ploptions=plargs$ploptions)
{
  ## draws a reference line (with extended range) and
  ##   band given by reflineyw (only inner range) if requested
  lf.irna <- function(x, rg) {
    x[x<rg[1]|x>rg[2]] <- NA
    x
  }
  lrfyb <- NULL
  if (is.function(refline))
    refline <- refline(y~x)
  if (is.list(refline)&&length(refline$coef)) refline <- refline$coef
  if (is.atomic(refline)) {
    if (length(refline)!=2) {
      warning(":plrefline: 'refline' not suitable. No reflines")
      return()
    }
    lusr <- par("usr")
    lrfx <- seq(lusr[1],lusr[2],
                length=i.getploption("functionxvalues"))
    lrfy <- refline[1]+refline[2]*lrfx ## needs correction if lIxir
  } else {
    lrfx <- refline$x
    lrfy <- refline$y
    if (length(lrfx)==0|length(lrfx)!=NROW(lrfy)) {
      warning(":plrefline: 'refline' not suitable. No reflines")
      return()
    }
    lrfyb <- refline$band
  }
  lIrfyb <- length(lrfyb)
  ## ---
  if (is.null(innerrange)) innerrange <- attr(x, "innerrange")
  ## haul refline to inner plotrange
  if (length(innerrange)>0) {
    if (is.null(x)) x <- plargs$pldata[,2] ## !!!
    if (is.null(y)) y <- plargs$pldata[,1]
    lIxir <- any(attr(x,"nmod")>0)
    lxir <- attr(x,"innerrange")
    lIyir <- any(attr(y,"nmod")>0)
    lyir <- attr(y,"innerrange")
    cutrange <- rep(i.def(cutrange, TRUE), length=2)
    if (lIxir) {
      if (i.def(cutrange[1], TRUE)) lrfx <- lf.irna(lrfx, lxir)
      else
        lrfx <- plcoord(lrfx, range=lxir, ploptions=ploptions)
    }
    if (lIyir) {
      if (i.def(cutrange[2], FALSE)) lrfy <- lf.irna(lrfy, lyir)
      else
      lrfy <- plcoord(lrfy, range=lyir, ploptions=ploptions)
      if (lIrfyb)  ## band: values outside inner range -> NA
        lrfb <- apply(as.matrix(lrfb),2, lf.irna, rg=lyir)
    }
  }
  ## draw reference lines
  llty <- rep(i.getploption("refline.lty"), length=2)
  llwd <- rep(i.getploption("refline.lwd"), length=2)
  lcol <- rep(i.getploption("refline.col"), length=2)
  matlines(lrfx, lrfy, lty=llty[1], col=lcol[1], lwd=llwd[1])
  if (lIrfyb)
    matlines(lrfx, as.matrix(lrfy+lrfyb), lty=llty[2], lwd=llwd[2], col=lcol[2])
}
## =========================================================================
plcoord <-
  function(x, range=NULL, innerrange.factor=NULL, innerrange.ext=NULL,
           plext=NULL, ploptions=NULL)
{
  ## Purpose:    values for plot with limited "inner" plot range
  ldtrg <- range(x, na.rm=TRUE)
  lirfunc <- i.getploption("innerrange.function")
  if (is.character(lirfunc)) lirfunc <- get(lirfunc)
  innerrange.factor <- i.getplopt(innerrange.factor)
  innerrange.ext <- i.getplopt(innerrange.ext)
  plext <- i.getplopt(plext)
  if (length(notna(range))>0) {
    lrg <- range(range, na.rm=TRUE)
    if (lrg[1]>ldtrg[2] | lrg[2]<ldtrg[1]) { ## ranges do not overlap
      warning(":plcoord: inadequate range. Not used")
      lrg <- NULL
    } 
  } else lrg <- NULL
  if (length(lrg)==0) lrg <- lirfunc(x, fac=innerrange.factor)
  ## inner range must not extend beyond data
  if (length(lrg)==0) lrg <- ldtrg
  if (diff(lrg)==0) lrg <- c(-1,1)*lrg
  lrgext <- i.extendrange(lrg, plext)
  if (ldtrg[1]>=lrgext[1]) lrg[1] <- ldtrg[1]
  if (ldtrg[2]<=lrgext[2]) lrg[2] <- ldtrg[2]
  ## --------- transformation of data into plcoord
  rr <- pmax(pmin(x,lrg[2]),lrg[1])
  lxd <- x-rr
  lnmod <- c(sum(lxd<0,na.rm=TRUE),sum(lxd>0,na.rm=TRUE))
  if (sum(lnmod)>0) rr <- rr+lxd/(1+abs(lxd)/(diff(lrg)*innerrange.ext))
  ## if data fits into extended inner range, then avoid inner range
  lrg <- ifelse(lnmod>0, lrg, ldtrg)
  ## ---------
  ## extend range to plotting range
  lplrg <- i.extendrange(lrg, ifelse(lnmod>0, innerrange.ext, plext))
  ## extend inner range if there are no modified points
  lrg <- ifelse(lnmod>0, lrg, lplrg) ## adjust to the needs of plotting 
  attr(rr,"innerrange") <- lrg
  attr(rr,"innerrange.ext") <- innerrange.ext
    ## needed for transforming further quantities
  attr(rr,"nmod") <- lnmod
  attr(rr,"plrange") <- lplrg
  class(rr) <- class(x)
  rr
}
## -------------------------------------------------------------------
gentimeaxis <-
  function(date=NULL, year=NULL, month=NULL, day=NULL, hour=NULL, data=NULL,
           ploptions=NULL)
{
  f.2char <-
    function(x)
      substring(ifelse(x<10, paste("0",x,sep=""), as.character(x)),1,2)
  ## --- arguments
  lcall <- match.call()
  ltrangelim <- i.getploption("timerangelim")
  ltickint <- i.getploption("tickintervals")
  ldate <- NULL
  ## arguments can be names of variables in  data
  lnm <- NULL
  if (length(data)) { ## get years, months, ... from  data
    data <- as.data.frame(data)
    lnm <- row.names(data)
    largs <- setdiff(names(lcall), c("", "data"))
    ## arguments that are names of variables
    largs <- largs[sapply(as.list(lcall[largs]),
                          function(x) is.character(x) & length(x)==1 )]
    if (length(largs)) {
      liargs <- match(largs, names(data))
      if (any(lj <- is.na(liargs)))
        stop("!gentimeaxis! variable ",paste(liargs[lj], collapse=", "),
             " not found")
      ldf <- data[,liargs, drop=FALSE]
      names(ldf) <- largs
      if ("date"%in%largs) {
        ldate <- lcall$date  ## ldate contains the name of the date variable
        date <- data[,ldate]
      }
      if ("year"%in%largs) year <- ldf[,year]
      if ("month"%in%largs) month <- ldf[,month]
      if ("day"%in%largs) day <- ldf[,day]
      if ("hour"%in%largs) hour <- ldf[,hour]
    }
  }
  ## identify months -> numeric
  if (is.factor(month)) month <- as.character(month)
  if (is.character(month)) {
    lmnum <- match(month, c.months)
    if (any(is.na(lmnum))) lmnum <- match(month, c.mon)
    if (any(is.na(lmnum))) stop("!gentimeaxis! inadequate argument 'month'")
    month <- lmnum
  }
  ##
  if (lIdate <- length(date)) {
    if (is.character(date)) date <- as.POSIXct(date)
    if (all(is.na(date))|!inherits(date, c("Date","POSIXt"))) {
      warning(":gentimeaxis: argument 'date' not suitable. It is returned as is.")
      return(date)
    }
    ldt <- format(date)
    year <- as.numeric(substring(ldt,1,4))
    month <- as.numeric(substring(ldt,6,7))
    day <- as.numeric(substring(ldt,9,10))
    if (inherits(date,"POSIXt")& length(hour)==0)
      hour <- as.numeric(substring(ldt,12,13))
  }
##-   else {
##-     if (length(year)) {
##-       date <-
##-         as.Date(paste(year,
##-                       if (is.null(month)) "01" else f.2char(month),
##-                       if (is.null(day)) "01" else f.2char(day), sep="-") 
##-     } }
  if (length(year)==0) year <- 2000 else
    if (is.null(lnm)) lnm <- names(year)
  if (length(month)==0) month <- 6.5 else
    if (is.null(lnm)) lnm <- names(month)
  if (length(day)==0) day <- 15  else
    if (is.null(lnm)) lnm <- names(day)
  if (length(hour)==0) hour <- 12  else
    if (is.null(lnm)) lnm <- names(hour)
  ## lcat contains the calendar elements that are available
  lcat <- list(year = unique(year), month = unique(month), day = unique(day),
               hour = unique(hour))
  ## 
  lx <-
    julian(as.Date(paste(year,f.2char(month),f.2char(day),sep="-")),
           origin = as.Date("2000-01-01")) + hour/24
  lndays <- diff(range(lx, na.rm=TRUE))
  ##
  lcvar <- sapply(lcat, length) > 1 ## those that vary
  lat2 <- NULL
  llab <- NULL
  llabel <- NULL ## variable label
  ## -- years vary
  if (lcvar["year"]) {
    lyr <- min(year):max(year)
    if (lcvar["month"]) { ## months vary
      if (!lcvar["day"]) day <- 1
      lystart <- julian(as.Date(paste(lyr,"01-01", sep="-")),
                        origin = as.Date("2000-01-01"))
      lmstart <-
        julian(as.Date(paste(rep(lyr, each=12),f.2char(1:12),"01", sep="-")),
               origin = as.Date("2000-01-01"))
      if (lndays > ltrangelim["year"]) { ## tick labels at years only
        latyr <-
          if (lndays > 365*ltickint[1]) pretty(year, n=ltickint[1])
          else lcat[["year"]] 
        liy <- match(latyr,lyr)
        lat1 <- lystart[liy] ## ticks
        llab <- as.character(latyr) ## labels
        llabat <- lat1 + 365/2 ## labels in middle of year
        ## secondary ticks
        lat2 <- if (length(lyr)>length(latyr)) lystart 
                else  lmstart[seq(1,length(lmstart),3)] ## ticks at quarters
      } else { ## tick labels at years and quarters
        lat1 <- lmstart[seq(1,length(lmstart),3)] ## ticks
        lyr <- paste(lcat[["year"]],":",sep="")
        llab <- paste(c(rbind(lyr,"","","")),
                      rep(c.quarters, length=length(lat1)), sep="")
        llabat=lat1 + 30*rep(c(0, 0.5, 0.5, 0.5), length=length(lat1) )
        lat2 <- lmstart
      }
    } else { ## years only: unit is year
      lx <- year
      latyr <- pretty(lyr, n=ltickint)
      lat1 <- latyr[latyr==round(latyr)] ## only integers
      llab <- as.character(lat1)
      llabat <- lat1+0.5
    }
  } else {
    ## months vary, years do not
    if (lcvar["month"]) {
      lmon <- lcat[["month"]]
##-       lx <-
##-         julian(
##-           as.Date(paste(lcat[["year"]],f.2char(month),f.2char(day),sep="-")),
##-                   origin = as.Date("2000-01-01"))
      lmstart <-
        julian(as.Date(paste(lcat[["year"]],f.2char(lmon),"01", sep="-")),
                       origin = as.Date("2000-01-01"))
      lndays <- diff(range(lx, na.rm=TRUE))
      ##if (lndays > timerangelim["year"]) { ## tick labels at years only
      lat1 <- lmstart-0.5 ## start of the month
      llab <- c.mon[lmon] 
      llabat <- lat1+15  ## middle 
      lat2 <- rep(lat1,each=2) + rep(c(10,20), length(lat1)) ## 2 additional t
    } else {
      if (lcvar["day"]) {
        if (lcvar["hour"]) {
          lx <- day+hour/24
          lndays <- diff(range(lx, na.rm=TRUE))
          lat1 <- lcat[["day"]]
          lat2 <- if(lndays<ltrangelim["day"])
                    rep(lat1, each=3)+rep(1:3/4, length(lat1))  else NULL
          llabat <- lat1+0.5
          llab <- as.character(lat1)
          llabel <- "day"
        } else {
          return(day)
        }
      } else { ## only hour varies
        lx <- hour
        lat1 <- c(0,6,12,18,24)
        lat2 <- 0:24
        llab <- paste(c(0,6,12,18,24))
        llabat <- lat1
        llabel <- NULL
      }
    }
  }
  lat1 <- c(lat1, max(lat1)+max(diff(lat1))) ## lazy for months...
  if (length(lnm)==length(lx)) names(lx) <- lnm
  ## rr <- lx ## if(length(date)==0) lx else structure(as.Date(date), plcoord=lx)
  structure(if (length(ldate)) data[,ldate] else lx,
            numvalues=if(lIdate) lx else NULL, 
            ticksat=structure(lat1, small=lat2),
            ticklabelsat=llabat, ticklabels=llab, varlabel=llabel)      
}
## =====================================================================
i.pchcens <-
  function(plargs, condquant)
    ##  Delta, nabla, >, <, quadrat : pch= c(24, 25, 62, 60, 32)
{
  if (is.null(condquant) | !is.null(lpc <- plargs$pldata$"(pch)"))
  return(lpc)
  ##
  ploptions <- plargs$ploptions
  lpch <- i.getploption("censored.pch")
  lpc <- rep(i.def(lpch[1],1), nrow(plargs$pldata))
  lpc[condquant[,"index"]] <- i.def(lpch[2],3)
  lpc
}

## ====================================================================
pllimits <-
  function(pllim, data, limfac = NULL, FUNC=NULL)
{ ## determine inner plot range
  ## if pllim is a list or a matrix, leave it alone
  lIcq <- inherits(data, "condquant")
  ldt <- cbind( if (lIcq) c(data[,1:3]) else data )
  pllim <- i.def(pllim, TRUE, TRUE, FALSE)
  if (length(pllim)>1)
    if (any(dim(cbind(pllim))!=c(2,NCOL(ldt)))) {
      warning(":plot.regr/pllimits: unsuitable argument  pllim ")
      pllim <- TRUE
    }
  lfunc <- i.getploption("innerrange.function")
  limfac <- i.getploption("innerrange.factor")
  ##
  if (length(pllim)==1 && is.logical(pllim))
    pllim <-
      if (pllim) apply(ldt, 2, lfunc, fac=limfac) else
      matrix(FALSE, 2, NCOL(ldt))
  if ((!is.list(pllim))&length(pllim)==2) pllim <- as.matrix(pllim)
  pllim
}
## ========================================================================
plsubset <-
  function(x, subset)
{
  if (!(is.logical(subset)|is.numeric(subset)))
    stop("!plsubset! 'subset' must be logical or numeric")
  if (!inherits(x, "data.frame"))
    stop("!plsubset! 'x' must inherit from data.frame")
  lx <- x[subset,,drop=FALSE]
  for (lj in seq_len(ncol(x)))
    if (length(lattr <- attributes(x[,lj]))) {
      if ("numvalues"%in%names(lattr))
        lattr$numvalues <- lattr$numvalues[subset]
      if ("plcoord"%in%names(lattr))
        lattr$plcoord <- lattr$plcoord[subset]
      attributes(lx[,lj]) <- lattr
    }
  lx
}

## ========================================================================
plyx <-
  function(x=NULL, y=NULL, group=NULL, data=NULL, type="p", panel=plpanel,
           xlab=NULL, ylab=NULL, 
           markextremes=0, rescale=TRUE, mar=NULL,
           plargs = NULL, ploptions = NULL, ...)
{
  ## --- intro
  lcl <- match.call()
  lcall <- sys.call() ## match.call modifies argument names
  lcnm <- i.def(names(lcall), names(lcl))
  names(lcall) <- ifelse(lcnm=="", names(lcl), lcnm)
  if (length(plargs)==0) {
    lcall$markextremes <- markextremes
    ldtnm <- substitute(data)
    ldtnm <- if (is.name(ldtnm))
               as.character(ldtnm) else if (length(ldtnm)) format(ldtnm)
    lcall$.subdefault <- if (length(ldtnm)<30) ldtnm else ""
    lcall[1] <- list(quote(pl.control))
    lcall <- as.call(lcall)
    plargs <- eval(lcall, envir=parent.frame())
    if (length(ploptions)) plargs$ploptions <- ploptions
  }
## else  eval(lcall[["plargs"]], sys.frame(1))
  ploptions <- plargs$ploptions
  ##
  pldata <- plargs$pldata
  lnr <- NROW(pldata)
  lnobs <- median(apply(pldata,2, function(x) sum(is.finite(x)) ))
  lxnm <- attr(pldata,"xvar")
  if (is.null(lynm <- attr(pldata,"yvar"))) {
    lynm <- last(lxnm)
    lxnm <- last(lxnm, -1)
  }
  lx <- pldata[,lxnm, drop=FALSE]
  lnx <- ncol(lx)
  ly <- pldata[,lynm, drop=FALSE]
  ## why so complicated? I need this when  group  is active
  for (lj in lynm) {
    lyj <- ly[[lj]]
    if (length(lplc <- attr(lyj, "plcoord"))) {
      ly[,lj] <- transferAttributes(
        if(inherits(lyj, "Surv")) Surv(lplc, lyj[,2], type=attr(lyj, "type")) else lplc ,
        lyj )
    }
  }
  lny <- ncol(ly)
  ly1 <- ly1g <- ly[,1]
  lrgy1 <- i.def(attr(ly1, "innerrange"), attr(ly1, "plrange") )
  ## style elements
  lIsmooth <- i.getploption("smooth")
  lIfirst <- TRUE
  lplab <- pldata[["(plab)"]]
  lpch <- pldata[["(pch)"]]
  lIpch <- length(lpch)>0
  lpcol <- pldata[["(pcol)"]]
  lmar <- i.getplopt(mar)
  loma <- i.def(i.getploption("oma"), par("oma"), valuefalse=rep(0,4)) ## ??? 
  lmgp <- i.getploption("mgp")
  lpsep <- i.getploption("panelsep")
  lcexgrp <- par("cex")
  lyaxcol <- 1
  ## group
  lgroup <- pldata[["(group)"]]
  lIgrp <- length(lgroup)>0
  if(lIgrp)   {
  ##  lgrpname <- i.def(attr(lgroup, "varname"), "group")
    lgrplab <- as.character(unique(lgroup))
  } else  lgroup <- rep(1, nrow(ly))
  if (is.factor(lgroup)) lgroup <- as.numeric(lgroup)
  lgrp <- unique(lgroup)
  lngrp <- length(lgrp)
  ## ranges
  lrg <-
    sapply(ly, function(y)
      i.def(i.def(attr(y, "innerrange"), attr(y,"plrange")),
            range(y, na.rm=TRUE))
      )
  ## inner plotting range
  lnmody <- as.matrix(sapply(ly, function(x) c(attr(x,"nmod"),0,0)[1:2]>0 ))
  lIinner <- apply(lnmody, 1, any)
  if (lny>1) {
    if (!rescale) {
      lrgy1 <- c(min(lrg[1,]),max(lrg[2,]))
      ly <- apply(ly, 2, plcoord, range=lrgy1)
        ##for (lj in 1:lny)  ly[,lj] <- plcoord(ly[,lj], lrgy1)
      ly1 <- ly[,1]
    } else
    attr(ly1,"plrange") <-     ## extend
      lrgy1 + diff(lrgy1)*c(-1,1)*
      ifelse(lIinner, i.getploption("innerrange.ext"), 0) ## ploptions$plext
  }
  ## mark extremes
  lmark <- i.getplopt(markextremes)
  if (is.function(lmark)) lmark <- lmark(lnobs)
  lmk <- unlist(lmark)
  lImark <- length(lmk)>0 && any(ifelse(is.na(lmk),TRUE,lmk>0))
  lImark <- is.na(lImark)||lImark
  lymark <- if (lImark & lny==1) ly  ## cannot mark extremes if  lny>1
  ##
  if (lny>1 & is.null(mar)) lmar[4] <- lmar[2]
  ploptions$mar <- lmar
  lsmcol <- i.getploption("smoothline.col")
  ## --- multiple figures
  lnpgc <- lnpgr <- 1
  lnr <- lnx
  lnc <- lngrp
  loldp <- NULL
  if (lnx>1 & lngrp>1) {
    ploptions$mar <- lmar <- rep(lpsep,4) ## c(lmar[1],0.5,0.5,0.5)
    loma[1:2] <- pmax(lmgp[1]+1, loma[1:2])
    if (lny>1) loma[4] <- max(lmgp[1]+1, loma[4])
    if (lnx>1) {
      lmfig <- plmfg(lnx, lngrp, oma=loma)
      loldp <- lmfig$oldpar
      lnr <- lmfig$mfig[1]
      lnc <- lmfig$mfig[2]
      lnpgr <- ceiling(lnx/lnr)
      lnpgc <- ceiling(lngrp/lnc)
    }
  } else if (lngrp>1) {
    ploptions$mar <- lmar <- lpsep + c(0,0,1.5,0)
    loma[-3] <- pmax(lmgp[1]+1, loma[-3])
    loldp <- par(oma=loma)
  }
  ## --- plot
  for (ipgr in 1:lnpgr) {
    lr <- (ipgr-1)*lnr ## start at row index  lr
    for (lj in 1:min(lnr, lnx-lr)) {
    lxjg <- lxj <- lx[,lr+lj]
    if (lImark) {
      lplab <- plmark(lxj, y=lymark, markextremes=lmark, plabel=plargs$plabel)
      pldata[["(plab)"]] <- lplab
      plargs$pldata <- pldata
    }
    lpchg <- lpch
    ## groups
    for (ipgc in 1:lnpgc) { ## columns on page
      lc <- (ipgc-1)*lnc
      for (lig in (lc+1):lnc) { ## groups
        lg <- lgrp[lig]
        if(lny>1)  {
          lyaxcol <- lsmcol <- lpcol <- attr(ly1, "col")
          plargs$ploptions$smoothline.col <- lsmcol
          pldata["(pcol)"] <- lpcol
          if(!lIpch) 
            lpchg <-
              if (lny>1) attr(ly1, "pch") else i.getploption("pch")[1]
        }
        ## frame
        lop <- plframe(lxj, ly1, xlab=xlab, ylab=ylab,
                      axcol=c(NA,lyaxcol,NA,NA), ploptions=ploptions)
        if (lIfirst) {
          loldp <- c(lop, loldp)
          on.exit(par(loldp[unique(names(loldp))]))
        }
        lIfirst <- FALSE
        lrgold <- lrgy1
        lyg <- ly
        ## group
        if (lIgrp) {
          if (lmar[3]>=1 | par("mfg")[1]==1)
            pltitle(main=NULL, sub=lgrplab[lg], ##paste(lgrpname, lgrplab[lg], sep=": "),
                    outer.margin=FALSE, xpd=TRUE,
                    line=min(1.2*lcexgrp,par("mar")[3]-0.5*lcexgrp),
                    cex=rep(lcexgrp,2))  ## !!! sub!
          li <- which(lgroup==lgrp[lg])
          lxjg <- lxj[li]
          lyg <- plsubset(ly, li) ## transferAttributes(ly[li,], ly)
          ly1g <- lyg[,1] ## transferAttributes(lyg[,1], lyg)
          plargs$pldata <- pldata[li,]
          if (length(lpch)==lnr) lpchg <- lpch[li]
          if (length(lplab)==lnr) lplabg <- lplab[li]
        }
        if (lny>1) plargs$pldata$"(pcol)" <- attr(ly1g,"col")
        panel(lxjg, ly1g, type=type, plargs=plargs)
        ## multiple y
        lusr <- par("usr")
        if (lny>1) {
          for (lj in 2:lny) {
            lyjg <- lyg[,lj]
            lrgj <- lrg[,lj]
            lpcol <- attr(lyjg, "col") ##  the color must reflect the variable
            if (!lIpch) lpchg <- attr(lyjg, "pch")
            if (rescale) {
              lusr[3:4] <- lrgj[1] + diff(lrgj)/diff(lrgold)*(lusr[3:4]-lrgold[1])
              par(usr=lusr)
            }
            if (lIsmooth) {
              plargs$ploptions$smoothline.col <- lpcol
              plsmooth(lxjg, lyjg, plargs=plargs)
            }
            plpoints(lxjg, lyjg, type=type, plargs=plargs,
                     plab=lplab, pch=lpchg, col=lpcol)
            ##!!!linecolor
            lrgold <- lrgj
            if (rescale & lj==2) {
              lmfg <- par("mfg")
              plaxis(4, ly[,lj], lab=lmar[4]>=lmgp[2]+1 | lmfg[2]==lmfg[4],
                     range=lrgj, col=lpcol,
                     tickintervals=i.getploption("tickintervals"))
            }
          }
        }
      }
    }
  }
  }
  pltitle(plargs=plargs, show=0.5)
  stamp(sure=FALSE, ploptions=ploptions)
}
## ==========================================================================
varattributes <-
  function(data, attributes = NULL)
{
  data <- as.data.frame(data)
  if (!is.list(attributes))
    stop("!varattributes! argument 'attributes' must be a list")
  if (is.null(names(attributes))) {
    if (length(attributes)!=ncol(data))
      stop("!varattributes! argument 'attributes' must have names ",
           "or be of appropriate length")
    names(attributes) <- names(data)
  }
  lnames <- names(attributes)
  if (any(linm <- lnames%nin%names(data)))
    stop("!varattributes! names of argument 'attributes' not in ",
         "names of 'data': ", paste(lnames[linm], collapse=", "))
  for (lnm in lnames) {
    lattr <- attributes(data[[lnm]])
    lattr[names(attributes[[lnm]])] <- attributes[[lnm]]
    attributes(data[[lnm]]) <- lattr
    if (lnm=="innerrange") {
    ## call plcoord !!!
    }
  }
  invisible(data)
}
## ==========================================================================
plotregr.control <-
  function(x, data = NULL, xvar = TRUE, transformed = FALSE,
           ## generate variables used for plotting
           weights = NULL, 
           ## specify some specific aspects / contents of plots
           glm.restype = "working", condquant = TRUE, smresid = TRUE, 
           partial.resid = TRUE, cookdistlines = 1:2,
           leveragelim = c(0.99, 0.5), condprobRange = NULL,
           ## smooth and refline
           refline = TRUE,
           reflineband = FALSE, testlevel = 0.05,
           smooth = 2, ## smoothPar=NA, smoothIter=NULL,
           smooth.sim=NULL,
           xlabs = NULL, reslabs = NULL, markextremes = NULL, 
           ## multiple frames
           mf = TRUE, mfcol = FALSE, multnrow = 0, multncol = 0, multmar = NULL,
           oma = NULL, ... )
{ ## get data for plotting, collect and check arguments
  ## do preparations that are common to plot.regr_ and plresx
  ## --------------------------------------------------------------
  lcall <- match.call()
  lnaaction <- x$na.action
  if (!is.null(lnaaction)) class(lnaaction) <- "exclude"
  x$na.action <- lnaaction
  ## --- family
  lfam <- c(x$distrname, x$family$family, "")[1]
  if (lfam=="" & inherits(x, "lm")) lfam <- "gaussian"
  lfamgauss <- lfam%in%c("gaussian","Gaussian")
  lfamcount <- (lfam%in%c("binomial","poisson")&&length(table(x$y)<=2)) | 
    inherits(x,"polr") ## lfam=="multinomial" | 
  ## --- na.action: always get full data
  ## residuals first because they fix the number of observations
  lres <- NULL
  lcq <- i.getplopt(condquant)
  rtype <- i.def(i.def(if (lfamcount & lcq)
                         rtype <- "condquant", glm.restype), "working")
  if (inherits(x, "survreg"))
    lres <- residuals.regrsurvreg(x, type = if (lcq) "condquant")
  if (inherits(x, "coxph"))
    lres <- residuals.regrcoxph(x, type = if (lcq) "condquant")
  if (inherits(x, "polr"))
    x$residuals <- lres <- residuals.regrpolr(x, type = if (lcq) "condquant")
  if (inherits(x, "lm"))
    lres <-
      if (lfamcount)  residuals.regrpolr(x, type=rtype)
      else  residuals(x, type = rtype)
  if (inherits(x, "regrMer")) lres <- naresid(lnaaction, resid(x))
  ## if (length(lres)==0) lres <- residuals(x)
  lres <- as.data.frame(lres)
  lmres <- ncol(lres)
  lnobs <- sum(is.finite(lres[,1]))
  lres0 <- all( apply(lres[,1:lmres, drop=FALSE],2,
                      function(x) all(x==notna(x)[1], na.rm=TRUE ) ) )
  if (lres0)
    stop("!plot.regr/plresx! all residuals are equal -> no residual plots")
  lnres <- nrow(lres)
  ## --- ldfres 
  ldfres <- df.residual(x)
  ldfmod <- i.def(x$rank, length(coef(x)))
  if (is.null(ldfres))  ldfres <- lnres-ldfmod
  ## --- sigma
  lsigma <- x$sigma
  if (length(lsigma)==0) lsigma <- c(x$scale, summary(x)$sigma)[1]
  if (length(lsigma)==0)
    lsigma <- if (lfamcount) 0 else sqrt(apply(lres^2,2,sum, na.rm=TRUE)/ldfres)
  x$sigma <- lsigma
  ## --- standardized residuals
  llev <- x$leverage
  lstres <- attr(lres, "stresiduals")
  if (length(lstres)==0) {
    llevlim <- i.getplopt(leveragelim)
    lrs <- i.stres(x, residuals=lres, leveragelim=llevlim)
    attributes(lres) <- c(attributes(lres), lrs)
    lstres <- attr(lres, "stresiduals")
    if (is.null(llev)) llev <- lrs$leverage
  } ## else lstres <- lres
  else llev <- naresid(lnaaction, llev)
  ##  if (class(lnaaction)=="exclude") lstres <- i.naresid.exclude(lnaaction, lstres)
  ## --- xvar
  lform <- if (inherits(x, "regrMer")) x$formula else formula(x)  ## formula(x) inappropriate for merMod
    ## ... if it inherits from  lm
  lmodvdupl <- u.varsin2terms(lform[-2])
  if (u.notfalse(xvar)) {
    lxf <- if (transformed) lform else u.asformula(all.vars(lform[-2]))
    if (!(is.null(xvar)|u.true(xvar))) {
      if (is.character(xvar)) {
        lxvarf <- u.asformula(setdiff(xvar,"."))
        lxf <- if ("."%in%xvar) update(lxf, lxvarf) else lxvarf
      }
      else {
        if (is.formula(xvar)) lxf <- update(lxf, xvar)
        else
          stop("!plotregr.control! Inadequate argument 'xvar'")
      }
    }
    lxvar <- getvarnames(lxf, transformed=TRUE)$xvar
    lxvraw <- u.allvars(lxvar)
  } else lxvar <- NULL
  ## residual names
  lyexpr <- deparse(lform[[2]])
  lynm <- if (nchar(lyexpr)>10) "Y" else lyexpr ## 
  lresname <- paste("res_", if (lmres>1) colnames(lres) else lynm, sep="")
  names(lres) <- lresname
  ## --- prepare  pl.control
  lcall <- as.list(match.call())[-1]
  ladrop <- c("xvar", "glm.restype", "smresid", "partial.resid",
              "cookdistlines", "leveragelim", "smooth", "smooth.sim",
              "refline", "reflineband", "testlevel", "xlabs", "reslabs",
              "mf", "mfcol", "multnrow", "multncol", "multmar", "oma")
  lcall$y <- lres
  lcall <- c(as.list(quote(pl.control)),
             as.list(lcall[setdiff(names(lcall),ladrop)]))
  if (length(lxvar)&u.notfalse(lxvar)) {
    lcall$x <- lxvar
    lcall$transformed <- transformed
    lcall$.subdefault <- i.form2char(lform)
    ## --- data argument
    if (length(lxvar)||any(names(lcall)%in%i.argPldata)) {
      ldata <- x$model
      if (length(ldata)&&length(lnaaction))
        ldata <- i.naresid.exclude(lnaaction, ldata)
      if (length(ldata)==0||!(transformed & all(lxvar%in%names(ldata)))) {
        ldata <- data
        if (length(ldata)==0) {
          if (length(x$allvars))
            ldata <- i.naresid.exclude(lnaaction, x$allvars)
          else ldata <- eval(x$call$data)
          ##-     } else {
##-       if (length(lav <- x$allvars)&&NROW(data)==NROW(lav)) 
##-        ldata <- cbind(data, lav[colnames(lav)%nin%colnames(data)])
##-       else stop("!plotregr.control! incompatible data and x$allvars")
    ##-     }
        }
        if (length(ldata)==0) {
          ldata <- i.naresid.exclude(lnaaction, x$model)
          if (any(lxvj <- lxvar%nin%names(ldata)))
            stop("!plot.regr/plresx! variable(s) ",
                 paste(lxvar[lxvj], collapse=", "), " not found")
        }
        if (length(ldata)==0)  stop("!plotregr.control! No data found")
      }
  ##
      if (lnres!=nrow(ldata)) {
        if (class(lnaaction)=="omit") ldata <- ldata[-lnaaction,]
        ##    if (lnr!=nrow(ldata)) ldata <- x$model ## needs at least a warning!
        if (lnres!=nrow(ldata))
          stop("!plotregr.control! nrow of residuals and data do not agree.")
      }
    }
  } else {
    ldata <- lres ## needed to get number of observations
    lcall$x <- FALSE
  } ## needed to get  nobs  in pl.control
  lcall$data <- ldata
  lftext <- i.form2char(lform)
  lcall$.subdefault <- lftext
  ## margins for multivariate regression
  if (lmres>1) {
    lmar1 <- i.def(i.getploption("mar")[1],3)
    lmar2 <- rep(i.getploption("panelsep"), length=4)[2:4]
    lcall$mar <- c(lmar1,lmar2)
  }
  mode(lcall) <- "call"
  ## -------------------------------------
  plargs <- eval(lcall, parent.frame())  ## pl.control
  ## -------------------------------------
  lpldata <- plargs$pldata
  xvar <- attr(lpldata, "xvar")
  ## attributes for residuals
  lres <- genvarattributes(lres, varlabels = lresname, ploptions=ploptions)
  if (lmres>1) { ## multivariate
    if (is.null(lcn <- colnames(lres))) lcn <- 1:ncol(lres)
    colnames(lres) <- lcn # paste("res", lcn, sep=".")
  }
  ## mark extreme  stres
  lmxdef <- markextremes(lnobs)
  if (is.atomic(markextremes)) {
    markextremes <- i.def(markextremes, NA)
    if (anyNA(markextremes)) markextremes <- lmxdef
    markextremes <- 
      list("(res)"=markextremes, "(fit)"=0, "(lev)"=c(0,max(markextremes)) )
  }
  lmxres <- i.def(markextremes$"(res)", lmxdef)
  lresplab <-
    if (lmxres>0)
      apply(lstres, 2,
            function(x) plmark(lmxres, x, plabel=plargs$plabel) )
    else NULL
## -------------------------------------------
  ## --- prepare objects needed for plotting
  ## --- smoothWeights, used for smooth calculation: get from  x  if needed
  lsmwgt <- lpldata[["(smoothWeights)"]]  ## possibly only logical
  lIsmweights <- is.logical(lsmwgt)&&all(lsmwgt) ## weights explicitly required
  lInosmweights <- is.logical(lsmwgt)&&!any(lsmwgt) ## weights explicitly denied
  if (lIsmweights | (length(lsmwgt)&&all(is.na(lsmwgt))))
      lsmwgt <- naresid(lnaaction, x$weights)
  lIsmwgt <-
    length(lsmwgt)>1 && any(lsmwgt!=notna(lsmwgt)[1],na.rm=TRUE) 
  if (lIsmweights&!lIsmwgt)
    warning(":plot.regr/plresx: no weights found for smooth calculation.")
  lpldata[["(smoothWeights)"]] <- lsmweights <-
    if (lIsmwgt)  lsmwgt / mean(lsmwgt, na.rm=TRUE) else NULL
  ## --- psize, used as sizes of plotting characters: same as weights
  lpsize <- lpldata[["(psize)"]]  ## possibly only logical
  lIpsize <- is.logical(lpsize)&&all(lpsize) ## psize expl. required
  lInopsize <- is.logical(lpsize)&&!any(lpsize) ## psize expl. denied
  if (is.null(lpsize) | lIpsize | (length(lpsize)&&all(is.na(lpsize))))
      lpsize <- if (length(lsmwgt)) lsmwgt else naresid(lnaaction, x$weights)
  lIpsz <-
    length(lpsize)>1 && any(lpsize!=notna(lpsize)[1],na.rm=TRUE)
  if (lIpsize&!lIpsz)
    warning(":plot.regr/plresx: no plot sizes found.")
  lpldata[["(psize)"]] <- 
    if (lIpsz)  lpsize / mean(lpsize, na.rm=TRUE) else NULL
  ## mahalanobis residuals
  lresmahal <- naresid(lnaaction, x$resmahal)
  if (lmres>1) 
    if (is.null(lresmahal)) lresmahal <- mahalanobis(lres,0,var(lres, na.rm=TRUE))
  ## --- simulated residuals
  ## when using  smooth.group , default is  0
  lnsims <- i.def(smooth.sim, 19, 19, 0)
  if (inherits(x, c("mlm", "polr", "survreg", "coxph", "regrMer"))) lnsims <- 0
  if (lmres>1) lnsims <- 0 # !!! not yet programmed for mlm
  if (lnsims>0 & !inherits(x, c("lm","glm","nls"))) {
    warning(":plot.regr/simresiduals: ",
            "I can simulate only for 'lm', 'nls' and 'glm' objects")
    lnsims <- 0
  }
  lsimres <- NULL
  if (lnsims>0) {
    lsimres <- if(u.debug())
               simresiduals(x, lnsims, glm.restype=glm.restype)  else
      try(simresiduals(x, lnsims, glm.restype=glm.restype), silent=TRUE)
    if (class(lsimres)=="try-error") {
      warning(":plot.regr/simresiduals: simresiduals did not work. ",
              "No simulated smooths")
      lsimres <- NULL
      lnsims <- 0
    }
  }
  ## --- some more arguments
##  xlabs <- i.def(xlabs, NULL, NULL, NULL) --> genvarattr
  reslabs <- i.def(reslabs, NULL, NULL, NULL)
  ## --- multiple frames
  mf <- i.def(mf, NULL, TRUE, NULL)
  oma <- i.def(oma, NULL, valuefalse=0)
  if (length(oma)==2) oma <- c(0,0,oma)
  ## --- more arguments
  smresid <- i.def(smresid, TRUE)
  if (lmres>1) smresid <- FALSE ## !!! muss noch gemacht werden
  refline <- i.def(refline, TRUE)
  reflineband <- i.def(reflineband, FALSE, TRUE, FALSE)
  testlevel <- i.getplopt(testlevel)
  if (testlevel<=0 | testlevel>=1)
    stop("!plotregr.control! invalid test level")
  partial.resid <- i.def(partial.resid, TRUE)
  cookdistlines <- i.def(1:2, 1:2, NULL)
  ## ------------------------------------------------------------
  ## result of plotregr.control
  plargs$pldata <- lpldata
  plargs$smooth <- i.def(smooth, 2, 2, 0)
  rr <- c(plargs,
    list(
      residuals = lres, rescol = lmres, 
      ## weights = x$weights,
      leverage = llev,
      ## linear.predictor = x$linear.predictor,
      resmahal = x$resmahal, 
      ## resmahal = lresmahal,
      simres = lsimres, ## simstres <- lsimstres,
      yexpr = lyexpr, resname = lresname, ## stresname = lstresname,
      ##    absresname = labsresname, fitname = lfitname,
      family = lfam, famgauss = lfamgauss, famcount = lfamcount,
      formula = lform, na.action = lnaaction, 
      sigma=lsigma, df.residual = ldfres, df.model = ldfmod, 
      ## -- return arguments
      glm.restype = glm.restype, smresid = smresid,
      partial.resid = partial.resid, cookdistlines = cookdistlines,
      leveragelim = leveragelim, ## condprobRange = condprobRange,
      ## smooth and refline
      smooth.sim = lnsims,
      refline = refline,
      reflineband = reflineband, testlevel=testlevel,
      resplab = lresplab,
      mf=mf, multnrow = multnrow, multncol = multncol, multmar = multmar, oma=oma #,
    ) )
  assign(".plargs", rr, pos=1)
  rr
}
## -----------------------------------------------------------------------
i.merprep <- function(x) {
  ldfr <- df.residual(x)
  rr <- list(
    ## na.aaction = NULL,
    family = family(x),
    fitfun = "mer",
    coefficients = coef(x),
    fitted = fitted(x),
    sigma = sqrt(sum(residuals(x)^2, na.rm=TRUE) / ldfr),
    ## leverage = NULL,
    model = x@frame,
    residuals = residuals(x), 
    call = x@call,
    formula = formula(x),
    na.action = if (length(lnaaction <- attr(x@frame, "na.action")))
                  structure(lnaaction, class="exclude")
    )
  class(rr) <-
    c("regrMer",
      switch(rr$family$family, gaussian="lm", binomial="glm", poisson="glm",
             "lm") )
  rr
}

## -----------------------------------------------------------------------
i.argPldata <- c("psize", "plab", "pch", "pcol",
                 "group", "smooth.group", "smooth.weights")
i.argPlcontr <-
  c("x", "y", "data", "transformed", "subset", ## "cex", "markextremes",
    "ycol", "ylty", "ypch", ##  "smooth",
    "main", "sub", ".subdefault", ## "mar",
    "xlab", "ylab", "varlabels",
    "ploptions", ".environment.")

i.argPlControl <- ##!!!
  c("x", "y", "data", "xvar", "transformed", i.argPldata,
    "refline", "ploptions", 
    "main", "sub", "cex.main", "varlabels",
    ## plotregr.control
    "xvar", "weights", "glm.restype", "smresid",
    "partial.resid", "cookdistlines", "leveragelim", "condprobRange", 
    "reflineband", "testlevel", "smooth.sim",
    "xlabs", "reslabs", 
    "mf", "mfcol", "multnrow", "multncol", "multmar", "oma"
    )
i.argPlregr <- c("plotselect", "sequence", "addcomp", "smooth.legend")

## ====================================================================
plot.regr <-
  function(x, data=NULL, plotselect = NULL, xvar = TRUE,
           transformed = FALSE, sequence=FALSE, weights=NULL, 
           addcomp = FALSE, smooth.legend = FALSE,
           markextremes = NA, ...)
{
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date:  7 May 93 / 2002
  argPlregr <- c("x", "data", "plotselect", ## "xvar", "transformed",
                 "sequence", "weights", "addcomp", "smooth.legend",
                 "reflinesband")
  lImer <- inherits(x, "merMod")
  if (lImer) x <- i.merprep(x)
  ## ----------------
  lcl <- match.call()
  lcall <- sys.call() ## match.call modifies argument names
  lcnm <- i.def(names(lcall), names(lcl))
  names(lcall) <- ifelse(lcnm=="", names(lcl), lcnm)
  lcall <- lcall[setdiff(names(lcall), i.argPlregr)]
  lcall$x <- x
  lcall[1] <- list(quote(plotregr.control)) 
##            lac[intersect(names(lac), c(i.argPlControl, names(.ploptions)))])
##            list(fit=TRUE, hat=TRUE, stresiduals=TRUE)  )
  mode(lcall) <- "call"
  plargs <- eval(lcall, parent.frame())
  ## -------------------------------------------------------------------
  lres <-  plargs$residuals
  lmres <- plargs$rescol
  lpldata <- plargs$pldata
  ploptions <- plargs$ploptions
  lsmgrp <- lpldata$"(smooth.group)"
  lsmooth <- i.getploption("smooth")
  lsmpar <- i.getploption("smooth.par")
  lsigma <- plargs$sigma
  llev <- plargs$leverage
  lresname <- plargs$resname
  lsimres <- plargs$simres
  lnsims <- if (length(lsimres)==0) 0 else ncol(lsimres)
  llevlim <- plargs$leveragelim
  lrefline <- i.def(plargs$refline, TRUE)
  ## lsimstres <- plargs$simstres
  x$na.action <- lnaaction <- plargs$na.action
  ## from x
  lform <- plargs$formula ## formula(x)
  lweights <- naresid(lnaaction, plargs$weights)
  lIwgt <- length(lweights)>0
  ## number of observations
  lnr <- nrow(lres)
  lnna <- apply(lres, 1, function(x) all(is.finite(x)))
  lnobs <- sum(lnna)
  ## smooth legend
  if (length(names(smooth.legend))==0) {
    lsmlegend <- i.def(smooth.legend, NULL, TRUE, NULL) 
    if (length(lsmlegend)==1)
      lsmlegend <-
        setNames(rep(lsmlegend,5),
                 c("yfit","resfit","absresfit","absresweight","(xvar)") )
  } else lsmlegend <- smooth.legend
  ## -------------------------
  ## family
  if (inherits(x,"mulltinom"))
    stop("!plot.regr! I do not know how to plot residuals of a multinomial regression")
  lglm <- inherits(x, "glm")
  lbinary <- lglm && length(unique(plargs$y))==2 ## binary binomial
  lcensored <- inherits(plargs$y, "Surv")
  lnnls <- !inherits(x, "nls")
  if (u.true(plargs$famcount)) plargs$ploptions$smooth.iter <- 0
  condquant <- NULL
  lIcq <- u.true(i.getploption("condquant")) & (lbinary|lcensored)
  ## plot selection
  lplsel <- unlist(
    i.plotselect(plotselect, smooth=plargs$smooth, Iwgt=lIwgt,
                 mult=lmres>1,
                 famgauss=plargs$famgauss, famglm=inherits(x,"glm"),
                 famcount=plargs$famcount)
    )
  ## -----------------------------------
  lplsel <- lplsel[is.na(lplsel)|lplsel>0]
  lnplsel <- sum(names(lplsel)%in%
                c("yfit","resfit","absresfit","qq","absresweight") )
  ## ---------------------------
  ## simulated residuals
##-   lsimres <- NULL
##-   lnsims <- plargs$smooth.sim
##-   if (lnsims>0) 
  ##-     lsimres <- simresiduals(x, lnsims, glm.restype=glm.restype)
  ## -----------------
  ## fit
  lfit <- x$linear.predictor
  lfitname <- "linear predictor"
  if (inherits(x, "polr")) lfit <- fitted.regrpolr(x, type="link")
  if(is.null(lfit)) {
    lfit <- fitted(x)
    lfitname <- "fitted value"
  } else lfit <- naresid(lnaaction, lfit)
  lfitname <- rep(lfitname, length=lmres)
  lir <- i.getploption("innerrange.fit") ## extra element of ploptions!
  lfit <- as.data.frame(lfit)
  lfit <- genvarattributes(as.data.frame(lfit), varlabels = lfitname,
                           innerrange.limits=lir)
  ## standardized residuals
  lstres <- as.data.frame(attr(lres, "stresiduals"))
  lstrratio <- as.data.frame(attr(lres, "strratio"))
  ## llev <- attr(lres, "leverage")
##-   if (length(lstres)==0) {
##-     lrs <- i.stres(x, leveragelim=llevlim)
##-     lstres <- lrs$stresiduals
##-     llev <- lrs$leverage
##-     lstrratio <- lrs$strratio
##-     {
##-     if (!plargs$famcount)
##-       warning(":plot.regr: I do not find standardized residuals.",
##-               " I will use raw residuals instead.")
##-     lstres <- lres
##-     lstrratio <- rep(1, lnr)
##-   } else {
##-     lstrratio <- i.def(attr(lres, "strratio"), rep(1,lnr))
##-   }
##-   lstrratio <- as.matrix(lstrratio)
##-   if (is.null(llev)) llev <- leverage(x)
  lresplab <- plargs$resplab
  lIrpl <- length(lresplab)>0
  if (any(c("absresfit","qq","leverage")%in%names(lplsel))) {
  ##  lstres <- as.data.frame(as.matrix(lres)*lstrratio)
    lstresname <- paste("st.", lresname, sep = "")
    labsresname <- paste("|st.",lresname,"|", sep="")
    ##
    for (lj in seq_len(lmres)) {
      ## residuals from smooth
      if (plargs$smresid) {
        lrsj <- lres[,lj]
        lfsm <- gensmooth(lfit[,lj], lrsj, plargs=plargs, ploptions=ploptions) 
        lfsmr <- residuals(lfsm)
        if (lna <- sum(is.na(lfsmr)&!is.na(lres[,lj])))
          warning(":plot.regr: residuals from smooth have ",
                  round(100*lna/lnobs,1), " % additional NAs")
        lstres[,lj] <- lstrsj <-
          lfsmr * lstrratio[,lj] ## !!! * f(leverage of smooth)
      }
##-       if (length(lcqj <- attr(lres[,lj], "condquant")))
##-         attr(lstres[,lj], "condquant") <- lcqj[,"index",drop=FALSE] 
      if (lnsims) lsimstres <- lsimres * lstrratio[,lj] ## index needed for multiv
    }
    if (plargs$smresid) {
      lstresname <- paste("st.sm.", lresname, sep = "")
      labsresname <- paste("|st.sm.",lresname,"|", sep="")
    }
    names(lstres) <- lstresname
    lstres <- genvarattributes(lstres, ploptions=ploptions)
  }
  ## --- multiple frames xxx
  lmf <- plargs$mf  ## i.def(lmf, TRUE, TRUE, FALSE)
  if (length(lmf)) {
    if (u.true(lmf)) ## is.logical(lmf)&&lmf)
      lmf <- if (lmres>1) {
               if (lmres<=4) c(NA, lmres) else lmres
             } else {
               if (lnplsel<=2) c(1,lnplsel) else c(NA, 2)
             }
  }
  if (length(lmf)==2 && is.na(lmf[1])) {
    lmf1 <- lnplsel
    if(lmf1>lmf[2]+1) lmf1 <- lmf[2]
    lmf[1] <- lmf1
  }
  loma <- i.def(plargs$oma, c(2,1)*(length(lmf)>0), valuefalse=NULL)
  if (length(loma)<4) loma <- c(0,0,loma,0)[1:4]
  loldpar <- 
    if (length(lmf)&(!is.logical(lmf))) {
      lop <- if (length(lmf)==1) attr(plmfg(mft=lmf, oma=loma),"oldpar") else
      attr(plmfg(lmf[1], lmf[2], oma=loma),"oldpar")
      lop[setdiff(names(lop), c("mfig","mrow","mcol"))]
    } else par(oma=loma) ## , ask=plargs$ask
  on.exit(par(loldpar), add=TRUE)
  lnewplot <- TRUE ## !!!
  ## --------------------------------------------------------------------------
  ## start plots
  if (length(lplsel))
    for (liplot in 1:length(lplsel)) {
      lpllevel <- lplsel[liplot]
      lpls <- names(lpllevel)
##  y on fit
  if(lpls=="yfit") {
    if (is.na(lpllevel)) lpllevel <- 3*(lplsel["resfit"]==0)
    if (is.na(lpllevel)) lpllevel <- 0
    if (lpllevel>0) { ## !!! condquant bereinigen 
      lsimy <- if(length(lsimres)) lfit + lsimres  else NULL ## !!! mult!
      lsml <- if (length(lsmlegend["yfit"]))
                setNames(rep(lsmlegend["yfit"],length(plargs$yname)),
                         plargs$yname)  else lsmlegend
      lyvar <- attr(lpldata, "yvar")
      ly <- lpldata[,lyvar]
      for (lj in seq_len(lmres)) {
        plframe(ly[,lj],lfit[,lj], ploptions=ploptions)
        if (lpllevel>1)
          plsmooth(lfit[,lj], cbind(ly[,lj], lfit[,lj]+lsimres),
                   plargs=plargs, band=lpllevel>2)
        plpoints(lfit[,lj], ly[,lj], plargs=plargs,
                 plab=if(lIrpl) lresplab[,lj])
        pltitle(plargs=plargs, show=NA)
      }
    }
  }
  ## ---
  if(lpls=="resfit") {
    for (lj in seq_len(lmres)) {
      plframe(lfit[,lj], lres[,lj], ploptions=ploptions)
      plargs$smooth <- lpllevel-1
##-      plpanel(lfit[,lj], lres[,lj], plargs=plargs, frame=TRUE, title=NA)
      if (lpllevel>1)
        plsmooth(lfit[,lj], cbind(lres[,lj], lsimres),
                 plargs=plargs, band=lpllevel>2)
      if (lrefline)
        plrefline(c(x=median(lfit[,lj], na.rm=TRUE),y=-1),
                   x=lfit[,lj], y=lres[,lj], plargs=plargs)
      plpoints(lfit[,lj], lres[,lj], plargs=plargs,
               plab=if(lIrpl) lresplab[,lj])
      pltitle(plargs=plargs, show=NA)
##      plargs$ploptions$refline <- NULL
    }
  }
## ---
  if(lpls=="absresfit")
    if(length(lstres)) {
      ploptsmod <- ploptions
      lplext <- rep(ploptions("plext"), length=4)
      lplext[3] <- 0
      ploptsmod$plext <- lplext
      lir <- i.def(i.getploption("innerrange"))
      for (lj in seq_len(lmres)) {
        lfitj <- lfit[,lj]
        lstrj <- lstres[,lj]
        lattrj <- attributes(lstres[,lj])
        if (lir && length(lirgj <- lattrj$innerrange)) {
          lirgj <- c(0, max(abs(as.numeric(lirgj))) )
          labssrj <-
            genvarattributes(as.data.frame(abs(lstrj)), innerrange.limits=lirgj)[,1]
        } else {
   ##       lirgj <- c(0, max(abs(lattrj$plrange)))
          labssrj <- genvarattributes(as.data.frame(abs(lstrj)))[,1]
   ##       attr(labssrj, "plrange") <- lirgj
        }
        attr(labssrj, "plrange")[1] <- 0
        attr(labssrj, "condquant") <- attr(lstrj,"condquant")
        ##
        plframe(lfitj, labssrj, ylab=labsresname[lj], ploptions=ploptsmod)
        if (lpllevel>1) {
          lsimabsres <- if (lnsims) abs(lsimstres) else NULL
          plsmooth(lfitj, cbind(as.matrix(labssrj), lsimabsres), power=0.5,
                   plargs=plargs, band=lpllevel>2)
          }
        plpoints(lfitj, labssrj, condquant=0, ## plargs=plargsmod, 
               plab=if(lIrpl) lresplab[,lj], plargs=plargs) 
        pltitle(plargs=plargs, show=FALSE)
      }
    }  else
    warning(":plot.regr: No standardized residuals found")
## --- plot abs. res vs. weights
  if(lpls=="absresweights") {
    if (length(lweights)!=lnr)
      warning(":plot.regr: no suitable weights found.",
              "cannot plot absres on weights")
    else {
      plargsmod <- plargs
      lplext <- rep(ploptions("plext"), length=4)
      lplext[3] <- 0
      plargsmod$ploptions$plext <- lplext
      lwg <- lweights
      lwg[lwg<=0] <- NA
      labsres <-
        if (length(lstres)) {
          lrlab <- labsresname
          abs(lstres)
        } else {
          lrlab <- paste("|",lresname,"| * sqrt(w)", sep="")
          abs(lres)*sqrt(lwg)
        }
      for (lj in seq_len(lmres)) {
        labssrj <- labsres[,lj]
        lsm <-
          if (lpllevel>=1)
            gensmooth(lwg, cbind(labssrj, lsimabsres), band=lsmooth>2,
                      plargs=plargs, ploptions=ploptions)  else  NULL
        plframe(lwg, labssrj, ylab=lrlab[lj], plargs=plargsmod)
        plsmoothline(lsm, plargs=plargs)
        plpoints(lwg, labssrj, plargs=plargsmod, condquant=0,
                 plab=if(lIrpl) lresplab[,lj])
        pltitle(plargs=plargs, show=FALSE)
  } } }
## --- normal plot qq plot
  if(lpls=="qq") {
    ##    if (lpllevel>0)
    lnsims <- plargs$smooth.sim
    for (lj in seq_len(lmres)) {
      llr <- lstres[,lj]
      lIcqj <- length(lcq <- attr(llr, "condquant"))>0
      lIcq <- lIcq | lIcqj
      lpch <-
        if (lIcqj) i.pchcens(plargs, lcq) else ploptions$pch[1]
      lxy <- qqnorm(llr, ylab = lstresname[lj], main="", type="n")
      lquart <- quantile(llr, c(0.25,0.75), na.rm=TRUE)
      plrefline(c(0, diff(lquart)/(2*qnorm(0.75))), plargs=plargs)
      if (lnsims>0) {
        lxx <- qnorm(ppoints(lnobs))
        for (lr in 1:lnsims) {
          llty <- last(ploptions$smoothline.lty)
          lines(lxx,sort(lsimstres[,lr]), lty=llty,
                lwd=ploptions$linewidth[llty],
                col=last(ploptions$smoothline.col))
        }
      }
      ## qq line
      li <- order(lxy$x)
      lxx <- lxy$x[li]
      lyy <- lxy$y[li]
##      lines(lxx,lyy, col=plargs$ploptions$col[1])
      lpa <- plargs
      lpa$pldata <- lpldata[li,]
      plpoints(lxx, lyy, pch=lpch, plargs=lpa)  ## ??? simplify?
      if(lIcq & lj==lmres)
        legend("bottomright",
               pch=c(rep(ploptions$censored.pch,length=2)),
               legend=c("uncensored","censored"))
      pltitle(plargs=plargs, show=FALSE)
    }
  }
## --- leverage plot. If weight are present, use "almost unweighted h"
  if(lpls=="leverage")
  if ((!is.na(lpllevel))&&lpllevel>0 && lnnls) {
    if (lIwgt) llev <- llev/lweights
    if (diff(range(llev,na.rm=TRUE))<0.001)
      warning(":plot.regr: all leverage elements equal, no leverage plot")
    else {
      llevpl <-
        genvarattributes(data.frame(leverages=plcoord(llev, c(0,0.5),
                                      ploptions=ploptions)) )[,1]
      lstres <- genvarattributes(lstres, varlabels = lstresname)
      ## mark extremes
      lmx <- ploptions$markextremes
      if (is.list(lmx)) lmx <- lmx[["(lev)"]]
      if (is.function(lmx)) lmx <- lmx(lnobs)
      lmx <- last(i.def(lmx, NA))
      if (lImxlev <- is.na(lmx)||lmx>0)
        lpllev <- plmark(llev, markextremes=lmx, plabel=plargs$plabel)
      ##
      llevtit <- paste("leverages", if(lIwgt) "(de-weighted)")
      ldfmod <- plargs$df.model
      lcookl <- plargs$cookdistlines
      if (i.def(ldfmod, 0)<=1) {
        warning(":plot.regr: model degrees of freedom <=1. No Cook-distance lines")
        lcookl <- NULL
      }
      if (lIcook <- length(lcookl)>0) {
        llx <- seq(0, ## min(c(lcookl,2),na.rm=TRUE)^2 / (4*(ldfmod-1)), ## r_i^*<=4
                   max(llev,na.rm=TRUE), length=50)
        ## see formula for curves of constant Cook's distance in terms of
        ##   standardized residuals
        llrcd <- outer(sqrt((1-llx)/((ldfmod-1)*llx)), c(lcookl,-lcookl)) 
      }
      for (lj in 1:lmres) {
        lstrj <- lstres[,lj]
        plframe(llevpl, lstrj, xlab=llevtit, plargs=plargs)
        if (lIcook) plrefline(list(x=llx, y=llrcd), x=llevpl, plargs=plargs)
        lplj <- if (lIrpl)  lresplab[,lj]  else  rep("", lnr)
        if (lImxlev) lplj <- ifelse(lpllev=="", lplj, lpllev)
        plpoints(llevpl, lstrj, psize=if(lIwgt) lweights, ## condquant=0,
                 plab=if(lIrpl|lImxlev) lplj, plargs=plargs)
        pltitle(plargs=plargs, show=FALSE)
      }
    }
  }
  ## -----------------------------------------------------------------
  ## multivariate:
  ## residual matrix for multivariate regr
  if(lpls=="resmatrix") {
    lpanel <- function(xx, yy, indx, indy, pch, col, plab, plargs=plargs, ...) {
      plpoints(xx,yy, plargs=plargs)
    }
    if (plargs$rescol>1) {
      lpa <- plargs
      lpa$pldata <- lpa$residuals
      lpa$mar <- NULL
      plmatrix(lres, panel=lpanel, plargs=lpa) #main=plargs$main, pch=plargs$pch
    }
  }
  ## mahalanobis residuals
  if(lpls=="qqmult")  ## qq plot of Mahalanobis lenghts for multivariate regr
    if ((!is.na(lpllevel))&&lpllevel>0) {
      lresmahal <- plargs$resmahal
      if (is.null(lresmahal)) lresmahal <- mahalanobis(lres,0,var(lres,na.rm=TRUE))
      lxx <- sqrt(qchisq(ppoints(lresmahal),ncol(lres)))
      lor <- order(lresmahal)
      lyy <- sqrt(lresmahal[lor])
      lop <- par(mfrow=c(1,1), oma=c(0,0,2,0))
      plframe(lxx,lyy, xlab="sqrt(Chisq.quantiles)",
              ylab = "Mahal.oulyingness", ploptions=ploptions)
      lines(lxx,lyy)
      ## !!! needs work!!!
      ##      if (ltxt) text(lxx,lyy, plab[lor]) # else points(lxx,lyy,pch=lplab[lor])
      abline(0,1,lty = ploptions$grid.lty, col=ploptions$grid.col)
      pltitle(plargs=plargs, show=FALSE)
      stamp(sure=FALSE, ploptions=ploptions)
      par(lop)
    }
} ## end lplsel
## ----------------------------------------------------------------
  ## plot residuals vs. explanatory variables by calling plresx
##-   ##  --- sequence
##-   lIseq <- i.def(sequence, FALSE, TRUE, FALSE)
##-   ## lIseq <- lIseq|length(sequence)>1
##-   if (lIseq) {
##-     ## is the seqence represented by any other variable?
##-     lseqvar <- if (length(plargs$xvar)>0)
##-       sapply(plargs$pldata[,plargs$xvar,drop=FALSE],function(x) {
##-            if (is.factor(x)||is.character(x)) FALSE else {
##-              ld <- diff(x)
##-              sum(ld==0)<0.1*length(x) && (all(ld<=0) | all(ld>=0)) }
##-            } ) else FALSE
##-     if (any(lseqvar)) {
##-       warning(paste(":plot.regr / lpresx: sequence represented by",
##-                     paste(plargs$xvar[lseqvar],collapse=", ")))
##-       lIseq <- FALSE
##-     ## otherwise, plresx_ will plot against
##-     }
##-   }
  ## -------------------------------------------
  lxvar <- attr(lpldata, "xvar")
  if (!sequence) lxvar <- setdiff(lxvar, "(sequence)")
  plargs$mf <- FALSE ## avoid a new page
  ## plargs$ylim <- lylim  ## no need to calculate again
  if (length(i.def(lxvar, NULL, valuefalse=NULL)))
    plresx(x, data=data, resid=lres, xvar=lxvar, 
           transformed = transformed, sequence=sequence,
           weights= if ("weights"%in%names(lplsel)) FALSE else NULL,
           addcomp = addcomp, smooth.legend=lsmlegend,
           plargs = plargs)
  ## --- end
  invisible(plargs)
}
## ==========================================================================
i.plotselect <-
  function(plotselect, smooth=2, Iwgt = FALSE, mult = FALSE, 
           famgauss = TRUE, famglm = FALSE, famcount = FALSE)
{
  ## plot selection
  lsmdef <- 1+smooth-famglm
  lplsel <- c( yfit=0, resfit=lsmdef, absresfit = NA,
              absresweights = NA, qq = NA,
              leverage = 1, resmatrix = 1, qqmult = 1)
  if (length(plotselect)>0) {
    lplotsel <- unlist(plotselect)
    if (length(lplotsel)!=length(plotselect)) {
      warning(":i.plotselect: unsuitable argument.",
              "'plotselect' should be a vector")
      return(lplsel)
    }
    lpls <- TRUE
    lplnm <- names(plotselect)
    if (length(lplnm)==0) {
      if (length(plotselect)==length(lplsel))
        lplnm <- names(lplsel)
      else {
        warning(":plot.regr: Inadequate argument plotselect")
        lpls <- FALSE}
    }
    if (lpls) {
      if ("default"%in%lplnm) {
        lplnm <- setdiff(lplnm, "default")
        lplsel[] <- if (plotselect["default"]==0) 0 else
        pmin(lplsel,plotselect["default"])
      }
      lina <- is.na(match(lplnm,c(names(lplsel),"default")))
      if (any(lina)) {
        warning(":plot.regr: Inadequate elements in plotselect: ",
              paste(names(plotselect)[lina], collapse=", "))
        lplnm <- lplnm[!lina] }
      lplsel[lplnm] <- plotselect[lplnm]
    }
  }
  if (!mult) lplsel[c("resmatrix","qqmult")] <- 0
  if (is.na(lplsel["yfit"])) lplsel["yfit"] <- 0
  if (is.na(lplsel["resfit"]))
    lplsel["resfit"] <- lsmdef * (lplsel["yfit"]==0)
  if (is.na(lplsel["absresfit"]))
    lplsel["absresfit"] <- 2 * !famcount
  if (is.na(lplsel["absresweights"]))
    lplsel["absresweights"] <- 2 * (famgauss&Iwgt)
  if (is.na(lplsel["qq"])) lplsel["qq"] <- famgauss # how about gamma? !!!
  lplsel
}
## ==========================================================================
plresx <-
  function(x, data = NULL, xvar = TRUE, transformed = FALSE,
           sequence = FALSE, weights = NULL, 
           addcomp = FALSE, smooth.legend = FALSE, markextremes = NA,
           plargs=NULL, ...)
## ------------------------------------------------------------
{ ## plresx
  lcall <- match.call()
  ## lIplargs <- is.null(lcall[["plargs"]])
  if (length(plargs)==0) {
  ## plargs <-
  ##  if (lIplargs) {
    lac <- as.list(lcall)[-1]
    ladrop <- c("sequence", "weights", "addcomp", "smooth.legend")
    lcall <-
      c(list(quote(plotregr.control)), lac[setdiff(names(lac), ladrop)])
    mode(lcall) <- "call"
    plargs <-eval(lcall, parent.frame())
  } ## else  eval(lcall[["plargs"]]) ## , sys.frame(1))
  ploptions <- plargs$ploptions
  ## -------------------------------------------------------------------
  if (inherits(x,"mulltinom"))
    stop("!plresx! I do not know how to plot residuals of a mulitnomial regression")
  lres <- plargs$residuals
  pldata <- plargs$pldata
##  lIcq <- inherits(lres, "condquant")
  lform <- plargs$formula
  lvars <- attr(pldata, "xvar")
  lrawv <- u.allvars(lvars)
  lsimres <- plargs$simres
  lnsims <- if (is.null(lsimres)) 0 else ncol(lsimres)
  lInnls <- !inherits(x, "nls")
  lnr <- nrow(pldata)
  lIwgt <- length(plargs$weights)>0 &
    !inherits(x, "glm") ## do not plot against weights for binom,...
  lnaaction <- plargs$na.action
  lweights <- naresid(lnaaction, plargs$weights)
  ##  --- sequence
  lIseq <- i.def(sequence, FALSE, TRUE, FALSE)
  if (lIseq) {
    if (length(lvars)) {
      ## is the seqence represented by any other variable?
      lseqvar <-
        if (length(lvars)>0)
          sapply(pldata[,lvars,drop=FALSE],function(x) {
            if (is.factor(x)||is.character(x)) FALSE else {
             ld <- diff(x)
             sum(ld==0)<0.1*length(x) && (all(ld<=0) | all(ld>=0)) }
           } ) else FALSE
      lIseq <- !any(lseqvar)
      if (!lIseq) warning(paste(":plresx: sequence represented by",
                                paste(lvars[lseqvar],collapse=", ")))
    }
    pldata$"(sequence)" <- structure(1:lnr, varlabel="sequence")
    lvars <- c(lvars,"(sequence)")
  }
  ##  --- weights  as x variable
  lIweights <- i.def(weights, lIwgt, TRUE, FALSE)
  if (lIweights)
    if (!lIwgt) 
      warning(":plresx; No weights found.",
              " Cannot plot residuals against weights")   else {
    pldata[,"(weights)"] <- naresid(lnaaction, plargs$weights)
    lvars <- c(lvars, "(weights)")
  }
  ## ------------------
  lnvars <- length(lvars)
  if (lnvars==0) {
    warning(":plresx: I did not find any x variables")
    return() }
  ## terminmodel
  lvmod <- all.vars(formula(x)[-2])
  if (transformed) {
    ltrms <- rownames(attr(terms(x),"factors"))
    lvrs <- u.allvars(ltrms) ## which terms contain which vars?
    litrms <- sapply(lvrs, function(x) any(x%in%c(lvars,ltrms)) )
    lvmod <- ltrms[litrms | ltrms%in%lvars]
  } else {
    lvi <- pmatch("(", lvars, nomatch=0)
    lvars <-
      if (any(lvi>0)) c(unlist(u.allvars(lvars[-lvi])),lvars[lvi])
      else  unlist(u.allvars(lvars))
  }
  terminmodel <- lvars%in%lvmod
  lvcomp <- intersect(lvars, lvmod)
  ##  reference lines
  refline <- i.def(plargs$refline, !inherits(x,"coxph"))
  ## type
  addcomp <- as.logical(i.def(addcomp, FALSE, TRUE, FALSE))
  ## lpa <- plargs
  lInnls <- !inherits(x, "nlls")
## -----------------------------------
  ## fit components for refline
  lIcomp <- addcomp|refline
  if (length(lvcomp)) {
    ##-     warning(":plresx: variables or terms not found -> no reference lines")
    ##-     refline <- FALSE
    if (any(terminmodel) && lInnls) {
      lcmp <- fitcomp(x, vars=lvcomp, transformed=transformed,
                      xfromdata=FALSE, se=i.def(plargs$reflineband, TRUE))
      lcompx <- lcmp$x
      lcompy <- if (addcomp) lcmp$comp else -lcmp$comp
      lcompse <- lcmp$se
      if (addcomp) {
        lcompdt <-
          fitcomp(x, pldata, vars=lvcomp, transformed=transformed,
                  xfromdata=TRUE)$comp
        ## !!! add to lres, careful for condq
      }
    } else refline <- FALSE
  }
  ## quantile 
  lqnt <-
    if (length(plargs$dfres)>0) {
      qt(1-plargs$testlevel/2, plargs$dfres)
    } else  qnorm(1-plargs$testlevel/2)
  ## --- smooth
  lIsmooth <- i.def(plargs$smooth, i.getploption("smooth"),TRUE)
  if (length(names(smooth.legend))==0) {
    lsmlegend <- i.def(smooth.legend, NULL, TRUE, NULL)
    if (length(lsmlegend)==1)
      lsmlegend <- setNames(lsmlegend, lvars[1])
  } else lsmlegend <-
           if("(xvar)"%in%names(smooth.legend))
             setNames(smooth.legend, lvars[1]) else smooth.legend
  ## --- multivariate ## !!! factors!
  if (inherits(x,"mlm")) {
    lpanel <- function(xx, yy, indx, indy, pch, col, plargs=NULL, ...) {
      lcmpx <- lcmpy <- NULL
      ltin <- terminmodel[indx]
      lvx <- lvars[indx]
      lcnt <- !is.factor(lvx)
      if (ltin&refline) {
        lcmpy <- lcompy[,lvx,indy]
        if (lcnt) lcmpx <- lcompx[,lvx]
      }
      lsm <- if(is.factor(xx)) FALSE else lIsmooth
      if (lsm) plsmooth(xx, yy, plargs=plargs)
      plpoints(xx, yy, plargs=plargs)
    }
    plmatrix(pldata[,lvars,drop=FALSE],lres, panel=lpanel,
             ##pch=plargs$plab, plcol=plargs$pldata$plcol,
             nrow = plargs$multnrow, ncol = plargs$multncol, mar = plargs$multmar,
             plargs=plargs) 
    return()
  }
  ## ------------------------------------------------------------------
  lmf <- plargs$mf
  if (length(lmf)) {
    if (is.logical(lmf)&&lmf)
      lmf <- if (lnvars<=6) lnvars else
        min(lnvars,ceiling(lnvars/((lnvars-1)%/%6+1)))
    }
  loma <- i.def(plargs$oma, c(2,1)*(length(lmf)>0), valuefalse=NULL)
  if (length(loma)<4) loma <- c(0,0,loma,0)[1:4]
  loldpar <- 
    if (length(lmf)&(!is.logical(lmf))) 
      attr(if (length(lmf)==1) plmfg(mft=lmf, oma=loma) else
          plmfg(lmf[1], lmf[2], oma=loma), "oldpar")
    ## else par(oma=loma)  ## , ask=plargs$ask
  if (length(loldpar)) on.exit(par(loldpar), add=TRUE)
  ##
  lmbox <- ploptions$factor.show=="mbox"
  lIjitter <- !lmbox ## ploptions$factor.show=="jitter"
  ljitfac <- i.getploption("jitter.factor")
  lrpl <- plargs$resplab
  if (lIrpl <- length(lrpl)>0) lrpl <- lrpl[,1]
  lpla <- plargs
  ## --- loop --- plresx
  for (lj in 1:lnvars) {
    lv <- unname(lvars[lj])
    lvr <- lvars[lj] ## if (transformed) lvars[lj] else lrawv[lj]
    lcmpj <- terminmodel[lj] && refline && lInnls
    lci <- if (lcmpj&refline) lcompy[, lvr] else 0 ## !!!
    rr <- lres
    if (inherits(rr, "condquant")) rr <- rr[,1,drop=FALSE]
    if (plargs$partial.resid) 
      if (refline && addcomp && lcmpj) {
        rr <- rr+lcompdt[, lvr]
      }
    lvv <- pldata[, lv]
    ##    lpa$pldata <- cbind(rr, lvv, pldata)
    mar <- i.def(plargs[["mar"]], c(NA, par("mar")[-1]))
    ## ---
    if (is.null(attr(lvv, "varlabel")))  attr(lvv, "varlabel") <- lv
    if (is.factor(lvv)) {
      ## factors
      ## lmar <- par("mar") ## might be adapted to nchar(levels)
      lvv <- structure(lvv, varlabel=attr(lvv, "varlabel"))
      ll <- levels(lvv)
      lnl <- length(ll)
      lrs <- rr[[1]]
      if (lmbox) {
        ##       if (lIcq) rr <- lres[,"random"]
        attr(lpla$pldata,"xvar") <- lv
        plmboxes(lvv, lrs, data=pldata, mar=mar, plargs=lpla)
      } else {
##        attr(lvv, "plcoord") <- jitter(as.numeric(lvv), factor=ljitfac)
        plframe(lvv,rr, ploptions=ploptions)
      }
      ## reference values
      if (lcmpj) {
        lx <- seq_along(ll)
        lcil <- lci[1:lnl]
        lrsrg <- attr(lrs,"innerrange")
        if ((!is.null(lrsrg))&&diff(lrsrg)) {
          lcilp <- plcoord(lcil, lrsrg, innerrange.ext=plargs$innerrange.ext)
          if (any(attr(lcilp,"nmod"))) {
            liout <- lcilp!=lcil
            lrsout <- pmax(pmin(lcilp[liout], lrsrg[2]), lrsrg[1])
            segments(lx[liout]-0.4, lrsout, lx[liout]+0.4, lrsout,
                     lty=ploptions$refline.lty, lwd = ploptions$refline.lwd,
                     col=ploptions$refline.col[1])
            lcil[liout] <- NA
          }
        }
        segments(lx-0.4, lcil, lx+0.4, lcil,
                 lty=ploptions$refline.lty[1],
                 lwd = ploptions$refline.lwd[1],
                 col=ploptions$refline.col[1])
        if (refline && plargs$reflineband) {
          wid <- lqnt * lcompse[1:lnl, lv]
          lcill <- pmax(lcil-wid,lrsrg[1])
          lcilu <- pmin(lcil+wid,lrsrg[2])
          lines(c(rbind(lx-0.1,lx+0.1,lx+0.1,lx-0.1,lx-0.1,NA)),
                c(rbind(lcill,lcill,lcilu,lcilu,lcill,NA)),
                lty=ploptions$refline.lty[2],
                lwd = ploptions$refline.lwd[2],
                col=last(ploptions$refline.col))
        }
      }
      if (!lmbox) plpoints(lvv,rr, plargs=plargs)
    } else { # ---
## --- continuous explanatory variable
      plframe(lvv,rr, ploptions=ploptions)
      if (lIsmooth) {
        if (lnsims>0)
          rr <- cbind(rr, if (refline && addcomp && lcmpj)
                            lsimres+lcompdt[, lvr] else lsimres )
##-         if (inherits(rr, "condquant"))
##-           rr <- rr[,1]
        plsmooth(lvv, rr, plargs=plargs, band=plargs$smooth>=2) ## was lpa
      }
      ## refline
      if (refline && lcmpj) {
        lrefx <- lcompx[,lvr]
        lrefyb <-
          if (plargs$reflineband)
            outer(lqnt*lcompse[,lvr], c(-1,1))  else  NULL
        plrefline(list(x=lrefx, y=lci, band=lrefyb), x=lvv, y=rr,
                   plargs=plargs)
      }
      ## points
      plpoints(lvv,rr, plargs=plargs, plab=if (lIrpl) lrpl)
    } ## ends  if factor else 
  }
  invisible(plargs)
}
## ==========================================================================
smoothRegr <-
  function(x, y, weights=NULL, par=5/length(x)^0.3, iterations=50,
           minobs=NULL, ...)
{
  minobs <- i.def(minobs, i.getploption("smooth.minobs"), valuefalse=NULL)
  if (length(x)-sumna(x)<minobs ) return(NULL)
  iterations <- max(iterations, 1)
  lform <- if (is.formula(x)) x else y~x
  ## ----------------------------------------------------------------
  lcall <- call("loess", formula=lform, data=data.frame(x=x, y=y),
              weights=weights, span=par, iterations=iterations,
              family=if (iterations>0) "symmetric" else "gaussian",
              na.action=na.exclude)
  if (is.null(weights)) lcall$weights <- NULL
  lsm <- if (u.debug()) eval(lcall, parent.frame())
         else try(eval(lcall), parent.frame(), silent=TRUE)
  if (class(lsm)=="try-error") {
    warning(":smoothRegr: span was too small. Using 0.99")
    lcall$span <- 0.99
    lsm <- eval(lcall, parent.frame())
  }
  fitted(lsm)
}
## ========================================================================
gensmooth <-
  function(x, y, band=FALSE, power=1, resid="difference",
           plargs=NULL, ploptions=plargs$ploptions, ...)
{
  ## Purpose:   smooth for multiple y : one column from data, the other sim
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  9 Feb 2016, 14:57
  ploptions <- plargs$ploptions
  lsmfunc <- i.getploption("smooth.function")
  if (is.character(lsmfunc)) lsmfunc <- get(lsmfunc)
  if (is.null(lsmfunc)) lsmfunc <- smoothRegr
  power <- i.def(power, 1,1,1)
  ## ---
  lnx <- NROW(x)
  ly <- as.matrix(y)
  if (nrow(ly)!=lnx) stop("!gensmooth! Incompatible dimensions of 'x' and 'y'")
  ## if (length(weights)<=1) weights <- rep(1, lnx)
  lweights <- plargs$pldata$"(smoothWeights)"
  lIwgt <- length(lweights)>0
  if (lIwgt&&length(lweights)!=lnx)
    stop("!gensmooth! Incompatible dimensions of 'x' and 'weights'")
  lgroup <- plargs$pldata$"(smooth.group)"
  if (is.null(lgroup)) lgroup <- plargs$pldata$"(group)"
  if (lInogrp <- length(lgroup)<=1) lgroup <- rep(1, lnx)
  if (length(lgroup)!=lnx)
    stop("!gensmooth! Incompatible dimensions of 'x' and 'group'")
  lgrp <- factor(lgroup)
  if (is.character(resid))
    resid <- pmatch(resid, c("difference","ratio"))
  if (is.na(resid)) {
    warning(":gensmooth: argument 'resid' not suitable.",
            " Difference residuals calculated")
    resid <- 1
  }
  lnobs <- median(table(lgrp))
  band <- i.def(band, FALSE, TRUE)
  lpar <- i.getploption("smooth.par")
  if (is.function(lpar)) lpar <- smoothpar(lnobs)
  lparband <- lpar[1]* i.def(lpar[2], 1.5, 1.5, 1)
  liter <- i.getploption("smooth.iter")
  ## data: look for numvalues
  lx <- i.def(attr(x,"numvalues"),x)
##  if (NCOL(ly)==1) ly <- as.matrix(i.def(attr(y,"numvalues"),ly))
##  if (any(apply(ly,2, function(y) length(attr(y, "numvalues"))>0)))
  ly <- apply(ly,2, function(y) i.def(attr(y, "numvalues"), y))
  lnna <- apply(cbind(lx,ly), 1, sumna)
  lx[lnna>0] <- NA
  lio <- order(as.numeric(lgrp), lx) ## order by group
  lio <- lio[!is.na(x[lio])]
  lxo <- lx[lio] # sorted without NA
  lyo <- ly[lio,,drop=F]
  lgrpo <- lgrp[lio]
  lgrpn <- as.numeric(lgrpo)
  lwgto <- if(lIwgt) lweights[lio] else NULL
  ## production
  oldopt <- options(warn=-1)
  on.exit(options(oldopt))
  lysm <- array(NA, dim=dim(lyo), dimnames=dimnames(lyo))
  ## presently only for matrices
  if (band) lysmband <- lsmrpos <- lysm[,1]
  for (lgr in seq_along(levels(lgrpo))) {  ## smooth within groups (if >1)
    lig <- which(lgrpn==lgr)
    for (j in ncol(lyo):1) {
      lsm <- lsmfunc(lxo[lig], lyo[lig,j]^power,
                     weights=if(lIwgt) lwgto[lig] else NULL,
                     par=lpar[1], iterations=ploptions$smooth.iter, ...)
      if (length(lsm)) lysm[lig,j] <- lsm^(1/power)
    }
    if (band & length(lsm)) {
      lysmb <- rep(0, length(lsm))
      lsmr <- lyo[lig,1]-lsm^(1/power) ## residual
      ## high end
      lii <- lsmr>=0
      lsmrh <- lsmr[lii]
      ligi <- lig[lii]
      lsmh <- lsmfunc(lxo[ligi], sqrt(lsmrh),
                      weights=if (lIwgt) lwgto[ligi] else NULL,
                      par=lparband, iterations=liter)
      if (length(lsmh)) lysmb[lii] <- lsmh^2
      ## low end
      lii <- lsmr<=0
      lsmrl <- - lsmr[lii]
      ligi <- lig[lii]
      lsml <- lsmfunc(lxo[ligi], sqrt(lsmrl),
                      weights=if (lIwgt) lwgto[ligi] else NULL,
                      par=lparband, iterations=liter)
      if (length(lsml)) lysmb[lii] <- - lsml^2
      ## resulting band
      lysmband[lig] <- lysmb + lsm
      lsmrpos[lig] <- !lii
    }
  }
  lysmin <- matrix(NA, lnx, ncol(lyo), dimnames=list(names(x),colnames(lyo)))
  lysmin[lio,] <- lysm
  lres <- if (resid==2) ly/lysmin else ly-lysmin
  rr <- list(x = lxo, y = lysm, group = if(!lInogrp) factor(lgrpo),
             index = lio, xorig = x, ysmorig = lysmin, residuals = lres)
  if (band) rr <- c(rr, yband = list(lysmband), ybandindex = list(lsmrpos) )
  rr
}
## ==========================================================================
smoothLm <- function(x, y, weights = NULL, ...) {
  if (is.null(weights)) lm.fit(cbind(1,x), y, ...)$fitted
  else  lm.wfit(cbind(1,x), y, weights, ...)$fitted
}
## smoothRegrrob <- function(x,y,weights,par=3*length(x)^log10(1/2),iter=50)
## =======================================================================
plmatrix <-
function(x, y=NULL, data=NULL, panel=plpanel, 
         nrow=0, ncol=nrow, reduce=TRUE, keeppar=FALSE,
         xaxmar=NULL, yaxmar=NULL, xlabmar=NULL, ylabmar=NULL,
         xlab=NULL, ylab=NULL, ## partial match!?!
         oma=NULL, mar=NULL, cex.diaglabels = 1.5, ...) 
{
## Purpose:    pairs  with different plotting characters, marks and/or colors
##             showing submatrices of the full scatterplot matrix
##             possibly on seeral pages
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: 23 Jul 93; minor bug-fix+comments:
  lf.axis <- function(k, x, axm, labm, txt, ...) {
    if (k %in% axm) plaxis(k, x, varlabel="", ploptions=ploptions)
    if (k %in% labm)
      mtext(txt, side=k, line=(1.5+(k %in% axm)), ...)
  }
  lf.eq <- function(v1,v2) {
    if (is.factor(v1)) is.factor(v2)&& all(notna(as.numeric(v1)==as.numeric(v2)))
    else all(notna(v1==v2))
  }
  ## ---
  oldpar <- par(c("mfrow","mar","cex","mgp")) ##, "ask"
##  lmfg <- par("mfg")
  if (!keeppar) on.exit(par(oldpar))
##---------------------- preparations --------------------------
##-   lcl <- match.call()
##-   lcall <- sys.call() ## match.call modifies argument names
##-   lcnm <- names(lcall)
##-   if (length(lcall)!=length(lcl)) stop("bug")
  ##-   names(lcall) <- ifelse(lcnm=="", names(lcl), lcnm)
  lcall <- match.call()
  lIplargs <- is.null(lcall[["plargs"]])
  if (lIplargs) {
    lcall[[1]] <- quote(pl.control)
    lcall$x <- x  ## needs evaluation
    lcall$.subdefault <- as.character(substitute(x))
    lcall$y <- y
    lcall$data <-
      if(length(data)) {
        if (is.name(substitute(data)))
          lcall$.subdefault <- as.character(substitute(data))
        data
      }
##    lcall$mar <- rep(i.def(mar, i.getploption("panelsep")), length=4 )
    lcall$gensequence <- FALSE
    lcall$.environment. <- parent.frame()
    ##
    plargs <- eval(lcall, envir=parent.frame())
##    plargs$main <- lmain
  } else    plargs <- eval(lcall[["plargs"]], sys.frame(1))
  ploptions <- plargs$ploptions
  ## margins, will be used by  plframe
  ploptions$mar <- rep(i.def(mar, i.getploption("panelsep")), length=4 )
  if (length(xlab)|length(ylab))
    warning(":plmatrix: Arguments 'xlab' and 'ylab' not used. ",
            "Set 'varlabels' instead!")
  ## -----------------------------------------------
  ## data and bookkeeping
  pldata <- plargs$pldata
  if (is.null(plargs$main)) plargs$main <- plargs$dataLabel
  xvar <- i.def(attr(pldata,"xvar"), names(pldata))
  nv1 <- length(xvar)
  lv1 <- lv2 <- 0
  if (is.null(y)) {
    xvar <- c(xvar, attr(pldata,"yvar"))
    nv1 <- length(xvar)
    if (reduce) { nv1 <- nv1-1; lv2 <- 1 }
    nv2 <- nv1
    ldata <- pldata[,xvar]
  } else { # cbind y to data for easier preparations
    reduce <- FALSE
    if (!is.null(dim(y))) {
      yvar <- colnames(y)
      ldata <- cbind(pldata[,xvar], as.data.frame(y))
    } else {
      yvar <- attr(pldata,"yvar")
      ldata <- pldata[,c(xvar,yvar)]
    }
    lformy <- as.formula(paste("~",paste(yvar, collapse="+")))
    nv2 <- length(yvar)
    lv2 <- length(xvar)
  }
  nvv <- ncol(ldata)
  lnr <- nrow(ldata)
  lnobs <- lnr-mean(sumna(ldata))
  lvsurv <- sapply(ldata, function(x) inherits(x, "Surv") )
  if (any(lvsurv)) { ## survival vars
    lf.surv <- function(dt) structure(dt[,1], pch=dt[,2]+1)
    ldata[lvsurv] <- lapply(ldata[lvsurv], lf.surv)
  }
  ## title !!!
  lsub <- plargs$sub
  lmain <- plargs$main
  lImain <- !( (length(lmain)==0||lmain=="") & is.logical(lsub)&&lsub )
  ## --- position of tick marks and axis labels, oma
  xaxmar <- i.def(xaxmar, 1+(nv1*nv2>1))
  xaxmar <- ifelse(xaxmar>1,3,1)
  yaxmar <- i.def(yaxmar, 2+(nv1*nv2>1))
  yaxmar <- ifelse(yaxmar>2,4,2)
  xlabmar <- i.def(xlabmar, if (nv1*nv2==1) xaxmar else 4-xaxmar )
  ylabmar <- i.def(ylabmar, if (nv1*nv2==1) yaxmar else 6-yaxmar )
  lcexmain <- i.getploption("title.cex")
  if (length(oma)!=4)
    oma <- c(2+(xaxmar==1)+(xlabmar==1), 2+(yaxmar==2)+(ylabmar==2),
             1.5+(xaxmar==3)+(xlabmar==3)+lcexmain[1]*lImain,
             2+(yaxmar==4)+(ylabmar==4))
  ## set par
  ## if (!keeppar)
  plmfg(nrow, ncol, nrow=nv2, ncol=nv1, oma=oma, mar=mar) # , mgp=c(1,0.5,0)
  ## else par(mfrow=par("mfrow")) ## new page! if mfcol was set, it will not work
  lmfig <- par("mfg")
  lnr <- lmfig[3]
  lnc <- lmfig[4]
  lnpgr <- ceiling(nv2/lnr)
  lnpgc <- ceiling(nv1/lnc)
  ## cex
  lcex <- ploptions$cex
  if (is.function(lcex)) lcex <- lcex(lnobs)
  plargs$ploptions$cex <- lcex
  plargs$ploptions$axes <- FALSE
##
##-   ## log
##-   if (length(grep("x",log))>0) ldata[ldata[,1:nv1]<=0,1:nv1] <- NA
##-   if (length(grep("y",log))>0) ldata[ldata[,lv2+1:nv2]<=0,lv2+1:nv2] <- NA
  ##----------------- plots ----------------------------
  for (ipgr in 1:lnpgr) {
    lr <- (ipgr-1)*lnr ## start at row index  lr
  for (ipgc in 1:lnpgc) {
    lc <- (ipgc-1)*lnc
    if (reduce&&((lr+lnr)<=lc)) break
  for (jr in 1:lnr) { #-- plot row [j]
    jd2 <- lr+jr  ##  index for  y  axis
    j2 <- lv2 + jd2
    if (jd2<=nv2)  v2 <- ldata[,j2]
    lylab <- i.def(attr(v2,"varlabel"), paste("V",j2,sep="."), valuefalse="")
    for (jc in 1:lnc) { #-- plot column  [j2-lv2] = 1:nv2
      jd1 <- lc+jc
      j1 <- lv1 + jd1
      if (jd2<=nv2 & jd1<=nv1) {
        v1 <- ldata[,j1]
        lxlab <- i.def(attr(v1,"varlabel"), paste("V",j1,sep="."), valuefalse="")
        if (!lf.eq(v1,v2)) { # not diagonal
          plframe(v1, v2, xlab="", ylab="", ticklabels = FALSE,
                  ploptions=ploptions) # plargs=plargs
          panel(v1,v2, indx=jd1, indy=jd2, plargs=plargs)
        }
        else {
          v1 <- as.numeric(v1)
          plot(v1,v1, type="n", axes=FALSE, xlab="",ylab="") ## frame()
          uu <- par("usr") # diagonal: print variable name
          text(mean(uu[1:2]),mean(uu[3:4]), lylab, cex=cex.diaglabels)
        }
      ## if (axes) {
        usr <- par("usr")
        lat=c(mean(usr[1:2]),mean(usr[3:4]))
        if (jr==lnr||jd2==nv2) lf.axis(1, v1, xaxmar, xlabmar, lxlab)
        if (jc==1) lf.axis(2, v2, yaxmar, ylabmar, lylab)
        if (jr==1) lf.axis(3, v1, xaxmar, xlabmar, lxlab)
        if (jc==lnc||jd1==nv1) lf.axis(4, v2, yaxmar, ylabmar, lylab)
        ## }
        ## box()
        if (lImain)   pltitle(plargs=plargs, show=NA)
      } else frame()
    }
  }}
    ## stamp(sure=FALSE, outer.margin=TRUE) 
  }
}

## ====================================================================
plpanel <-
  function(x = NULL, y = NULL, indx=NULL, indy=NULL, type="p", frame = FALSE,
           title = FALSE, plargs = NULL, ploptions = plargs$ploptions,
           ...) ## data=plargs$pldata
{
  if (is.null(plargs)) 
    plargs <- get(".plargs", ".GlobalEnv") ## list in calling fn
  pldata <- plargs$pldata
  plargs$ploptions$stamp <- FALSE
  mbox <- i.getploption("factor.show")=="mbox"
  lIsm <- i.getploption("smooth")
  ## intro, needed if formulas are used or data is given or ...
  lcl <- match.call()
  lcall <- sys.call() ## match.call modifies argument names
  lcnm <- i.def(names(lcall), names(lcl))
  names(lcall) <- ifelse(lcnm=="", names(lcl), lcnm)
  if (is.formula(x)|is.formula(y)|any(c("data","pcol","psize")%in%lcnm)) {
    lcall$assign <- FALSE
    lcall$ploptions <- ploptions
    if (is.null(lcall$data)) lcall$data <- pldata
    lplargs <- do.call(pl.control, as.list(lcall[-1]), envir=parent.frame())
    ploptions <- lplargs$ploptions
    plargs$pldata <- pldata <- lplargs$pldata
    x <- pldata[,2]
    y <- pldata[,1]
  } else {
  ## ---
    if (length(x)==0) x <- pldata[,2]
    if (length(y)==0) y <- pldata[,1]
  }
  if (is.character(x)) x <- factor(x) ## !!! attributes!
  if (is.character(y)) y <- factor(y)
  if (is.factor(x)) { lIsm <- FALSE
    if (!is.factor(y) & mbox) {
      plargs$pldata <- data.frame(x=x,y=y)
      plmboxes.default(x, y, plargs=plargs, add=TRUE, ...)
      return()
    }
  } else { 
    if (is.factor(y) & mbox) {
      plargs$pldata <- data.frame(x=y, y=x)
      plmboxes.default(y, x, plargs=plargs, add=TRUE, horizontal=TRUE, ...)
      return()
    }
  }
  ## ---
  if (frame) plframe(x,y, ploptions=ploptions)
  if (lIsm)  plsmooth(x,y, plargs=plargs)
  ## refline
  if (length(lrfl <- i.getploption("refline")))
    plrefline(lrfl, x=x, y=y, plargs=plargs)
  plpoints(x, y, type=type, ## plab = plab, pch = pch, col = col,
           ##lty = lty, lwd = lwd, psize = psize, condquant = condquant,
           plargs=plargs, ploptions=ploptions, ...)
  if (title) pltitle(plargs=plargs, show=title)
}
## ====================================================================
panelSmooth <-
  function(x, y, indx, indy, plargs=.plargs, ploptions=plargs$ploptions, ...)
    graphics::panel.smooth(x, y, pch=plargs$pch, col=plargs$pcol,
                           cex=i.getploption("cex.pch"), ...)
## ====================================================================
plmbox <-
  function(x, at=0, probs=NULL, outliers=TRUE, na.pos=NULL,
           horizontal = FALSE,
  width=1, wfac=NULL, minheight= NULL, adj=0.5, extquant=TRUE, 
  widthfac=c(max=2, med=1.3, medmin=0.3, outl=NA),
  colors=c(box="lightblue2",med="blue",na="gray90"), lwd=c(med=3, range=2),
  warn=options("warn") )
{
  ## Purpose:   multi-boxplot
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 14 Dec 2013, 10:09
  f.box <- function(wid, quant, col, lwmax) {
  ##  if (is.na(col)||col==0) col="white"
    if (wid>lwmax) {
      if (horizontal) {
        polygon(quant, at+lpos*lwmax, col="black")
        polygon(quant, at+lpos*lwmax^2/wid, col=col)
      } else {
        polygon(at+lpos*lwmax, quant, col="black")
        polygon(at+lpos*lwmax^2/wid, quant, col=col)
      }
    } else {
      if(wid>0)
      {
        if (horizontal)
          polygon(quant, at+lpos*wid, col=col) else
          polygon(at+lpos*wid, quant, col=col)
      }
    }
  }
  ## -----------------------------------
  lq <- lwid <- NULL
  lfac <- 0
  lx <- x[!is.na(x)]
  if (length(lx)==0) {
    if (warn>=0) warning(":plmbox: no non-missing data")
    return(structure(numeric(0), q=lq, width=lwid))
  } else { ## ----
  stopifnot(length(width)==1,length(wfac)<=1)
  if (is.null(probs))
      probs <- if (sum(!is.na(x))<20) c(0.1,0.5,1)/2 else
               c(0.05,0.1,0.25,0.50,0.75,1)/2
  lprobs <- if (all(probs<=0.5))  c(probs,1-probs)  else c(probs)
  lprobs <- sort(unique(lprobs))
  colors <- as.list(colors)
  box.col <- colors[["box"]]
  if (length(box.col)==1)
    box.col <- ifelse(0.25<=last(lprobs,-1) & lprobs[-1]<=0.75, box.col, NA)
  ## values for degenerate case
  lxsd <- IQR(lx)
  lfac <- if (is.null(wfac)) width*2*lxsd else wfac*length(lx)
                                        # was mad/dnorm(0)
  lmed <- median(lx)
  lwmed <- width
  lrg <- range(lx)
  lirg <- i.def(attr(x, "innerrange"), lrg)
  ljrg <- any(lirg!=lrg)
  lirgd <- diff(lirg)
  loutl <- lx
  lwoutl <- widthfac["outl"]
  if (diff(lrg) > 0) { ## non-degenerate
    if (is.null(minheight))
      minheight <- if (lxsd==0) lirgd*0.02  else lxsd*0.01
    lqpl <- lq <- quinterpol(lx, probs=lprobs, extend=extquant)
    lirgext <- attr(x, "innerrange.ext")
    if (ljrg) {  ## transformed coord
      lx <- plcoordtrsf(lx, lirg, lirgext)
      lrg <- plcoordtrsf(lrg, lirg, lirgext)
      lqpl <- plcoordtrsf(lqpl, lirg, lirgext)
    }
    loutl <- lx[lx<min(lq)|lx>max(lq)]
  ## ---
    lwid <- lfac*diff(lprobs)/pmax(diff(lq), minheight)
##    lxsd <- IQR(x, na.rm=TRUE)
    lwmax <- widthfac["max"]*lfac*0.5/ifelse(lxsd>0, lxsd, 1)
    lwmed <- max(widthfac["med"]*min(lwmax,max(nainf.exclude(lwid))),
                 widthfac["medmin"],na.rm=TRUE)
    lpos <- c(-adj,-adj,1-adj,1-adj)
    if (is.na(lwoutl)) lwoutl <- 0.1*lwmax
  ## ---
    for (li in 1:(length(lprobs)-1)) 
      f.box(lwid[li], lqpl[li+c(0,1,1,0)], box.col[li], lwmax)
  }  ## 
  ## median
  if (horizontal) {
    lines(rep(lmed,2), at+lwmed*c(-adj,1-adj), col=colors[["med"]],
          lwd=lwd["med"])
    lines(lrg, c(at,at), # +linepos*0.01*diff(par("usr")[1:2])*(0.5-adj),
          lwd=lwd["range"])
  } else {
    lines(at+lwmed*c(-adj,1-adj), rep(lmed,2), col=colors[["med"]],
          lwd=lwd["med"])
    lines(c(at,at), # +linepos*0.01*diff(par("usr")[1:2])*(0.5-adj),
          lrg, lwd=lwd["range"])
  }
  if (outliers&&length(loutl)) {
    lat <- rep(at,length(loutl))
    if (horizontal)
      segments(loutl, lat-lwoutl*adj, loutl, lat+lwoutl*(1-adj))  else
      segments(lat-lwoutl*adj, loutl, lat+lwoutl*(1-adj), loutl)
  }
}
  if (!is.null(na.pos)) {
    lmna <- mean(is.na(x))
    if (lmna) {
      ldna <- diff(na.pos)
      if (length(ldna)==0 || is.na(ldna) || ldna==0)
        stop("!plmbox! argument 'na.pos' not suitable")
      lwidna <- lfac*lmna/abs(ldna)
      f.box(lwidna, na.pos[c(1,2,2,1)], colors[["na"]], lwmax) 
    }
  }
  invisible(structure(lfac/length(x), q=lq, width=lwid))
}
## ====================================================================
plmboxes <- function(x, ...)
  UseMethod("plmboxes")

plmboxes.formula <-
  function(x, data, ...)
{
  ldt <- genvarattributes(getvariables(x, data))
##-   l1asymbox <- length(x[[3]])>1 && as.character(x[[3]][[2]])=="1"
##-   if (l1asymbox)
##-     ldt <- data.frame(transferAttributes(ldt[,1]),0,
##-                       transferAttributes(ldt[,2]))
  ldtnm <- substitute(data)
  ldtnm <- if (is.name(ldtnm)) as.character(ldtnm) else format(ldtnm)
  plmboxes.default(ldt[,attr(ldt,"xvar"),drop=FALSE],
           ldt[,attr(ldt,"yvar"),drop=FALSE], ldt,
           .subdefault = if (length(ldtnm)<30) ldtnm else "", ...)
}
## ------------------------------------------------------------
plmboxes.default <-
  function(x, y, data, width=1, at=NULL, horizontal = FALSE,
    probs=NULL, outliers=TRUE, na=FALSE, asymbox=NULL,
    refline=NULL, add=FALSE, 
    xlim=NULL, ylim=NULL, axes=TRUE, xlab=NULL, ylab=NULL, 
    labelsperp=FALSE, xmar=NULL, 
    widthfac=NULL, minheight=NULL, colors=NULL, lwd=NULL,
    plargs=NULL, .subdefault=NULL,
    ...)
{
  ## Purpose:    multibox plot
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 14 Dec 2013, 23:38
  f.ylim <- function(ylm, ext)
    c((1+ext)*ylm[1]-ext*ylm[2], (1+ext)*ylm[2]-ext*ylm[1])
  lcall <- match.call()
  lIplargs <- is.null(lcall[["plargs"]])
  if (lIplargs) {
    lcall[[1]] <- quote(pl.control)
    lcall$x <- x  ## needs evaluation
    ## lcall$.subdefault <- .subdefault ## as.character(substitute(x))
    lcall$y <- y
    lcall$data <-
      if(length(data)) {
        if (is.null(.subdefault) & is.name(substitute(data)))
          lcall$.subdefault <- as.character(substitute(data))
        data
      }
    lcall$.environment. <- parent.frame()
    ##
    plargs <- eval(lcall, parent.frame())
  }
  ploptions <- plargs$ploptions
  title <- i.def(ploptions$title, TRUE)
  ## ----------------------------------------------------------------------
  pldata <- plargs$pldata
  x <- pldata[,i.def(attr(pldata,"xvar"),1), drop=FALSE] ## may have two columns
  y <- pldata[,i.def(attr(pldata,"yvar")[1],2)]
  lhoriz <- as.logical(i.def(horizontal, FALSE, valuetrue=TRUE))
  ## widths
  lwfac <- modarg(widthfac, c(max=2, med=1.3, medmin=0.3, outl=NA, sep=0.003))
  ## colors, line widths
  lcol <- modarg(colors,
                 c(box="lightblue",med="blue",na="gray90",refline="magenta") )
  llwd <- modarg(lwd, c(med=3, range=2))
  ## data
  ## preliminary 
  lx <- transferAttributes(factor(x[,1]),x[,1]) # unused levels are dropped
  llr <- ncol(x)>=2 ## asymmetrix mboxes required for binary factor
  if (llr && length(unique(x[,2]))!=2) {
    warning(":plmboxes: second x-variable must be binary. I ignore it.")
    x <- x[,1, drop=FALSE]
    llr <- FALSE
  }
  llist <- split(y,x)
  asymbox <-
    i.def(asymbox, length(unique(x[,1]))==2, valuetrue=TRUE) && ncol(x)==1
  if (asymbox) {
    llev <- ""
    llev2 <- c(levels(x[,1]),"","")[1:2]
    llr <- TRUE
  } else llev <- levels(lx)
  lng <- length(llev)
  lnn <- sapply(llist,length)
  lsd <- mean(sapply(llist,mad,na.rm=TRUE),na.rm=TRUE)
  width <- rep(width, length=lng)
  lfac <- width*lsd/(max(lnn)*(1+llr))
  if (is.null(minheight)) {
    lscales <- sapply(llist, IQR, na.rm=TRUE)
    minheight <- median(lscales)*0.02
  }
  ## labels
  xlab <- i.def(xlab, attr(x[,1],"varlabel"), valuefalse="")
  if (length(xlab)>1) xlab <- xlab[2]
##  if (xlab=="1") xlab <- ""
  ylab <- i.def(ylab, attr(y, "varlabel"), valuefalse="")
  ## position
  if (is.null(at)) at <- 1:lng else
    if (length(at)!=lng) {
      warning(":plmboxes: 'x' has wrong length")
      at <- at[1:lng] ## may produce NAs
  }
  ## probabilities
  if (is.null(probs))
      probs <- if (sum(!is.na(y))/(lng*(1+llr))<20) c(0.1,0.5,1)/2 else
               c(0.05,0.1,0.25,0.50,0.75,1)/2
  ## box for NA's?
  ##if (is.null(na)||is.na(na)||(is.logical(na)&&!na)) na.pos <- NULL else 
  ##    if (is.logical(na))
  na <- i.def(na, NA)
  ly <- i.def(attr(y, "plcoord"), y)
  na.pos <- i.def(na, c(min(ly, na.rm=TRUE)*(1-0.3)-0.3*max(ly, na.rm=TRUE)),
                  valuefalse=NULL)
  if (length(na.pos)==1)
    na.pos <- na.pos+ 0.03*diff(range(ly, na.rm=TRUE))*c(-1,1)
  lusr <- par("usr")
  ## plot range
  ##  lrg <- if (add) lusr[3:4] else attr(ly, "plrange")
  lirg <- attr(y, "innerrange")
  ljlim <- any(c(attr(y, "nmod"),0)>0)
  lyat <- i.def(attr(y,"ticksat"),
                pretty(ly, n=i.def(ploptions$tickintervals, 7)) )
  if (add) {
    xlim <- lusr[1:2]
    ylim <- lusr[3:4]
  } else {
    if(u.nuna(xlim)) xlim <- ## better: NAs -> default value
      range(at, na.rm=TRUE)+ max(width[c(1,length(width))])*c(-1,1)*0.5
    if (lIna <- !is.null(na.pos)) {
      lyat <- lyat[lyat>max(na.pos)]
      if (length(lyat)<3)
        lyat <- pretty(ly, n=i.def(ploptions$tickintervals, 7))
      attr(y, "plrange") <- range(c(attr(y, "plrange"), na.pos))
    }
    ylim <- i.def(ylim, attr(y, "plrange"))
    ## margins
    lmar <- i.getploption("mar")
    if (anyNA(lmar)) lmar <- ifelse(is.na(lmar), par("mar"), lmar)
    lxmardef <- lmar[1+lhoriz]
    ## the next statement defines the maximal label length as 10
    lmaxnchar <-
      ifelse(is.numeric(labelsperp), min(max(0,labelsperp),20), 10)
    lxmar <-
      i.def(xmar, ifelse(labelsperp,
                    2 + 0.6*min(max(nchar(llev)), lmaxnchar), lxmardef) )
    lmar[1+lhoriz] <- max(lxmar, lmar[1+lhoriz])
    plargs$ploptions$mar <- lmar
    plargs$ploptions$axes <- FALSE
    if (u.true(ploptions$grid))
      plargs$ploptions$grid <-
        if(lhoriz) list(TRUE, at) else list(at, TRUE)
    ## ---------------------------------
    if (lhoriz)
      plframe(y, xlim, plargs=plargs)
      else plframe(xlim, y, plargs=plargs)
    if (axes) {
      axis(1+lhoriz, at=at, labels=llev, las=2*as.logical(labelsperp))
      if(lIna && anyNA(y)) 
        mtext("NA", 2-lhoriz,1, at=mean(na.pos), las=2) 
      if (asymbox) {
        mtext(llev2[1], 1+lhoriz,1, at=0.75)
        mtext(llev2[2], 1+lhoriz,1, at=1.25)
      }
      mtext(xlab, side=1+lhoriz, line=ploptions$mar[1+lhoriz]-1)
      plaxis(2-lhoriz, x=y)
    }
  } # if (!add)
  ## ---
  if (!is.null(refline)) 
    if(lhoriz)
      abline(v=refline, col=ploptions$refline.col,
             lty=ploptions$refline.lty, lwd=ploptions$refline.lwd)
    else
      abline(h=refline, col=ploptions$refline.col,
             lty=ploptions$refline.lty, lwd=ploptions$refline.lwd)
  ## ---
  lusrd <- diff(par("usr")[1:2])
  lsep <- lwfac["sep"]*llr*lusrd
  lwoutl <- lwfac["outl"]
  if (is.na(lwoutl)) {
      lwoutl <- 0.05*lusrd
      lwfac["outl"] <- lwoutl/lng
  }
  if (llr) lwfac[c("medmin","outl")] <- lwfac[c("medmin","outl")] /2
  ## ------------
  lattr <- attributes(y)
  for (li in 1:lng) {
    if (is.na(at[li])) next
    lli <- llist[[li]]
    if (length(notna(lli) )) {
      attributes(lli) <- lattr
      plmbox(lli,at[li]+lsep, probs=probs, outliers=outliers,
             horizontal=lhoriz,
             wfac=lfac[li], adj=0.5*(1-llr), na.pos=na.pos, extquant=TRUE,
             widthfac=lwfac, colors=lcol, lwd=llwd, warn=-1)
      if (llr) { ## second half of asymmetrix  mbox
        if (length(llir <- llist[[li+lng]]))
          plmbox(llir,at[li]-lsep,probs=probs, outliers=outliers,
                 horizontal=lhoriz,
                 wfac=lfac[li], adj=1, na.pos=na.pos, extquant=TRUE,
                 widthfac=lwfac, colors=lcol, warn=-1)
      }
    }
  }
  pltitle(plargs=plargs, show=title)
  stamp(sure=FALSE)
  invisible(at)
}
## ---------------------------------------------------------------------------
plcoordtrsf <- function(x, innerrange, innerrange.ext)
{
  llx <- pmax(pmin(x,innerrange[2]),innerrange[1])
  lxd <- x-llx
  if (any(lxd!=0, na.rm=TRUE))
    x <- llx + lxd/(1+abs(lxd)/(diff(innerrange)*innerrange.ext))
  x
}

## ===========================================================================
plres2x <-
  function(formula=NULL, reg=NULL, data=NULL, restrict=NULL, size = 2,
  xlab = NULL, ylab= NULL, plextext = NULL, pale = 0.2, ...)
{
## Purpose:  plot residuals vs. two x`s
## Author:   ARu , Date:  11/Jun/91
## Aenderungen: MMae, 30/Jan/92, Dez.94
## --------------------------------------------------------------------------
## Arguments:
##   formula    z~x+y, where
##              x, y  coordinates of points given by two vector arguments.
##              z     gives orientation (by sign)
##                   and size (by absolute value) of symbol.
##   reg        regression results
##   data       data
##              you must specify either  reg  or  data
##   restrict   absolute value which truncates the size.
##              The corresponding symbols are marked by stars.
##   size       the symbols are scaled so that "size" is the size of
##              the largest symbol in cm.
##   main       main title, defaults to the formula
##   ...        additional arguments for the S-function `plot`
## the function currently only plots  z  for the first two terms of the
## right hand side of  formula
## --------------------------------------------------------------------------
  lcall <- match.call()
  lIplargs <- is.null(lcall[["plargs"]])
  plargs <-
    if (lIplargs) {
      lac <- as.list(lcall)[-1]
      lac$x <- lac$reg
      lac$xvar <- lac$formula
      ladrop <- c("formula", "restrict", "size")
      lcall <- c(list(quote(plotregr.control)),
                 lac[setdiff(names(lac), ladrop)])
      mode(lcall) <- "call"
      plargs <- eval(lcall, parent.frame())
    } else  eval(lcall[["plargs"]], sys.frame(1))
  ## ------------------------------------------------------------------
  ploptions <- plargs$ploptions
  lz <- plargs$residuals
  pldata <- plargs$pldata
  lxvar <- attr(pldata, "xvar")
  if (length(lz)==0) lz <- attr(pldata, "yvar")
  ##--- restrict z values: --- !!! use innerrange of z
  lrestr <-
    if (any(attr(lz, "nmod")>0)) max(abs(attr(lz, "innerrange"))) else NULL
  restrict <- i.def(i.def(restrict, lrestr), NULL, valuefalse=NULL)
  if(length(restrict)==0)   restr <- FALSE else {
    restr <- abs(lz) > restrict
    lz <- pmin( pmax( lz, -restrict), restrict) }
## size
  size <- i.def(size, 3)
  lratio <- size * par("cin")[1] / par("pin")[1]
  llwd <- i.getploption("lwd")
  lpcol <- pldata[["(pcol)"]]
  if (length(lpcol)) {
    lgrpcol <- i.getploption("group.col")
    if (is.factor(lpcol)) lpcol <- as.numeric(lpcol)
    if (is.numeric(lpcol)) lpcol <- rep(lgrpcol, length=max(lpcol))[lpcol]
  } else lpcol <- i.getploption("col")[1]
  pldata$"(pcol)" <- colorpale(lpcol, pale=pale)
  ## for ...
  lzj <- lz[,1]
  lzj <- lzj/max(abs(lzj), na.rm = TRUE)
  lpanel <-
    function(xx, yy, indx, indy, pch, col, plab, lwd, zz, ...) {
      lusr <- par("usr")
      lfx <- lratio * diff(lusr[1:2])
      lfy <- lratio * diff(lusr[3:4])
      lsxz <- c(lfx * abs(zz))
      lsyz <- c(lfy * zz)
      lx <- attr(xx, "plcoord")
      if (is.null(lx)) lx <- xx
      ly <- attr(yy, "plcoord")
      if (is.null(ly)) ly <- yy
      plpoints(xx, yy, plargs=plargs)
      segments(lx - lsxz, ly - lsyz,  lx + lsxz, ly + lsyz,
               lwd = llwd, col=col, xpd=TRUE)
      ##--- mark restricted observations: ---
      if(any(restr)) {
        points((xx - lsxz)[restr], (yy - lsyz)[restr], pch= 8, mkh = 1/40,
               xpd=TRUE)
        points((xx + lsxz)[restr], (yy + lsyz)[restr], pch= 8, mkh = 1/40,
               xpd=TRUE)
      }
    } 
  ljx <- lxvar[1]
  lx <- pldata[,ljx]
##-   llrg <- attr(lx, "plrange")
##-   attr(lx, "plrange") <- llrg + c(-1,1)*diff(llrg)* 2*lratio
  ljy <- lxvar[2]
  ly <- pldata[,ljy]
##-   llrg <- attr(ly, "plrange")
##-   attr(ly, "plrange") <- llrg + c(-1,1)*diff(llrg)* 2*lratio
  ##--
  ##---------------
  lplxx <- i.getploption("plextext")
  plframe(lx, ly, xlab=xlab[1], ylab=ylab[1], plextext=lplxx,
          ploptions=ploptions)
  ##--- draw symbols: ---
  lpanel(lx, ly, zz=lzj, lwd=llwd, col=lpcol)
  lmain <- plargs$main
  if (length(lmain)==0)
    lmain <- paste(attr(lzj, "varlabel"), "~", plargs$formula[2])
  pltitle(lmain, sub=plargs$sub, outer.margin=FALSE, cex=plargs$cex.main)
  stamp(sure=FALSE, ploptions=ploptions)
## "plres2x done"
}
## ==========================================================================
plfitpairs <- function(object, ssize=0.02, main=NULL, ...) 
{
  ## Purpose:   pairs plot of fitted values for multinomial regression
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  5 Aug 2004, 10:54
  if (is.null(main)) main <- paste("fitted prob.",object$formula)
  lpr <- object$fitted.values
  lny <- ncol(lpr)
  ly <- object$y
  if(length(ly)==0) stop("!plfitpairs! no response values found")
  ly <- as.numeric(factor(object$y))
##-   if (is.factor(ly)) ly <-  as.numeric(factor())
  if (max(ly)!=lny)
    stop("!plfitpairs! ncol of fitted values != number of levels in y")
##  if (length(pch)<lny) pch <- 1:lny
  lmx <- max(lpr)
  l.panel <- function(x,y,indx,indy,ly,col, ssize) {
    lix <- indx==ly
    liy <- indy==ly
    x[!(lix|liy)] <- NA
    segments(x-ssize*lix,y-ssize*liy,x+ssize*lix,y+ssize*liy,col=col)
    abline(1,-1,lty=3)
  }
  plmatrix(lpr, panel=l.panel, pch=ly, range.=c(0,lmx), main=main, ssize=ssize)
##  "plfitpairs done"
}
## ============================================================================
plmfg <-
  function(mfrow=NULL, mfcol=NULL, mft=NULL, nrow=NULL, ncol=NULL, row=TRUE,
           oma=NULL, mar=NULL, mgp=NULL, ...)
{
## Purpose:    par(mfrow...)
  ## Author: Werner Stahel, 1994 / 2001
  lf.fgt2fg <- function(mft, mfrow, din) {
    if (mfrow==0)
      mfrow <- max(1, ceiling(sqrt(mft*ldin[2]/ldin[1])) )
    lmcol <- ceiling(mft/mfrow)
    c(ceiling(mft/lmcol), lmcol)
  }
  ## number of rows and cols
  lmfg <- if (length(mfrow)==2) mfrow else c(i.def(mfrow, 0), i.def(mfcol,0))
  ldin <- par("din")
  if (length(mft)) lmfg <- lf.fgt2fg(mft, lmfg[1], ldin)
  ## nrow, ncol
  lnfig <- lf.fgt2fg(i.getploption("mfgtotal"), 0, ldin) 
  if (length(nrow)) 
    if(lmfg[1]==0) {
      lmfg[1] <- min(nrow, lnfig[1]+1)
      if (nrow>lmfg[1])
        lmfg[1] <- min(nrow, ceiling(nrow/((nrow-1)%/%lnfig[1]+1)))
      }
  if (length(ncol)) 
    if(lmfg[2]==0) {
      lmfg[2] <- min(ncol, lnfig[2]+1)
      if (ncol>lmfg[2])
        lmfg[2] <- min(ncol, ceiling(ncol/((ncol-1)%/%lnfig[2]+1)))
      }
  lmfg <- pmax(lmfg,1)
  ## mar
  mar <- rep(i.getplopt(mar), length=4)
  if (anyNA(mar)) mar <- ifelse(is.na(mar), par("mar"), mar)
  mgp <- i.getplopt(mgp)
  oma <- if (prod(lmfg)>1) i.getplopt(oma) else i.def(oma, rep(0,4))
  oldpar <- if(row)
              par(mfrow=lmfg, oma=oma, mar=mar, mgp=mgp, ...)
            else par(mfcol=lmfg, oma=oma, mar=mar, mgp=mgp, ...)
  loldo <- ploptions(c("mar","mgp"))
  if (length(mar)+length(mgp)) {
    if (length(mar)) ploptions(mar=mar)
    if (length(mgp)) ploptions(mar=mgp)
  }
  invisible(
    structure(list(mfig = lmfg, mrow = if (row) lmfg, mcol = if(!row) lmfg,
                   mar=mar, mgp=mgp, oma=oma), old=loldo, oldpar=oldpar)
  )
}
## ==========================================================================
stamp <- function(sure=TRUE, outer.margin = NULL,
                  project=getOption("project"), step=getOption("step"),
                  stamp=NULL, line=NULL, ploptions=NULL, ...)
{
## Purpose:   plot date and project information
## -------------------------------------------------------------------------
## Arguments:
##   sure     if F, the function only plots its thing if  getOption("stamp")>0
##   outer    if T, the date is written in the outer margin
##   project  project title
##   step     title of step of data analysis
##   ...      arguments to  mtext , e.g., line=3
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: 13 Aug 96, 09:00
  stamp <- i.getploption("stamp")
  if (length(outer.margin)==0) outer.margin <- par("oma")[4]>0
  t.txt <- date()
  t.txt <- paste(substring(t.txt,5,10),",",substring(t.txt,22,23),"/",
                 substring(t.txt,13,16),sep="")
  if (length(project)>0 && project!="")
    t.txt <- paste(t.txt, project, sep=" | ")
  if (length(step)>0 && step!="")
    t.txt <- paste(t.txt, step, sep=" | ")
  if( sure | stamp==2 | ( stamp==1 & (
    ##     last figure on page :
    { t.tfg <- par("mfg") ; all(t.tfg[1:2]==t.tfg[3:4]) }
    || (is.logical(outer.margin)&&outer.margin) ))  ) {
    lline <-
      i.def(line, ( if(outer.margin) par("oma") else par("mar") )[4] - 1 )
    mtext(t.txt, 4, line=lline, cex = 0.6, adj = 0, outer = outer.margin, ...)
  }
  invisible(t.txt)
}
## ===========================================================================
legendr <- function(x=0.05,y=0.95,legend, ...) {
  lusr <- par("usr")
  lx <- lusr[1] + x*diff(lusr[1:2])
  ly <- lusr[3] + y*diff(lusr[3:4])
  legend(lx,ly,legend, ...)
}

## ======================================================================
colorpale <- function(col=NA, pale=0.3, ...)
{
  pale <- i.def(pale, 0.3)
  lcolna <- is.na(col)
  if (any(lcolna)) {
    col[lcolna] <- palette()[2]
    warning(":colorpale: Argument 'col' is NA. I assume  ", col)
  }
  crgb <- t(col2rgb(col)/255)
  rgb(1-pale*(1-crgb), ...)
}

## ======================================================================
ploptions <-
  function (x=NULL, default=NULL, list=NULL, ploptions = NULL,
            assign=TRUE, ...)
{ ## 
  lnewo <- loldo <-
    if (is.null(ploptions)) {
    if (exists(".ploptions", where=1)) get(".ploptions", pos=1)
    else  ploptionsDefault
    } else ploptions
  largs <- c(list, list(...))
  lop <- c(names(largs))
  ##
  if (!is.null(x)) {  ## asking for options
    if(!is.character(x)) {
      warning(":ploptions: First argument 'x' must be of mode character")
      return(NULL)
    }
    return(if(length(x)==1) loldo[[x]] else loldo[x])
  }
  ## ---
  if (u.notfalse(default) & !is.null(default)) { ## get default values
    if (u.true(default)) default <- "all"
    if (!is.character(default))
      stop("!ploptions! Argument 'default' must be of mode character or logical")
    if (default[1]=="all") return(ploptions(list=ploptionsDefault, assign=assign))
    ## resets all available components
    if (default[1]=="unset")
      return(ploptions(list=ploptionsDefault[names(ploptionsDefault)%nin%
                                             names(loldo)],
                       assign=assign) )
    if (any(default!="")) {
      llopt <- ploptionsDefault[default[default%in%names(ploptionsDefault)]]
      return( ploptions(list=llopt) )
    }
  }
  ## set options
  ## check
  largs <- check.ploption(list=largs)
  if (length(largs))  lnewo[names(largs)] <- largs
  if (assign) assign(".ploptions", lnewo, pos=1)
  ## assignInMyNamespace does not work
  invisible( structure(lnewo, old=loldo[intersect(lop,names(loldo))] ) )
}
## ====================================================================
i.getploption <- function(opt, plo=NULL) {
  ## opt is character, plo list or NULL
  if (is.null(plo))
    plo <- get("ploptions", envir=parent.frame()) ## list in calling fn
  if (is.function(plo)) plo <- NULL
  lopt <- plo[[opt]]  
  if (is.null(lopt)||(!is.function(lopt)&&all(is.na(lopt)))) ## NULL or NA
      lopt <- ploptions(opt)
  else {lopt <- check.ploption(opt, lopt)
    if (length(lopt)) lopt <- lopt[[1]]
  }
  if (length(lopt)==0) lopt <- ploptionsDefault[[opt]]
##  names(lopt) <- opt
  lopt
}
i.getplopt <- function(opt, plo=ploptions) {
  lnam <- as.character(substitute(opt))
  lopt <- opt 
  if (is.function(plo)) plo <- NULL
  if (is.null(lopt)||(is.atomic(lopt)&&all(is.na(lopt))))
    lopt <- plo[[opt]]
  if (is.null(lopt)||(is.atomic(lopt)&&all(is.na(lopt)))) 
    lopt <- ploptions(lnam)
  else unlist(check.ploption(lnam, lopt))   ## check
  if (is.null(lopt)) lopt <- ploptionsDefault[[lnam]]
##  names(lopt) <- opt
  lopt
}
## -----------------------------------------------------------------------
cexSize <- function(n)  min(1.5/log10(n),2)
markextremes <- function(n) ceiling(sqrt(n)/2)/n
smoothpar <- function(n) c(5*n^log10(1/2), 1.5)
smoothxtrim <- function(n, c=1.5) 2^(log(n)/c)/n
## -----------------------------------------------------------------------
check.ploption <- function(optname, value, list=NULL) {
  if (is.null(list)) list <- setNames(list(value), optname)
  lnl <- length(list)
  loptnames <- names(list)
  for (lil in seq_len(lnl)) {
    lnm <- loptnames[lil]
    lvalue <- list[[lnm]]
    lcheck <- ploptionsCheck[[lnm]]
    if (length(lcheck)) {
##    if (is.list(lcheck)) {
      if (!is.list(lcheck[[1]])) lcheck <- list(lcheck)
      lnopt <- length(lcheck)
      lmsg <- rep("", lnopt)
      for (lj in seq_len(lnopt)) {
        lch <- lcheck[[lj]]
        lfn <- get(lch[[1]])
        lmsg[lj] <- lmsgj <-
          switch(paste("v",length(lch),sep=""), v0="", v1=lfn(lvalue),
                 v2=lfn(lvalue, lch[[2]]), v3=lfn(lvalue, lch[[2]], lch[[3]]),
                 "")
        if (lmsgj=="") break
      }
##    }
      if (all(lmsg!="")) {
        warning(":check.ploption: argument '", lnm,
                "' not suitable. It should\n    ",
                paste(lmsg, collapse=" -- or \n  "),
                "\n  instead of (str())\n    ", format(str(lvalue)))
        list[lnm] <- NULL
      }
    }
  }
  list
}
## -----------------------------------------------------------
check.color <- function(x, dummy) {
  if (is.atomic(x) && is.character(x)) {
    lpal <- palette()
    lx <- try(palette(c(x,"black")), silent=TRUE)
    ## palette asks for at least 2 colors
    palette(lpal) ## restore palette
    if (class(lx)=="try-error")
      return("consist of known color names")
    else return("")
  }
  else {
    if(is.atomic(x) && is.numeric(x)) {
      if(all(x>=0 & x<=length(palette()))) return("")
         else return("if numeric, be >=0 and <=length(palette())")
    }
    else {
      if (is.matrix(x)) {
        if (!any(li <- apply(x, 2, function(x) x<0 | x>255))) return("")
        return(
          paste("be a matrix with 3 rows with numbers in [0, 255]", 
                if (length(x)>1) paste("\n  columns ",paste(li, collapse=", "),
                                       " out of range")) )
      }
    }
  }
  "a (vector of) color name(s) or an rgb matrix"
}
##---------
check.numrange <- function(x, range, na.ok=TRUE, length=NA) {
  if (!is.na(length)) {
    if (length > (lnx <- length(x)))
      return(paste("have length at least ",length))
  ##  if (length < lnx) 
  }
  if (na.ok && (is.null(x) || all(is.na(x))) ) return("")
  if (!(is.numeric(x)|is.logical(x))) return("be numeric")
  if ((!na.ok) && any(is.na(x))) return("not contain NAs")
  if (all(is.na(range))) return("")
  range <- ifelse(is.na(range), c(-Inf,Inf), range)
  if (!any(li <- x<range[1]|x>range[2], na.rm=TRUE)) return("")
  paste("be within [",paste(range, collapse=", "),"]",
        if (length(x)>1) paste("\n  violated for element(s) ",
                           paste(which(li), collapse=", ")))
}
check.numvalues <- function(x, values=NA, na.ok=TRUE) {
  if (na.ok && (is.null(x) || all(is.na(x))) ) return("")
  if (!(is.numeric(x)|is.logical(x))) return("be numeric")
  if ((!na.ok) && any(is.na(x))) return("not contain NAs")
  if (all(is.na(values))) return("")
  if (any(li <- (!is.na(x) & (x %nin% values))))
    return(paste("have values in [",
                 paste(values[1:min(5, length(values))], collapse=", "),
                 if (length(values)>5) ",...", "]",
                 if (length(x)>1) paste("\n  violated for elements ",
                                        paste(li, collapse=", "))) )
  ""
}
check.char <- function(x, values, na.ok=TRUE) {
  if (na.ok && (is.null(x) || all(is.na(x))) ) return("")
  if (!is.character(x)) return("be of mode character")
  if ((!na.ok) && any(is.na(x))) return("not contain NAs")
  if (all(is.na(values))) return("")
  if (!any(li <- (!is.na(x) & (x %nin% values)))) return("")
  paste("have values in [",
                 paste(values[1:min(5, length(values))], collapse=", "),
                 if (length(values)>5) ",...", "]",
                 if (length(x)>1) paste("\n  violated for elements ",
                                        paste(li, collapse=", ")))
}
check.logical <- function(x, values, na.ok=TRUE) {
  if (na.ok && (is.null(x) || all(is.na(x))) ) return("")
  if (is.logical(x) | (is.numeric(x)&length(x)==1) )  return("")
  "be of mode logical (or interpretable as such)"
}
check.listnum <- function(x, values=NA, na.ok=TRUE) {
  if (is.list(x)) {
    lchk <- lapply(x, function(xx) check.numvalues(xx, values, na.ok) )
    if (all(lchk=="")) return("")
    return(paste("if a list, all components must be numeric"))
  }
  "be a list"
}
check.function <- function(x, values, na.ok=TRUE) {
  if (is.function(x)) return("")
  if (na.ok && (is.null(x) || all(is.na(x))) ) return("")
  if (is.character(x))  {
    lfn <- try(get(x), silent=TRUE)
    if (class(lfn)=="try-error")
      return(paste("be a function or the name of an existing function.\n   '",
                   lfn, "' is not available.") )
    else return("")
  }
  "be a function or the name of an existing function."
}
## ----------------------------------------------------------------------
cnr <- function(range=NA, na.ok=TRUE, length=NA)
  list("check.numrange", range=range, na.ok=na.ok, length=length)
cnv <- function(values=NA) list("check.numvalues", values=values)
ccv <- function(values=NA) list("check.char", values=values)
ccl <- function() list("check.color", NULL)
clg <- function() list("check.logical", NULL)
cfn <- function() list("check.function", NULL)
cln <- function(values=NA) list("check.listnum", values=values)
## ---------------------------------------------------------------------
c.pchvalues <- 0:180
c.ltyvalues <- 1:6
## ==========================================================================
c.colors <- c("black","red","blue","darkgreen","brown","orange","purple",
              "olivedrab", "burlywood", "violet")
## ----------------------------------------------------------------------
.ploptions <- ploptionsDefault <-
  list(
    colors = c.colors,
    linewidth = c(1,1.3,1.7,1.3,1.2), cex = cexSize,
    ## basic
    pch = 1, cex.pch=1, cex.plab=1,
    lty=1, lwd=1, col=c.colors, lcol=c.colors,
    ## group
    group.pch=2:18, group.col=c.colors[-1], group.lty=2:6,
    group.lcol=c.colors[-1],
    ## variables
    variables.pch=1:18, variables.col=c.colors, variables.lty=1:6,
    variables.lcol=c.colors,
    ## censored
    censored.pch =  c(62, 60, 24, 32, 32, 25, 32, 32),
    ##                 >,  <, Delta, q,q, nabla, q,quadrat 
    censored.size=1.3, censored.pale = 0.3,
    ## frame
    axes = 1:2, mar=c(3.1,3.1,3.1,1.1), oma=c(2.1,2.1,3.1,2.1), mgp=c(2,0.8,0),
    panelsep = 0.5, 
    tickintervals = c(7,3), xlab = "", ylab = "", stamp=1, mfgtotal = 30, 
    innerrange = TRUE, innerrange.factor=4, innerrange.ext=0.1,
    innerrange.function = "robrange", 
    plext=0.05, plextext=0.03,
    markextremes = markextremes, ## is a function...
    ## title (mtext)
    title.cex=c(1.2,1,1), title.cexmin=0.6,
    ## grid
    grid = TRUE, grid.lty = 1, grid.lwd = 1,
    grid.col = "gray85",
    zeroline = TRUE, zeroline.lty = 1, zeroline.lwd = 1,
    zeroline.col = "gray50",
    ## refline
    refline.lty = c(4,6), refline.lwd = c(1,0.7),
    refline.col = "darkgreen",
    ## smoothline
    smoothline.lty = 2, smoothline.lwd = c(2, 0.7),
    smoothline.col = "blue", smoothline.pale = 0.3,
    ## smooth
    smooth = TRUE, 
    smooth.function = "smoothRegr", smooth.par = NA, smooth.iter = 50,
    smooth.minobs = 8,
    ## bars
    bar.lty = 1, bar.lwd = c(2,0.5), bar.col = "burlywood4",
    ## factors
    factor.show = "mbox", jitter = TRUE, jitter.factor = 2,
    ## time axes
    timerangelim = c(year=365*4, month=30*12, day=4),
    ## condquant
    condquant = TRUE, condquant.pale = c(0.5, 0.5), condprob.range = c(0,1),
    ## plotregr
    functionxvalues = 51, smooth.xtrim = smoothxtrim, leveragelim = c(0.99,0.5)
  )
.plargs <- list(ploptions=.ploptions)
  ## makes sure that  .plargs  extists when starting
## -----------------------------------------------------------------------
ploptionsCheck <-
  list(
    colors=ccl(),
    linewidth = cnr(c(0.1,5)), cex = list(cfn(),cnr(c(0.1,5))),
    cex.pch=cnr(c(0.1,5)), cex.plab=cnr(c(0.1,5)),
    ## basic
    pch = cnv(c.pchvalues),
    lty=cnv(c.ltyvalues), lwd=cnr(c(0.1,5)),
    col=ccl(), lcol=ccl(),
    ## group
    group.pch=cnv(c.pchvalues),
    group.col=ccl(), group.lty=cnv(c.ltyvalues),
    group.lcol=ccl(),
    ## variables
    variables.pch=cnv(c.pchvalues), variables.col=ccl(),
    variables.lty=cnv(c.ltyvalues),
    variables.lcol=ccl(),
    ## censored
    censored.pch = cnv(c.pchvalues),
    censored.size=cnr(c(0.1,5)), censored.pale = cnr(c(0,1)),
    ## frame
    mar=cnr(c(0,20)), oma=cnr(c(0,5)), mgp=cnr(c(0,5), na.ok=FALSE, length=3),
    panelsep=cnr(c(0,3)),
    tickintervals = cnr(c(2,20)), stamp=list(clg(),cnr(c(-1,2))),
    mfgtotal = cnr(c(4,100)), 
    innerrange = list(clg(),cnr()), innerrange.factor=cnr(c(0.5,10)),
    innerrange.ext=cnr(c(0,0.5)),
    plext=cnr(c(0,0.5)), plextext=cnr(c(0,0.5)),
    ## title (mtext)
    title.cex=cnr(c(0.1,5)), title.cexmin=cnr(c(0.1,5)),
    ## grid
    grid = list(clg(),cnr(),cln(NA)),
    grid.lty = cnv(c.ltyvalues),
    grid.lwd = cnr(c(0.1,5)),
    grid.col = ccl(),
    ## refline
    refline.lty = cnv(c.ltyvalues), refline.lwd = cnr(c(0.1,5)),
    refline.col = ccl(),
    ## smoothline
    smoothline.lty = cnv(c.ltyvalues), smoothline.lwd = cnr(c(0.1,5)),
    smoothline.col = ccl(), smoothline.pale = cnr(c(0,1)),
    smooth = clg(), 
    smooth.function = cfn(), smooth.minobs = cnr(c(3,20)),
    ## bars
    bar.lty = cnv(c.ltyvalues), bar.lwd = cnr(c(0.1,5)), bar.col = ccl(),
    bar.midpointwidth = cnr(c(0.1,5)),
    ## factors
    factor.show = ccv(c("mbox","jitter","asis","")), jitter = clg(),
    jitter.factor = cnr(c(0.1,5)),
    ## time axes
    ## timerangelim = c(year=365*4, month=30*12, day=4),
    ## condquant
    condquant = clg(), condquant.pale = cnr(c(0,1)),
    condprob.range = cnr(c(0,1)),
    ## plotregr
    functionxvalues = cnr(c(5,500)),
    smooth.xtrim = list(cfn(), cnr(c(0,0.4), na.ok=FALSE)),
      leveragelim = cnr(c(0.1,1))
  )
## ==========================================================================
i.def <- function(arg, value = TRUE, valuetrue = value, valuefalse = FALSE)
{
  rr <- arg
  if (length(arg)==0 ||
      (mode(arg)%in%c("numeric","character","logical","complex")&&
       all(is.na(arg)))
      )  rr <- value
  else {
    if (length(arg)==1 && is.logical(arg))
      rr <- if (arg) valuetrue else valuefalse
  }
  rr
}
## ---------------------------------------
i.col2hex <- function(col) {
  ## convert colors given in any form to rgb
  rgb <- if (is.character(col)) col2rgb(col)
  else
    if (is.matrix(col)&&nrow(col)==3) col
  else
    if (is.numeric(col)&&all(col>=0)) col2rgb(c.colors[col])
  else
    matrix(0,3,length(col))
  lrgb <- rgb/255
  structure(rgb(lrgb[1,],lrgb[2,],lrgb[3,]), names=names(col), rgb=rgb)
}
## ==========================================================================
i.form2char <- function(formula) {
  if (length(formula)==2) paste("~",formula[2])
  else paste(formula[2],"~",formula[3])
}
## ====================================================================
u.allvars <- function(x)
  setNames( lapply(as.list(x),
                   function(lterm) all.vars(as.formula(paste("~",lterm))) ),
           x)
u.varsin2terms <- function(formula) {
  ## which raw variables appear in more than 1 term?
  ltrm <- rownames(attr(terms(formula[c(1,length(formula))]), "factors"))
  lraw <- unlist(u.allvars(ltrm))
  unique(lraw[duplicated(lraw)])
}

i.extendrange <- function(range, ext=0.05)  range + c(-1,1)*ext*diff(range)
clipat <- function(x, range=NULL, clipped=NULL) {
  ## truncate
  if (length(range)==0) return(x)
  lrg <- i.extendrange(range(range), 0.000001) ## make sure limits are not excluded
  li <- which(x>=lrg[1]&x<=lrg[2])
  if (length(clipped)==0) return (x[li])
  x[-li] <- clipped
  x
}


