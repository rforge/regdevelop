## ==================================================================
## Plotting functions, low and high-level
## -------------------------------------------------------------------------
pl.control <-
  function(x=NULL, y=NULL, data = NULL, subset = NULL, transformed = TRUE,
           gensequence = NULL, cex = NULL, 
           psize = NULL, plab = FALSE, pch = NULL, pcol = NULL, cex.pch = NULL,
           smooth.weights = NULL, smooth.weight = NULL,
           markextremes = NULL, smooth = NULL,
           xlab = NULL, ylab = NULL, varlabels = NULL,
           vcol = NULL, vlty = NULL, vpch = NULL, plscale = NULL, log = NULL,
           main = NULL, sub = ":", .subdefault = NULL, mar = NULL, 
           ## needed because it hides  markextremes  otherwise
           ploptions = NULL, .environment. = parent.frame(), assign = TRUE, ... )
  ## get data for plotting, collect and check arguments
  ## do preparations that are common to all plots
  ## --------------------------------------------------------------
{
  ## ---
  lcall <- lcl <- match.call()
  if (length(lcl$smooth.weight)) {
    lcl$smooth.weights <- lcl$smooth.weight
    lcl$smooth.weight <- NULL
  }
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
    if (is.matrix(data)) data <- as.data.frame(data)
    ldtl <- if (length(ltit <- tit(data))) ltit
    if (is.null(ldtl) & is.name(substitute(data))) ldtl <- substitute(data)
  }
  ## --- get variables and data
  lftext <- lform <- NULL
  lformarg <- lynames <- lxnames <- NULL
  if (length(x)==0 & length(y)==0) x <- data
  ## argument x
  if (length(i.def(x, NULL, TRUE, NULL))) {
    if (is.atomic(x)&&is.character(x))  ## names of variables
      lxnames <- x
    else {  ## matrix or data.frame
      if (is.atomic(x)) x <- data.frame(x=x)
      if (is.matrix(x)) x <- as.data.frame(x)
      if (is.data.frame(x)) {
        lxnames <- names(x)
        if (is.null(ldtl))
          ldtl <- if (is.null(ltit <- tit(x)))
                    as.character(attr(x, "dname")) else ltit
        data <- if (is.null(data)) x else {
          if (nrow(x)!=nrow(data))
            stop("!pl.control! arguments 'x' and 'data' have different numbers of rows")
          cbind(x, data)
        }
      } else {  ## formula
        if (is.formula(x)) lform <- x
        else stop("!pl.control! Unsuitable argument 'x'")
      }
    }
    if (is.null(lxnames)) {
      lvars <- getvarnames(lform, data=data, transformed=transformed)
      lxnames <- lvars$xvar
      lynames <- lvars$yvar
    }
    if (is.null(lform))
      lform <-
        as.formula(paste(paste(lynames, collapse="+"), "~",
                         paste(lxnames, collapse="+")))
    environment(lform) <- environment() ## ???
    lftext <- format(lform)
  }
  ## y
  lyform <- NULL
  if (length(y)) {
    if (is.atomic(y)&&is.character(y)) 
      lynames <- y
    else {
      if (is.atomic(y)) y <- cbind(y=y)
      if (is.matrix(y)) y <- as.data.frame(y)
      if (is.data.frame(y)) {
        lynames <- names(y)
        data <- if (is.null(data)) y else {
          if (nrow(y)!=nrow(data))
            stop("!pl.control! arguments 'y' and 'data' have different numbers of rows")
          cbind(y, data)
        }
      }
      else {
        if (is.formula(y)) lyform <- y
        else stop("!pl.control! Unsuitable argument 'y'")
      }
    }
    if (is.null(lynames)) lynames <-
      getvarnames(lyform, data=data, transformed=transformed)$varnames
    if (is.null(lyform))
      lyform <- as.formula(paste("~", paste(lynames, collapse="+")))
    environment(lyform) <- environment() ## ?
    lftext <- paste(shorten(substring(format(lyform),2,100),50), lftext)
    ## fixes maximal length of formula text
  }
  if (length(data)==0)
    stop("!pl.control! arguments 'x', 'y', and 'data' are all empty")
  ## only  x  or  y : plot against sequence
  if ((is.null(lxnames)|is.null(lynames))&&u.notfalse(gensequence)) { 
    data <-
      cbind(".sequence."= 1:NROW(data),data)
    lynames <- c(lxnames, lynames)
    lxnames <- ".sequence."
  }
  ## ---
  lformarg <- list(lxnames,lynames)
  lvarnames <- c(lxnames, lynames)
  ## --- data
  ## ltransformed <- i.def(lcl$transformed, TRUE)
  if (!is.null(attr(data,"terms"))) { ## data is model.frame
    if (!transformed) {
      warning(":pl.control! Raw data not available.",
              "I can only use transformed data.")
      transformed <- TRUE
    }
  }
  ## --- get variables from  lform  and  data
  largs <- c(i.argPldata, "transformed")
  ## variables
  if (length(varlabels)) 
    if (length(names(varlabels))==0) {
      warning(":pl.control: 'varlabels' must have names")
      lcl$varlabels <- varlabels <- NULL
    }
  ## data
  if (length(lvarnames)||length(varlabels) ||
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
    if (length(lpldata)==0)
      stop("!pl.control! No data found")
    ## set varlabels in a special case
    if (length(lxnames)==1 && is.formula(lcx <- lcall$x) && length(lcx)==2)
        attr(lpldata[[lxnames]], "varlabel") <- sub("~","",format(lcx))
    if (length(lynames)==1 && is.formula(lcy <- lcall$y) && length(lcy)==2)
        attr(lpldata[[lynames]], "varlabel") <- sub("~","",format(lcy))
    ##
    if (length(lgroup <- lpldata[[".group."]]))
      if (length(llb <- as.character(lcall$group))<=20)
        attr(lpldata[[".group."]], "varname") <- llb
    ## --- subset
    if (length(lcall$subset)) {
      lsub <- eval(lcall$subset, data, enclos=parent.frame())
      if (sum(lsub)==0)
        stop("!plsubset! No data fulfills the selection criteria")
      lpldata <- plsubset(lpldata, lsub)
    }
    ## --- attributes of variables
    ## plrange
    if(any(c("plrange","xlim","ylim")%in%names(lcall))) 
      lpldata <- i.setplrange(lpldata, ...) 
    ## plscale
    if (length(log)) plscale <-
      ifelse(c(length(grep("x",log)>0),length(grep("y",log)>0)), "log","")
    if (lns <- length(plscale)) {
      if (length(names(plscale))==0)
        if (lns==2)
          plscale <-
            setNames(c(rep(plscale[1],length(lxnames)),
                       rep(plscale[2],length(lynames))), lvarnames)
        else {
          warning(":pl.control: argument 'plscale' not suitable")
          plscale <- NULL
        }
    }
    lvarnames <- attr(lpldata,"variables")
    if (length(lvarnames)) {
      lpldata[,lvarnames] <- 
        genvarattributes(lpldata[,lvarnames, drop=FALSE], vnames = lynames,
                         vcol = vcol, vlty = vlty, vpch = vpch,
                         plscale = plscale, varlabels = varlabels,
                         ploptions=ploptions)
      lnr <- nrow(lpldata)
      lnobs <- lnr-median(sumNA(lpldata[,lvarnames]))
    } else lnr <- 0
  } else {
    lpldata <- NULL
    lvarnames <- NULL
    lnr <- 0
  }
  if (lnr==0) {
    if (length(data)) {
      lnr <- NROW(data)
      lnobs <- lnr-median(sumNA(data))
    } else {
      warning(":pl.control: data not found. I set 'nobs' to 100")
      lnobs <- 100
    }
  }
  ## -------------------------
  lvnm <-
    setdiff(intersect(names(data),
                      c(".pch.",".plab.",".pcol.",".psize.",".smooth.weight.")),
            names(lpldata))
  if (length(lvnm))
    lpldata[,lvnm] <- data[,lvnm, drop=FALSE]
##-     lpldata <- transferAttributes(cbind(lpldata, data[,lvnm, drop=FALSE]),
##-                                   lpldata)
  ## labels anmd plotting character
  ## priorities:  plab , pch , row.names
  lpch <- lpldata$".pch."
  ## default plotting character
  if (length(lpch)) {
    if (is.factor(lpch)) lpch <- as.numeric(lpch)
    if (is.numeric(lpch)) {
      if(length(lpna <- setdiff(lpch,c.pchvalues))) {
        warning(":plcontrol: pch", paste(lpna, collapse=", "), " do not exist")
        lpchunused <- i.def(setdiff(c.pchvalues, c(lpch,0)), 13)
        li <- match(lpch, lpna, nomatch=0)
        lpch[li>0] <- lpchunused[(li[li>0]-1)%%length(lpchunused) +1] ## recycle!
      }
    }  
    if (is.character(lpch) && any(nchar(lpch)>1)) {
      warning(":pl.control: 'pch' must be an integer or a single character.",
              " Use 'plab' to label points with strings")
      lpch <- substr(lpch,1,1)
    }
    if (is.logical(lpch)) lpch <- as.numeric(lpch)
    lpldata$".pch." <- lpch
  }
  lplab <- lpldata$".plab."
  ## --- row.names
  lrown <- substring(row.names(data),1,3)
  if (length(lrown)==0) lrown <- as.character(1:lnr)
  lplabel <- lrown
  lIplab <- length(lplab)>0 && is.logical(lplab) && all(lplab) ## original arg T
  if (lIplab)  lpldata$".plab." <- lplab <- lplabel
  else
    if (length(lplab)) lplabel <- as.character(lplab)
  ## now, lplabel is always useful
  ## --- group
  lgroup <- lpldata$".group."
  if (length(lgroup)) {
    if (is.logical(lgroup)) lpldata[,".group."] <- lgroup+1
    else if (is.factor(lgroup)) lpldata[,".group."] <- factor(lgroup) ## drop levels
  }
  ## pcol
  lpcol <- lpldata$".pcol."
  if (length(lpcol)) {
    if (is.logical(lpcol)) lpldata[,".pcol."] <- lpcol+1
    else if (is.factor(lpcol)) {
      lclr <- i.getploption("color")
      lpldata[,".pcol."] <- lclr[(as.numeric(lpcol)-1)%%length(lclr)+1] ## recycle!
    }
  }
  ## smooth.group
  lsmgrp <- lpldata$"(smooth.group)"
  if (length(lsmgrp)) {
    if (is.logical(lsmgrp)) lpldata[,"(smooth.group)"] <- lsmgrp+1
    else if (is.factor(lsmgrp)) lpldata[,"(smooth.group)"] <- factor(lsmgrp) 
  }
  ## ----------------------------------------------------
  ## more ploptions
  ##
  if (length(ploptions$smooth.col)==1)
    ploptions$smooth.col[2] <-
      colorpale(ploptions$smooth.col, i.getploption("smooth.pale"))
  ## ---
  lmardf <- i.getplopt(mar, ploptions)
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
  ploptions$smooth <- i.getplopt(smooth, ploptions)
  ## i.def(ploptions$smooth, 0, 2, 0) 
  ## --- main
  main <- i.def(main, "", "", "")
  sub <- i.def(sub, NULL, ":", NULL)
  ## --- more arguments
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
i.setplrange <- function(data, plrange=NULL, xlim=NULL, ylim=NULL, ...)
{
  lf.checklim <- function(lim, varnames) {
    if (!is.list(lim)) lim <- list(lim)
    lxv <- attr(data, varnames)
    if (length(lnm <- names(lim))) {
      lim <- lim[intersect(lnm, lxv)]
    } else {
      if (length(lim)==1) lim <- rep(lim, length(lxv))
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
      warning(":plregr/innerrange: unsuitable argument  innerrange ")
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
  function(data, vnames = NULL, vcol = NULL, vlty = NULL, vpch = NULL,
           varlabels = NULL, innerrange = NULL, innerrange.limits = NULL,
           plscale = NULL, zeroline = NULL, replace = FALSE,
           ploptions = NULL)
{
  if (!is.data.frame(data))
    stop("!genvarattributes! 'data' must be a data frame")
  ## ---
  tickintervals <- i.getploption("tickintervals")
  ## variable names
  lnmdata <- colnames(data)
    if (anyNA(lnmdata)) colnames(data) <- lnmdata <-
      ifelse(is.na(lnmdata), paste("V",1:NCOL(data), sep=""), lnmdata)
  ##
  lrown <- row.names(data)
  ## varlabels
  llb <- NULL
  if (length(varlabels))
    if(!is.character(varlabels))
      warning("!genvarattributes! 'varlabels' must be of mode 'character'")
    else {
      if (length(names(varlabels))==0) {
        if (length(varlabels)==ncol(data))
          llb <- setNames(rep(varlabels, length=ncol(data)), colnames(data))
        else warning("!genvarattributes! 'varlabels' must be of length  ",
                     ncol(data), "  or have names")
      }  else llb <- varlabels
    }
  ## jitter
  ljt <- i.getploption("jitter")
  lljt <- length(ljt)
  if (lljt)
    ljt <- if (is.list(ljt)) ljt
            else setNames(rep(ljt, ncol(data)),
                          colnames(data))
  jitter.factor <- i.getploption("jitter.factor")
  ## innerrange
  lir <- i.def(i.getplopt(innerrange), FALSE, TRUE)
  lirg <- setNames(rep(list(lir), ncol(data)), colnames(data)) ## set T/F
  lirglim <- innerrange.limits  
  if (length(lirglim)) {
    if (!is.list(lirglim)) ## convert to list
      lirglim <- if (is.logical(lirglim)) as.list(lirglim) else list(lirglim) 
    if (length(lnm <- names(lirglim))) { ## named list
      lnm <- intersect(lnm, names(lirg))
      lirg[lnm] <- as.list(lirglim)
    } else lirg <- setNames(rep(lirglim, ncol(data)), colnames(data))
  } ## lirg  is now a list. components either T/F or a range.
  lirfactor <- i.getploption("innerrange.factor")
  ## Date
  for (lj in 1:ncol(data))
    if (inherits(ld <- data[,lj], c("Date", "times")) &&
        is.null(attr(ld, "numvalues")))
      data[,lj] <- gendateaxis(setNames(ld,lrown))
  ## plscale
  if (length(plscale)) {
    lnm <- names(plscale)
    if (any(lina <- lnm%nin%lnmdata)) {
      warning(":genvarattributes: variable(s)  ",paste(lnm[lina], collapse=", "),
              "  of argument 'plscale' not in the data frame.",
              "They are left unscaled")
      plscale <- plscale[!lina]
    }
    if (length(plscale)) 
      for (lv in names(plscale))
        data[,lv] <- plscale(data[,lv], plscale[lv])
  }
  ## line color and type
  if (is.null(vnames))
    vnames <- union(union(names(vcol),names(vlty)),names(vpch))
  lny <- length(vnames)
  if (lImulty <- lny>0) { ##
    ldt <- data[,vnames, drop=FALSE]
    lvpch <- 
      i.getvarattribute("vpch", vpch, ldt, ploptions$variables.pch, drop=lny>1)
    ## i.getvarattribute  avoids duplicates
    lvlty <- 
      i.getvarattribute("vlty", vlty, ldt, ploptions$variables.lty, drop=lny>1)
    lvcol <-
      i.getvarattribute("vcol", vcol, ldt, ploptions$variables.col, drop=lny>1)
  }
  ## --- loop
  for (lv in lnmdata) {
    lvv <- data[,lv]
    lcls <- class(lvv)[1]
    attr(lvv, "varname") <- lv
    lnv <- sum(!duplicated(dropNA(lvv)))
    ##
    if (replace || is.null(attr(lvv, "nvalues")) )
      attr(lvv, "nvalues") <- lnv
    ## turn character into factor
    if (lcls=="character") lvv <- factor(lvv)
    ##    if (lv %in% lfacgen)  class(lvv) <- c(class(lvv, "usedAsFactor"))
    if (lImulty) { ## line color and type, pch
      if(lv%in% names(lvcol)) attr(lvv, "vcol") <- lvcol[lv]
      if(lv%in% names(lvlty)) attr(lvv, "vlty") <- lvlty[lv]
      if(lv%in% names(lvpch)) attr(lvv, "vpch") <- lvpch[lv]
    }
    if (inherits(lvv, c("factor", "usedAsFactor", "character"))) {
      data[[lv]] <- lvv <- factor(lvv)
      ## factor
      lat <- seq_along(levels(lvv))
      if (replace || is.null(attr(lvv, "plrange")))
        attr(lvv, "plrange") <- c(0.5, max(lat)+0.5)
      if (replace || is.null(attr(lvv, "ticksat")))
        attr(lvv, "ticksat") <- lat
      if (replace || is.null(attr(lvv, "ticklabels")))
        attr(lvv, "ticklabels") <- levels(lvv)
      ## jitter
      if(replace || is.null(attr(lvv, "numvalues")) && (lij <- ljt[lv])) {
        attr(lvv, "numvalues") <-
          jitter(as.numeric(lvv), factor=jitter.factor,
                 amount=if(is.numeric(lij)) lij else NULL)
        attr(lvv, "plrange") <- c(0.5, length(levels(lvv))+0.5)
      }
      attr(lvv, "zeroline") <- FALSE
    } else { ## continuous variable
      lvvv <- if(inherits(lvv,"Surv")) lvv[,1] else lvv
      lvvv <- if (length(lnv <- attr(lvvv, "numvalues"))) lnv else c(lvvv)
##      class(lvvv) <- setdiff(class(lvv), c("matrix","data.frame"))
      names(lvvv) <- lrown
      lnouter <- c(0,0)
      if (is.numeric(lvvv)) {
##        lrg <- NULL
        lirgv <- lirg[[lv]]
        lplrg <- attr(lvv, "plrange")
        lIplrg <- is.null(lplrg) ## new range
        if (lir && u.notfalse(lirgv) ) {
          ## inner range
          lirga <- attr(lvv, "innerrange")
          lIirg <- is.null(lirga) ## new innerrange
##-           lrg <- if (length(lirgv)==2 && is.numeric(lirgv)) lirgv
##-                  else if (u.notfalse(lirgv))
##-                    plinnerrange(lirgv, lvvv, factor=lirfactor)
          if (!(length(lirgv)==2 && is.numeric(lirgv)))
            lirgv <- plinnerrange(lirgv, lvvv, factor=lirfactor)
          if (replace | lIirg) {
            attr(lvv, "innerrange") <- lirgv
            lpc <- plcoord(lvvv, lirgv, ploptions=ploptions)
            ## attributes: avoid a level...
            lpca <- attributes(lpc)
            if (!lIplrg) lpca$plrange <- lplrg 
            attributes(lvv)[names(lpca)] <- lpca
            attributes(lpc) <- NULL
            attr(lvv, "plcoord") <- lpc
            lnouter <- lpca$nouter
            if (is.null(lnouter)) lnouter <- c(0,0)
            lirgv <- attr(lvv, "innerrange") ## innerrange may have changed
            ## to avoid unnecessary inner bounds when plotting
          }
          ## end inner range
        }
        ## set plrange
        lplrg <- attr(lvv, "plrange")
        if (replace || is.null(lplrg))
          attr(lvv, "plrange") <- lplrg <- 
            i.extendrange(if (length(lirgv)) lirgv else range(lvvv, na.rm=TRUE),
                          i.getploption("plext"))
        ## ticks
        lrg <- if (length(lirgv)==2) lirgv else lplrg
        lnouter <- i.def(attr(lvv, "nouter"), c(0, 0))
        lrg <- ifelse(lnouter>0, lrg, lplrg)
        ltint <- min(max(tickintervals[1],3)+sum(lnouter>0),20)
        lat <- clipat(pretty(lrg, n=ltint, min.n=ltint-2), lrg)
        llabat <-
          if (length(tickintervals)>1)
            clipat(pretty(lat, n=tickintervals[2], min.n=1),lrg)  else lat
        llabat <- clipat(llabat, lat)
        if (length(llabat)<2) {
          llabat <- clipat(pretty(lrg, n=4, min.n=1), lat)  
        } else lat
        if (replace || is.null(attr(lvv, "ticksat"))) {
          attr(lvv, "ticksat") <- lat ## lat[lat>lrg[1]&lat<lrg[2]]
          ## do not replace  ticklabesat  if  ticksat  is set
          if (replace || is.null(attr(lvv, "ticklabelsat")))
            attr(lvv, "ticklabelsat") <- llabat ## lat[lat>lrg[1]&lat<lrg[2]]
        }
        if (replace || is.null(attr(lvv, "zeroline"))) 
          attr(lvv, "zeroline") <- i.getplopt(zeroline)
      }
    } ## end of continuous variable
    attr(lvv,"varlabel") <- if (lv%in%names(llb)) unname(llb[lv]) else
      i.def(attr(lvv, "varlabel"), lv)
    attributes(data[[lv]]) <- attributes(lvv)
  }
  data
}
## =====================================================================
plscale <-
  function(x, plscale = "log10", ticksat = NULL, n = NULL, logscale = NULL)
{
  lfn <- plscale
  if (is.character(plscale)) {
    if (lfn=="") return(x)
    lfnn <- lfn
    lfn <- get(lfnn)
  } else {
    lfnn <- as.character(substitute(plscale))
    if (length(lfnn)>1) lfnn <- "noname"
  }
  ## ---
  lx <- i.def(if (is.null(attr(x,"plscale"))) attr(x, "numvalues"), x)
  ln <- i.def(n, i.getploption("tickintervals"))
  if (is.null(ticksat)) {
    ticksat <- prettyscale(lx, lfnn, n=ln, logscale=logscale)
  } else {
    llab <- attr(ticksat, "ticklabels")
    if (is.null(llab)) llab <- sapply(as.list(ticksat), format)
    ticksat <- structure(lfn(ticksat), ticklabels=llab)
  }
  structure(x, numvalues = lfn(lx),
            ticksat=c(ticksat), ## ticklabelsat=c(ticksat),
            ticklabels=attr(ticksat, "ticklabels"),
            plscale = if (lfnn=="unknown") lfn else lfnn)
}
## =====================================================================
prettyscale <-
  function(x, plscale = "log10", inverse = NULL, range = NULL,
           range.transformed = NULL, n = NULL, logscale = NULL)
{
  ## generate ticks for transformed scale
##-   lf.pretty <- fn(i, x) { ## format: avoid common formatting
##-     lpr <- pretty(x[i+c(-1,1)], n=7)
##-     lpr[which.min(abs(lpr-x[i]))]
##-   }
  lf.format <- function(x) sapply(as.list(x), format)
  lf.cens <- function(x, rg) if (is.null(rg)) x else pmax(pmin(x, rg[2]),rg[1])
  lfn <- if (is.character(plscale)) plscale else
  as.character(substitute(plscale)) 
  if (is.character(plscale)) plscale <- get(plscale)
  ln <- i.def(n, i.getploption("tickintervals"))
  ln1 <- ln[1]+1
  ln2 <- if (length(ln)>=2) ln[2] else 3 ## ceiling(n/2)
  ## --- log according to axTicks
  if (lfn%in%c("log","log10","logst") & is.null(logscale)) {
    lx <- plscale(pmax(x, 0)) ## the transforming function given by the argument
    lrg <- range(lx, na.rm=TRUE)
    if (lfn=="log") lrg <- lrg/log(10)
    if (lrg[2]-lrg[1]>1) { ## use R function
      ltatlab <- 
        axTicks(1, usr=lrg,
                axp=c(10^floor(lrg), max(3-floor((2.5*diff(lrg)+1)/ln1), 1)),
                nintLog=ln1, log=TRUE)
      ltat <- log10(ltatlab)
      if (lfn=="log") ltat <- ltat*log(10)
      return(structure(ltat, ticklabels=lf.format(ltatlab)))
    } ## else - for small range - treat logs by the following
  }
  ## --- back function
  if (lfn=="log")  inverse <- exp
  if (lfn%in%c("log10", "logst"))  inverse <- function(x) 10^x
  if (lfn%in%c("log", "log10", "logst"))  range <- c(0,Inf)
  if (lfn=="sqrt")  {
    inverse <- function(x) x^2
    range <- range.transformed <- c(0, Inf)
  }
  if (lfn=="qnorm") {
    inverse <- pnorm
    range.transformed <- c(0, Inf)
  }
  range <- i.def(range, attr(plscale, "range"))
  range.transformed <-
    i.def(range.transformed, attr(plscale, "range.transformed"))
  if (is.null(inverse)) {
    inverse <- attr(plscale, "inverse")
    if (is.null(inverse))
      stop("!prettyscale! Inverse transformation function missing")
  }
  if (is.character(plscale)) plscale <- get(plscale)
  if (is.character(inverse)) inverse <- get(inverse)
  ## ---
  lrg <- range(plscale(lf.cens(x, range)), na.rm=TRUE) ## transformed scale
  ld <- diff(lrg)/(ln1-1)
  lk <- c(-2,0:(ln1-1),ln1+1)
  lk <- lk+0.01*(lk-mean(lk))/ln1 ## correction helps to break ties
  ltr <- lf.cens(lrg[1]+ld*lk, range.transformed)
  lat <- rep(NA, ln1)
  lrw <- inverse(ltr) ## original scale
  for (li in 1:ln1) {
    lrwx <- lrw[li+c(0,2)]
    lpr <- pretty(lrwx+diff(lrwx)*c(-1,1)*0.5, n=ln2)
    if (length(range)) lpr <- lf.cens(lpr, range)
    ## nearest in transformed scale:
    lat[li] <- lpr[which.min(abs(plscale(lpr)-ltr[li+1]))]
  }
  lat <- unique(lat)
  structure(plscale(lat), ticklabels = lf.format(lat))
}
## =====================================================================
plframe <-
  function(x, y, xlab=NULL, ylab=NULL, ticklabels = TRUE, plextext=NULL,
           axcol = rep(1,4), mar = NULL, 
           plargs = NULL, ploptions = NULL, setpar = TRUE)
{
  ## -------------------
  lf.getcoord <- function(x, lext) {
    lx <- if (is.data.frame(x)) setNames(x[,1], row.names(x)) else x 
    if (length(lx)==0)
      stop("!plframe! unsuitable argument 'x' or 'y'")
    lxx <- i.def(i.def(attr(lx,"numvalues"), attr(lx,"plcoord")), lx)
    ltat <- attr(lx,"ticksat")
    lrg <- attr(lx,"plrange")
    if (is.factor(lxx)) {
      llv <- levels(lx)
      lxx <- c(1, length(llv))
      attr(lx, "ticksat") <- seq_len(length(llv))
      attr(lx, "ticklabels") <- llv
      if (length(lrg)==0) lrg <- lxx+(0.5+lext[1:2])*c(-1,1)
    } else {
      if (length(ltat)==0) 
        attr(lx, "ticksat") <- pretty(lrg, i.getploption("tickintervals")[1])
      if (length(lrg)==0)
        lrg <- i.extendrange(range(lxx, na.rm=TRUE), lext[1:2])
    }
    list(x = lx, xx = lxx, range = lrg)
  }
  ## ---
  if (length(ploptions)==0) ploptions <- plargs$ploptions
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
  ## margins
  lmar <- rep(i.getplopt(mar, ploptions), length=4)
  if (anyNA(lmar)) lmar <- ifelse(is.na(lmar), par("mar"), lmar)
  lmgp <- i.getploption("mgp")
  loldp <- NULL
  if (setpar)
    loldp <- c(par(cex=i.getploption("cex")*par("cex")),
             par(mar=lmar, mgp=lmgp)) 
  ## do not combine with previous line. mar will be in wrong (?) units
  ## on.exit(par(loldp))  ## produces artifact!
  plot(i.extendrange(lrgx, plextext[1:2]), i.extendrange(lrgy, plextext[3:4]),
       xlab = "", ylab = "", type="n", axes=FALSE, xaxs="i", yaxs="i")
  lmfg <- par("mfg")
  ## inner range
  lnouterx <- c(attr(lx, "nouter"),0,0)[1:2]
  lirgx <- c(attr(lx, "innerrange"), lrgx)[1:2]
  abline(v=unique(c(lirgx[lnouterx>0], lrgx)),lty=3)
  lnoutery <- c(attr(ly, "nouter"),0,0)[1:2]
  lirgy <- c(attr(ly, "innerrange"), lrgy)[1:2]
  ## --- grid
  lgrl <- i.getploption("grid")
  ## need attr("ticksat") -> otherwise, generate it!
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
  lzlx <- attr(x, "zeroline")
  lzly <- attr(y, "zeroline")
##-   if (length(lzl <- i.getploption("zeroline"))) {
##-     lzl <- i.def(lzl, 0, TRUE, NA)
##-     if (is.atomic(lzl)) lzlx <- lzly <- lzl
##-     else {
##-       lzlx <- lzl[[1]]
##-       lzly <- lzl[[2]]
##-     }
  if (is.logical(lzlx)) lzlx <- if(lzlx) 0 else NULL
  if (is.logical(lzly)) lzly <- if(lzly) 0 else NULL
  if (length(lzlx)|length(lzly))
    abline(v=lzlx[lzlx>lirgx[1]&lzlx<=lirgx[2]],
           h=lzly[lzly>lirgy[1]&lzly<=lirgy[2]],
           lty=i.getploption("zeroline.lty"), lwd=i.getploption("zeroline.lwd"),
           col=i.getploption("zeroline.col"))
##-  }
  ## bounding boxes
  abline(h=unique(c(ifelse(lnoutery>0, lirgy, lrgy), lrgy)), lty=3)
  abline(v=unique(c(ifelse(lnouterx>0, lirgx, lrgx), lrgx)), lty=3)
  lxrg <- ifelse(lnouterx>0, lirgx, lrgx)
  lyrg <- ifelse(lnoutery>0, lirgy, lrgy)
  lines(lxrg[c(1,2,2,1,1)], lyrg[c(1,1,2,2,1)])
  ## axes
  laxes <- i.getploption("axes")
  if (length(xlab)) attr(lx, "varlabel") <- xlab
  if (length(ylab)) attr(ly, "varlabel") <- ylab
  if (u.notfalse(laxes)) {
    ltint <- i.getploption("tickintervals")[1]
    if (1%in%laxes)
      plaxis(1, lx, lab = lIlab[1] && (lmar[1]>=lmgp[2]+1 | lmfg[1]==lmfg[3]),
             range=lxrg, col=axcol[1], tickintervals=ltint,
             plargs=plargs, setpar=FALSE)
    if (2%in%laxes)
      plaxis(2, ly, lab = lIlab[2] && (lmar[2]>=lmgp[2]+1 | lmfg[2]==1), 
             range=lyrg, col=axcol[2], tickintervals=ltint,
             plargs=plargs, setpar=FALSE)
    if (3%in%laxes)
      plaxis(3, lx, lab = lIlab[3] && (lmar[3]>=lmgp[2]+1 | lmfg[1]==1), 
             range=lxrg, col=axcol[3], tickintervals=ltint,
             plargs=plargs, setpar=FALSE)
    if (4%in%laxes)
      plaxis(4, ly, lab = lIlab[4] && (lmar[4]>=lmgp[2]+1 | lmfg[2]==lmfg[4]),
             range=lyrg, col=axcol[4], tickintervals=ltint,
             plargs=plargs, setpar=FALSE)
  }
  invisible(loldp)
} ## end plframe
## =========================================================================
plaxis <-
  function(side, x=NULL, lab=TRUE, range=NULL, varlabel=NULL, col=1,
           tickintervals=NULL,
           plargs = NULL, ploptions = NULL, setpar=TRUE, ...)
{ ## ------------------------------------------------------------------
  if (length(plargs)==0) plargs <- get(".plargs", globalenv()) ## !!! list in calling fn ?
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  lcex <- par("cex") * if (setpar) i.getploption("cex") else 1
  lmar <- i.getploption("mar")
  loldp <- NULL
  if (setpar) {
    loldp <- c(par(cex=lcex), par(mar=lmar, usr=par("usr")))
    on.exit(par(loldp))
  ##  par(mar=lmar, cex=lcex)
  }
  lmgp <- par("mgp")
  ## ---
  lx <- i.def(i.def(attr(x,"numvalues"), attr(x,"plcoord")), x)
  range <- i.def(range, i.def(attr(x, "innerrange"), attr(x, "plrange")),
              valuefalse = range(lx,na.rm=TRUE))
  lat <- attr(x,"ticksat")
  latsmall <- attr(lat, "small")
  llabat <- attr(x,"ticklabelsat")
  llab <- attr(x,"ticklabels")
  varlabel <- i.def(i.def(varlabel, attr(x,"varlabel"), valuefalse=""),
                    attr(x,"varname") )
  lmfg <- par("mfg")
  col <- i.def(col, 1)
  lIouter <- switch(side, lmfg[1]==lmfg[3], lmfg[2]==1,
                    lmfg[1]==1, lmfg[2]==lmfg[4])
  if(lab && (lmar[side]>=lmgp[1]+1 | lIouter) )
    mtext(varlabel, side=side, line=lmgp[1], xpd=TRUE, col=col,
          cex=lcex*par("cex.lab"))
  ## ticks and tick labels
  if (length(lat)>1) {
    if (length(llabat)) {
      if (!is.numeric(llabat)) {
        warning(":plaxis: 'ticklabelsat' must be numeric")
        llabat <- NULL
      }
    } else llabat <- lat
    if (length(llab)) {
      if (NCOL(llabat)>1) { ## labels placed in the middle of intervals
        lc <- rbind(apply(llabat, 2, clipat, range=range, clipped=NA))
        li <- apply(lc, 1, sumNA) ==2
        llabat <- apply(llabat, 2, clipat, range=range, clipped=range)
        llabat <- apply(rbind(llabat), 1, mean, na.tm=TRUE)
        ##-                    pmin(pmax(llabat[,2],range[1]),range[2]))/2
        llabat[li] <- NA
      }
      if (length(llab)!=length(llabat)) {
        if (length(llab)==1) llab <- rep(llab,length(llabat))
        else {
          warning(":plaxis: 'ticklabel' has wrong length")
          llab <- NULL
        }
      }
    }
    if (length(llab)==0) llab <- if (lab) format(llabat)
    lat <- clipat(lat, range)
    llabat <- clipat(llabat, range, NA)
    li <- is.na(llabat)|is.na(llab)
    llabat <- llabat[!li]
    llab <- llab[!li]
  } else { ## lat not given
    lat <- llabat <-
      clipat(pretty(range, i.getplopt(tickintervals, ploptions)), range)
    llab <- if (lab) format(lat) else FALSE
  }
  if (length(llab)==0 | !lab) llab <- rep("",length(llabat))
  ## axes
  ## tick lengths
  ltcl <- rep(i.getploption("ticklength"))
  if (length(ltcl)<2) ltcl <- rep(c(ltcl,0.5),length=2)*c(1,0)
  if (length(ltcl)<4) ltcl <- c(ltcl, c(ltcl[-(1:2)],0.2)[1]*c(1,-1))
  axis(side, at=lat, labels=rep("",length(lat)), col=col, ...)
  if (ltcl[2]) axis(side, at=lat, labels=rep("",length(lat)),
                    tcl=ltcl[2], col=col, ...)
  if (length(llabat))
    mtext(llab, side, line=i.getploption("mgp")[2], at=llabat, col=col,
          cex=lcex*par("cex.axis"), ...)
  if (length(latsmall)) {
    if (ltcl[3]) axis(side, at=latsmall, labels=rep("",length(latsmall)),
                      tcl=ltcl[3], col=col, ...)
    if (ltcl[4]) axis(side, at=latsmall, labels=rep("",length(latsmall)),
                      tcl=ltcl[4], col=col, ...)
  }
  invisible(loldp)
} ## end plaxis
## ===========================================================================
pltitle <-
  function(main=NULL, sub=NULL, cex=NULL, cexmin=NULL, 
           side=3, line=NULL, adj=NULL, outer.margin=NULL, col="black",
           doc=NULL, show=TRUE, plargs = NULL, ploptions = NULL, ...)
{
  ## Purpose:   title
  ## ----------------------------------------------
  if (length(plargs)==0) plargs <- get(".plargs", globalenv())
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  lparcex <- par("cex")
  lcex <- max(i.def(cex, 1))
  if (lcex!=1) {
    lop <- c(par(cex=lcex*lparcex), par(mar=par("mar"))) ## i.getploption("cex")
    on.exit(par(lop))
  }
  ltadj <- rep(i.def(adj, i.getploption("title.adj"), valuefalse=0.5),
               length=3)
  main <- i.def(main, plargs$main)
  sub <- i.def(sub, plargs$sub, valuefalse=NULL) ## paste(":",plargs$sub,sep="")
  if (is.logical(sub)) sub <- if(sub) ":" else NULL
  sub <- if (length(sub) && substr(sub,1,1)==":")
            paste(substring(sub,2,30), plargs$.subdefault, sep="")
          else as.character(sub)
  if (sub=="") sub <- NULL
  if (length(main)==0 & length(sub)==0) return()
  if (length(main) && all(sub==main)) sub <- NULL
  rr <- c(main=main, sub=sub)
  ## -------------------------
  lf.text <- function(text, cex, cexdef, adj, ...) {
    ## calculate text size and write text
    lcex <-
      i.def(cex, max(cexmin, min(cexdef, lfac/max(nchar(text)))),
            valuefalse = 1 ) 
    lmaxchar <- lfac/lcex
    if (any(li <- nchar(text)>lmaxchar))
      text[li] <- paste(substr(text[li], 1, lmaxchar-3),"...")
    ladj <- i.def(adj, max(minadj,0.5*(lcex>cexmin)), 0.5, minadj)
    mtext(text, side, line, cex = lcex*lparcex, adj=ladj,
          outer = outer.margin, col=col, ...)
    lcex
  }
  ## ---------------------------
  show <- i.def(show, 0.5)
  outer.margin <- i.def(outer.margin, par("oma")[3]>0)
  if (show<=0 || (show<1 && outer.margin && 1!=prod(par("mfg")[1:2])) )
    return(rr)
  ## mar= needed to obtain consistency
  lcex <- cex ##*ifelse(outer.margin, 1, lparcex)
  lcexdef <- rep(i.getploption("title.cex"), length=3) ## *i.getploption("cex")
  title.cexmin <- cexmin
  cexmin <- i.def(cexmin, i.getploption("title.cexmin"), valuefalse=0.1)
  scale <- 2.5 ## scale factor to get width of a character
  minadj <- 0.2
  ##
  lImain <- length(main) && main!=""
  lIsub <- length(sub) && sub%nin%c("",":")
  lIdoc <- length(doc) && doc!=""
  lmarmax <- if (outer.margin) par("oma")[side] else par("mar")[side]
##-   lmarmax <- lmarmax-lcexdef[2-lImain]
##-   llinedef <- sum(lcexdef*c(lImain, lIsub, lIdoc))
##-   line <- min(i.def(line, llinedef), lmarmax)
  line <- lmarmax-0.3-lcexdef[2-lImain] ## -0.3: leave some space above title
  lside24 <- side%in%c(2,4)
  lwid <- if (outer.margin) par("mfg")[4-lside24] else 1
  lfac <- lwid * scale*par("pin")[1+lside24]/(par("cin")[1]*par("cex"))
  ##
  if (lImain) {
    lcex <- lf.text(main, cex=lcex[1], cexdef=lcexdef[1], adj=ltadj[1], ...)
    line <- line-lcex*lparcex
  }
  if (lIsub) {
    lcex <- lf.text(sub, cex=lcex[2], cexdef=lcexdef[2], adj=ltadj[2], ...)
    line <- line-lcex*lparcex
  }
  if (line>=0 && (!is.null(doc)) && doc && length(tit(main)))
    lf.text(tit(main), cex=lcex[3], cexdef=lcexdef[3], adj=ltadj[3], ...)
  invisible(rr)
}
## -----------------------------------------------------------------
plpoints <-
  function(x = NULL, y = NULL, type = "p", plab = NULL, pch = NULL,
           col = NULL, lcol = col, lty = NULL, lwd = NULL, psize = NULL,
           plargs = NULL, ploptions = NULL, setpar=TRUE, ...)
{
  if (length(plargs)==0) plargs <- get(".plargs", globalenv())
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  if (setpar) { ## get  cex, mar,  from ploptions
##-    lusr <- par("usr")
    loldp <- c(par(cex=par("cex") * i.getploption("cex")),
               par(mar=i.getploption("mar"), usr=par("usr")))
    on.exit(par(loldp))
  }
  pldata <- plargs$pldata
  condquant <- i.getploption("condquant")
  ## intro, needed if formulas are used or data is given or ...
  lcall <- match.call()
  if (is.formula(x)|is.formula(y) | any(c("data","pcol")%in%names(lcall))) {
    ## !!! this does not (always) work if plpoints is called by a call
    lcall$assign <- FALSE
    lcall$ploptions <- ploptions
    lplargs <- do.call(pl.control, as.list(lcall[-1]), envir=parent.frame())
    ploptions <- lplargs$ploptions
    plargs$pldata <- pldata <- lplargs$pldata
    x <- pldata[,attr(pldata, "xvar")]
    y <- pldata[,attr(pldata, "yvar")]
  } else {
  ## ---
    if (length(x)==0) x <- pldata[,1]
    if (length(y)==0) y <- pldata[,2]
  }
  ## ---
  if (is.data.frame(x))  x <- x[,1]
  lattrx <- attributes(x)
  lx <- i.def(lattrx$plcoord, i.def(lattrx$numvalues, as.numeric(x)))
  if (is.data.frame(y))  y <- y[,1]
  lattry <- attributes(y)
  ly <- i.def(lattry$plcoord, i.def(lattry$numvalues, as.numeric(y)))
  ## condquant
  lIcq <- condquant>0  ## condquant representation by bars
  lIcqx <- length(lcqx <- lattrx$condquant) >0
  lIcqy <- length(lcqy <- lattry$condquant) >0
  lIcensx <- inherits(x, "Surv")
  lIcensy <- inherits(y, "Surv")
  ##
  lnr <- length(x)
  pldata <- as.data.frame(plargs$pldata) ## 'as.data.frame' only needed for
  ## (not) finding the () columns if plargs is NULL
  psize <- i.def(psize, pldata[[".psize."]])
  lpsize <-
    if (is.null(psize)) 1  else  sqrt(psize/median(psize, na.rm=TRUE))
  if (is.null(plab)) plab <- pldata[[".plab."]]
  if (is.null(pch))
    pch <- i.def(i.def(pldata[[".pch."]], lattry$vpch),
                 i.getploption("pch"), valuefalse="")
  pch <- rep(pch, length=lnr)
  lpcol <- i.def(col, pldata[[".pcol."]])
  if (length(lpcol)) {
    lgrpcol <- i.getploption("group.col")
    if (is.factor(lpcol)) lpcol <- as.numeric(lpcol)
    if (is.numeric(lpcol)) lpcol <- rep(lgrpcol, length=max(lpcol))[lpcol]
    if (is.logical(lpcol)) lpcol <- lgrpcol[lpcol+1]
  } else lpcol <- i.getploption("col")[1]
  pcol <- rep( lpcol, length=lnr)
  lty <- i.def(lty, i.getploption("lty"))
  lwd <- i.def(lwd, i.getploption("lwd")) * i.getploption("linewidth")[lty]
  lcol <- i.getplopt(lcol, ploptions)[1]
  lnobs <- sum(is.finite(lx) & is.finite(ly))
  cex.pch <- i.getploption("cex.pch")
  cex.pch <- if (is.function(cex.pch)) cex.pch(lnobs)
             else i.def(cex.pch, cexSize(lnobs))
  cex.plab <- cex.pch*i.getploption("cex.plab")
  lsplab <- rep(abs(cex.plab*lpsize), length=lnr) 
  lspch <- rep(cex.pch*lpsize, length=lnr)
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
  if (lIcq & (lIcqx|lIcqy)) {
##-     lpale <- rep(c(i.getploption("condquant.pale"), 0.5), length=2)
    lpale <- i.def(i.getploption("condquant.pale"), 0.5, valuefalse=1)[1]
    lcqpch <- i.getploption("condquant.pch")
    lix <- if(lIcqx) lcqx[,"index"]
    liy <- if(lIcqy) lcqy[,"index"]
    lixy <- union(lix,liy)
    if (length(lcqpch)) pch[lixy] <- lcqpch[1]  ## pch
##-     if (length(lixy)==lnr) lpale <- c(1,lpale[1]) ## all observations are cq
##-     pcol[lixy] <- colorpale(pcol[lixy], lpale[1])
    if (length(lixy)<lnr)  pcol[lixy] <- colorpale(pcol[lixy], lpale)
    ##
    if (lIcqx) {
      li <- lix %nin% liy ## why?
      if (any(li)) {
        lsg <- if(length(lrgx <- lattrx$innerrange))
                 plcoord(lcqx[li,2:3], range=lrgx, ploptions=ploptions)
               else lcqx[li,2:3]
        lyy <- ly[lix[li]]
        segments(lsg[,1], lyy, lsg[,2], lyy,
                 col=colorpale(pcol[lix[li]], lpale) )
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
                 col=colorpale(pcol[liy[li]], lpale) )
      }
    }
  }
  ## ---
  lIpl <- (length(plab)>0) && any(!is.na(plab))
  plab <- as.character(plab)
  ## --- plot!
  if (lIpl) {
    text(lx, ly, plab, cex=lsplab*lpsize, col=pcol)
    lipch <- ifelse(is.na(plab), TRUE, plab=="") ## point is not labelled
  } else lipch <- rep(TRUE,length(x))
  if (any(lipch)) {
    if (type%in%c("l","b")) {
      lx[is.na(ly)] <- NA
      ly[is.na(lx)] <- NA
      lxy <- cbind(lx,ly)
      lnox <- attr(x, "nouter")
      if(any(lnox>0)) {
        lrgx <- attr(x, "innerrange")
        if(lnox[1]) 
          lxy <- i.lineout(lxy[,1], lxy[,2], lrgx[1], FALSE)
        if(lnox[2]) 
          lxy <- i.lineout(lxy[,1], lxy[,2], lrgx[2], TRUE)
        lxy[lxy[,1]<lrgx[1]|lxy[,1]>lrgx[2],] <- NA
      }
      lnoy <- attr(y, "nouter")
      if(any(lnoy>0)) {
        lrgy <- attr(y, "innerrange")
        if(lnoy[1]) 
          lxy <- i.lineout(lxy[,2], lxy[,1], lrgy[1], FALSE)[,c(2,1,3)]
        if(lnoy[2]) 
          lxy <- i.lineout(lxy[,2], lxy[,1], lrgy[2], TRUE)[,c(2,1,3)]
        lxy[lxy[,2]<lrgy[1]|lxy[,2]>lrgy[2],] <- NA
      }
      points(lxy[,1], lxy[,2], type=type, pch=NA, col=attr(y, "vcol"), ...)
      if (type=="l") return()
      type <- "p"
    }
    points(lx[lipch], ly[lipch], type=type, pch=pch[lipch], cex=lspch[lipch],
           col=pcol[lipch])
  }
  NULL
} ## end plpoints
pllines <- function(x=NULL, y=NULL, type="l", ...) {
  plpoints(x=x, y=y, type=type, ...)
}
## ---------------------------------------------------------
i.lineout <- function(z, zz, lim, upper) {
  ## determines indices before which a point on the inner range limit is needed
  ## and the "other" coordinate
  lout <- if(upper) z>lim  else  z<lim
  li0 <- which(c(lout,FALSE) & !c(TRUE,lout))
  lzout <- zz
  lz <- z[li0-1]
  lr <- (lim-lz)/(z[li0]-lz)
  lzz0 <- zz[li0-1]+lr*(zz[li0]-zz[li0-1])
  li1 <- which(c(FALSE,lout) &! c(lout, TRUE))
  lz <- z[li1-1]
  lr <- (lim-lz)/(z[li1]-lz)
  lzz1 <- zz[li1-1]+lr*(zz[li1]-zz[li1-1])
  lii <- order(c(1:length(z), li0-0.6, li1-0.4))
  zz[c(li0,li1-1)] <- NA
  cbind(z=c(z,rep(lim,length(c(li0,li1))))[lii], zz=c(zz,lzz0,lzz1)[lii],
        zout=c(ifelse(lout, lzout, NA), lzz0, lzz1)[lii])
}
## -----------------------------------------------------------------
plbars <-
  function(x = NULL, y = NULL, midpointwidth = NULL,
           plargs = NULL, ploptions = NULL, setpar = TRUE)
{
  if (length(plargs)==0) plargs <- get(".plargs", globalenv())
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  if (setpar) { ## get  cex, mar,  from ploptions
    loldp <- par(mar=i.getploption("mar"))
    on.exit(par(loldp))
  }
  if (length(x)==0) x <- plargs$pldata[,1]
  if (length(y)==0) y <- plargs$pldata[,2]
  if (NCOL(x)==1) x <- c(x)
  if (NCOL(y)==1) y <- c(y)
  if (!(is.null(dim(x))&NCOL(y)>=3 | is.null(dim(y))&NCOL(x)>=3)) {
    warning(":plbars: unsuitable arguments 'x', 'y' -> no bars.")
    return()
  }
  lnobs <- sum(apply(is.finite(cbind(x,y)),1,all))
  lcol <- i.getploption("bar.col")
  llty <- rep(i.getploption("bar.lty"), length=2)
  llwd <- rep(i.getploption("bar.lwd"), length=2)
  lmpw <- i.def(midpointwidth, i.getploption("bar.midpointwidth")) * 0.3 *
    diff(par("usr"))[c(1,3)]/ lnobs
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
  function(x=NULL, y=NULL, markextremes=NULL, plabel=NULL,
           plargs = NULL, ploptions = NULL)
{ ## ---------------------------------
  lf.pmark <- function(mprop, x) {
    lrk <- (rank(x, na.last="keep")-0.5)/sum(is.finite(x))
    lrk<mprop[1] | lrk>1-mprop[2]
  }
  ##
  if (length(plargs)==0) plargs <- get(".plargs", globalenv())
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  plabel <- i.def(i.def(plabel, plargs$plabel), .plargs$plabel)
  if (is.null(plabel)) {
    warning(":plmark: no plabels found")
    return(NULL)
  }
  if (length(x)==0) x <- plargs$pldata[,1]
  lx <- if (is.data.frame(x)) x[,1] else x 
  ly <- if (is.data.frame(y)) y[,1] else y 
  lnobs <-
    if (is.null(ly)) sum(is.finite(lx)) else sum(is.finite(lx)&is.finite(ly))
  lmxdef <- ceiling(sqrt(lnobs)/2)/lnobs
  lmx <- i.getplopt(markextremes)
  if (is.function(lmx)) lmx <- markextremes(lnobs)
  if (is.atomic(lmx)) lmxx <- lmxy <- i.def(lmx, lmxdef)
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
} ## end plmark
## -----------------------------------------------------------------
plsmooth <-
  function(x, y = NULL, ysec = NULL, band=NULL, power = NULL,
           plargs = NULL, ploptions=NULL)
{
  if (length(plargs)==0) plargs <- get(".plargs", globalenv())
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  lIsm <- i.getploption("smooth")
  band <- i.def(band, i.getploption("smooth.band"))
  lsm <- NULL
  if (lIsm) {
    power <- i.def(power, 1,1,1)
    band <- i.def(i.def(band, lIsm>=2), FALSE, TRUE, FALSE)
    ly <- if(length(ysec)==0) y else {if (length(y)==0) ysec else cbind(y,ysec)}
    lsm <- gensmooth(x, ly, band=band, power=power,
                     plargs=plargs, ploptions=ploptions)
    ysecsm <- NULL
      if (length(ysec)) { ## separate  y  and  ysec
      lny <- NCOL(y) * (length(y)>0) ## ... NCOL(NULL) is 1 !
      ysecsm <- lsm[["y"]]
      if (lny==0) lsm$y <- NULL else {
        lsm$y <- ysecsm[,1:lny]
        ysecsm <- ysecsm[,-(1:lny)]
      }
    }
    plsmoothline(lsm, x, ysec = ysecsm, plargs=plargs, ploptions=ploptions)
  }
  invisible(lsm)
}
## --------------------------------------------------------------
plsmoothline <-
  function(smoothline, x, y = NULL, ysec = NULL,
           plargs = NULL, ploptions = NULL, ...)
{
  if (length(plargs)==0) plargs <- get(".plargs", globalenv())
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  if (!is.list(smoothline)) {
    warning(":plsmoothline: 'smoothline' is not suitable. No smooth line")
    return()
  }
  lxtrim <- i.getploption("smooth.xtrim")
  lrgx <- attr(x,"innerrange")
  x <- i.def(smoothline[["x"]], x)
  if (length(x)==0) stop("!plsmoothline! no x values")
  y <- i.def(smoothline[["y"]], y)
  if (length(y)) y <- as.matrix(y)
  if (length(ysec)) ysec <- as.matrix(ysec)
  if (length(y)==0 & length(ysec)==0) {
    warning(":plsmoothline: neither 'y' nor 'ysec' given. No smooth line")
    return()
  }
  lIy <- length(y)>0
  lIysec <- length(ysec)>0
  lgrp <- as.numeric(smoothline$group)
  if (lInogrp <- length(lgrp)==0 || length(unique(dropNA(lgrp)))<=1) {
    lgrp <- rep(1, NROW(x))
    lngrp <- 1
    lcol <- i.getploption("smooth.col")
    llty <- i.getploption("smooth.lty")
  } else {
    lngrp <- max(lgrp)
    lcol <- i.getploption("group.col")
    llty <- rep(i.getploption("group.lty"), length=lngrp)
  }
  lbd <- smoothline$yband
  lIband <- length(lbd)>0
  ## check if ordered
  lio <- order(lgrp, x) 
  if (any(lio!=1:length(lio))) {
    x <- x[lio]
    if (lIy) y <- y[lio,, drop=FALSE]
    if (lIysec) ysec <- ysec[lio,, drop=FALSE]
    lgrp <- lgrp[lio]
    lngrp <- max(lgrp)
    if (lIband) {
      lbd <- lbd[lio]
      smoothline$ybandindex <- smoothline$ybandindex[lio]
    }
  }
  x <- clipat(x, lrgx, clipped=NA) 
  ## lny <- ncol(ly)
  ## may be a 2-vector or  a matrix of 2 rows
  llwd <- rep(c(i.getploption("smooth.lwd"), 0.7), length=2)
  llwid <- i.getploption("linewidth")
  lpale <- i.getploption("smooth.pale")
  for (lgr in seq_len(lngrp)) {
    lig <- which(lgrp==lgr)
    lxg <- x[lig]
    if (1< (lng <- length(lig))) {
      lndr <- i.def(smoothline$xtrim, 0) * 
        round(lng * if(is.function(lxtrim)) lxtrim(lng)
                    else i.def(lxtrim, 0, smoothxtrim(lng), 0) )
      if (lndr) lxg[- ((lndr+1):(lng-lndr)) ] <- NA
      lcl <- lcol[min(lgr,length(lcol))]
      llt <- llty[lgr]
      llw <- llwd[1]*llwid[llt]
      if (lIysec) 
        matlines(lxg, ysec[lig,], lty=llt, lwd=llwd[2]*llw, 
                 col = colorpale(lcl, lpale), ...)
      if (lIy) matlines(lxg, y[lig,], lty=llt, lwd=llw, col=lcl, ...)  ## xxx
      if (lIband) {
        li <- smoothline$ybandindex[lig] ## separate upper and lower smooth
        if (any(li))
          lines(lxg[li], lbd[lig[li]], lty=llt, lwd=llw/2, col = lcl, ...) 
        if (any(!li))
          lines(lxg[!li], lbd[lig[!li]], lty=llty, lwd=llw/2, col = lcl, ...) 
      }
    }
  }
  NULL
}
## -----------------------------------------------------------------
plrefline <-
  function(refline, x=NULL, innerrange=NULL, y=NULL,
           cutrange = c(x=TRUE, y=FALSE),
           plargs = NULL, ploptions = NULL, ...)
{
  ## draws a reference line (with extended range) and
  ##   band given by reflineyw (only inner range) if requested
  lf.irna <- function(x, rg) {
    x[x<rg[1]|x>rg[2]] <- NA
    x
  }
  if (length(plargs)==0) plargs <- get(".plargs", globalenv())
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  if (is.character(x)) x <- plargs$pldata[,x]
  if (is.character(y)) y <- plargs$pldata[,y]
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
    if (is.null(x)) x <- plargs$pldata[, i.def(attr(plargs$pldata,"xvar"),1)[1]]
    if (is.null(y)) y <- plargs$pldata[, i.def(attr(plargs$pldata,"yvar"),1)[2]]
    lIxir <- any(attr(x,"nouter")>0)
    lxir <- attr(x,"innerrange")
    lIyir <- any(attr(y,"nouter")>0)
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
  NULL
}
## =========================================================================
plcoord <-
  function(x, range=NULL, innerrange.factor=NULL, innerrange.ext=NULL,
           plext=NULL, ploptions=NULL)
{
  ## Purpose:    values for plot with limited "inner" plot range
  lx <- structure(as.numeric(x), dim=dim(x))
  ldtrg <- range(lx, na.rm=TRUE)
  lirfunc <- i.getploption("innerrange.function")
  if (is.character(lirfunc)) lirfunc <- get(lirfunc)
  innerrange.factor <- i.getplopt(innerrange.factor, ploptions)
  innerrange.ext <- i.getplopt(innerrange.ext)
  plext <- i.getplopt(plext)
  if (length(dropNA(range))>0) {
    lrg <- range(range, na.rm=TRUE)
    if (lrg[1]>ldtrg[2] | lrg[2]<ldtrg[1]) { ## ranges do not overlap
      warning(":plcoord: inadequate range. Not used")
      lrg <- NULL
    } 
  } else lrg <- NULL
  if (length(lrg)==0) lrg <- lirfunc(lx, fac=innerrange.factor)
  ## inner range must not extend beyond data
  if (length(lrg)==0) lrg <- ldtrg
  if (diff(lrg)==0) lrg <- c(-1,1)*lrg
  lrgext <- i.extendrange(lrg, plext)
  if (ldtrg[1]>=lrgext[1]) lrg[1] <- ldtrg[1]
  if (ldtrg[2]<=lrgext[2]) lrg[2] <- ldtrg[2]
  ## --------- transformation of data into plcoord
  rr <- pmax(pmin(lx,lrg[2]),lrg[1])
  lxd <- lx-rr
  lnouter <- c(sum(lxd<0,na.rm=TRUE),sum(lxd>0,na.rm=TRUE))
  if (sum(lnouter)>0) rr <- rr+lxd/(1+abs(lxd)/(diff(lrg)*innerrange.ext))
  ## if data fits into extended inner range, then avoid inner range
  lrg <- ifelse(lnouter>0, lrg, ldtrg)
  ## ---------
  ## extend range to plotting range
  lplrg <- i.extendrange(lrg, ifelse(lnouter>0, innerrange.ext, plext))
  ## extend inner range if there are no modified points
  lrg <- ifelse(lnouter>0, lrg, lplrg) ## adjust to the needs of plotting 
  attr(rr,"innerrange") <- lrg
  attr(rr,"innerrange.ext") <- innerrange.ext
    ## needed for transforming further quantities
  attr(rr,"nouter") <- lnouter
  attr(rr,"plrange") <- lplrg
  class(rr) <- class(lx)
  rr
}
## -------------------------------------------------------------------
## =====================================================================
gendate <-
  function (date=NULL, year=2000, month=1, day=1, hour=0, min=0, sec=0,
            data=NULL, format="y-m-d", origin=NULL)
{
  ## generate time -> chron
  ## -------------------------------------------------------------
  lcall <- match.call()
  ltimeargs <- c("date", "year", "month", "day", "hour", "min", "sec")
  ldt <- list(date=NA, year=2000, month=1, day=1, hour=0, min=0, sec=0)
  origin <- i.def(origin, c(month=1, day=1, year=i.getploption("date.origin")) )
  lnm <- NULL
  ## --- arguments
  largs <- setdiff(names(lcall)[-1], c("data", "format","origin"))
  data <- as.data.frame(data)
  if (length(largs)==0)    stop("!gendate! no date and time variables specified")
  inp <-
    parse(text = paste("list(", paste(as.list(lcall[largs]), collapse = ","),")"),
                 keep.source = FALSE)
  vars <- try(eval(inp, data, enclos=parent.frame()), silent=TRUE)
  if (class(vars)=="try-error") {
    lvnmiss <- setdiff(largs, names(data))
    stop(sub("object", "!gendate! variable (or object)",
               attr(vars, "condition")$message),
         if (length(lvnmiss)>1)
           paste(". \n    All of ",
                 paste(lvnmiss, collapse=", "), "may be unavailable.")
         )
  }
  names(vars) <- largs
  ldt[largs] <- vars
  ldt <- as.data.frame(ldt)
  ##
  lln <- unique(lapply(ldt, length))
  if (length(lln[lln>1])>1)
    stop("!gendate! differing lengths of arguments")
  ## month
  if (is.factor(ldt$month)) ldt$month <- as.character(ldt$month)
  if (is.character(ldt$month)) {
    lmnum <- match(ldt$month, c.months, nomatch=0) ## distiguish NA and nommatch
    if (any(lmnum==0, na.rm=TRUE)) lmnum <- match(ldt$month, c.mon, nomatch=0)
    if (any(lmnum==0, na.rm=TRUE))
      stop("!gendate! inadequate argument 'month'")
    ldt$month <- lmnum
  }
  ## --- argument 'date'
  date <- ldt$date
  if (length(date)&&!all(is.na(date))) { ## get  timeargs  from  date
    if (inherits(date, "Date"))
      date <- chron::dates(format(date), format="y-m-d", origin=origin)
    else {
      if (inherits(date,"factor")) date <- as.character(date)
      if (is.character(date))
        date <- chron::dates(date, format=c(format, "h:m:s"), origin=origin)
    }
    if (!inherits(date, "times")) 
      stop("!gendate! argument 'date' not suitable.")
##    if (length(setdiff(largs, "date"))==0) return(chron::chron(date))
  } else 
    date <- chron::dates(
      paste(ldt$year, ldt$month, floor(as.numeric(ldt$day)), sep="-"),
      format="y-m-d", origin=origin)
  ## ---
##-   if (inherits(date, "chron")) ## date already has 'times'
##-     return(structure(date, format=i.getploption("date.format")) )
  ## --- hours, min, sec
  ## convert too large numbers
  if ("sec"%in%largs && any(li <- ldt$sec>=60, na.rm=TRUE)) {
    li <- which(li)
    ldt$min[li] <- ldt$min[li]+ldt$sec[li]%/%60
    ldt$sec[li] <- ldt$sec[li]%%60
  }
  if ("min"%in%largs && any(li <- ldt$min>=60, na.rm=TRUE)) {
    li <- which(li)
    ldt$hour[li] <- ldt$hour[li]+ldt$min[li]%/%60
    ldt$min[li] <- ldt$min[li]%%60
  }
  if ("hour"%in%largs && any(li <- ldt$hour>=24, na.rm=TRUE)) {
    li <- which(li)
    date[li] <- date[li]+ ldt$hour[li]%/%24
    ldt$hour[li] <- ldt$hour[li]%%24
  }
  ## convert decimals to lower units
  lf.dec <- function(x)  round(x,2)-floor(x+0.005) ## not general enough
  if ("day"%in%largs && any(0!= (ldec <- lf.dec(ldt$day)), na.rm=TRUE)) {
    ldt$day <- floor(ldt$day) 
    if ("hour"%nin%largs) ldt$hour <- round(ldec*24 + ldt$hour - 0.000005, 5)
  }
  if (any(0!= (ldec <- lf.dec(ldt$hour)), na.rm=TRUE)) {
    ldt$hour <- floor(ldt$hour+0.005)
    if ("min"%nin%largs) ldt$min <- round(ldec*60 + ldt$min - 0.0005)
  }
  if (any(0.001<abs(ldec <- lf.dec(ldt$min)), na.rm=TRUE)) {
    ldt$min <- floor(ldt$min+0.005)
    if ("sec"%nin%largs) ldt$sec <- round(ldec*60 - 0.5) ## + ldt$sec  must be zero
  }
  ##
  date <-
    chron::chron(date, paste(ldt$hour, ldt$min, ldt$sec, sep=":"),
                 format=c(dates=format, times="h:m:s"), origin=attr(date,"origin"))
  if (length(lnm)) names(date) <- lnm
  structure(date, format=i.getploption("date.format"))
} ## end gendate
## =====================================================================
gendateaxis <-
  function(date=NULL, year=2000, month=1, day=1, hour=0, min=0, sec=0,
           data=NULL, format="y-m-d", origin=NULL, ploptions=NULL)
{
  ## generate time axis.
  ## resulting tick labels may exceed data range by quite a bit
  ## -------------------------------------------------------------
  lf.2char <-
    function(x)
      substring(ifelse(x<10, paste("0",x,sep=""), as.character(x)),1,2)
  lf.seq <- function(x) seq(min(x),max(x))
  lf.tickat <-
    function(tickunit, tickint, llev, llvlg, ystart, mstart, lnlev) {
      ## generate ticks in  tickint [tickunit]  intervals
      if (tickunit=="y") return(ystart)
      llunit <- match(tickunit, names(lnlev))
      llev[[llunit]] <- ltatu <-
        seq(llunit%in%2:3, lnlev[llunit], tickint) ## m and d start at 1
##        seq(0, lnlev[llunit], tickint) + (llunit%in%2:3) ## m and d start at 1
      ## keep llev component of highest category if higher than tickunit,
      ## generate levels of lower ones for the whole span PLUS 1 for the end
      lv1 <- llev[[llvlg]]
      if (llvlg < llunit-1) {
        for (ll in (llvlg+1):llunit-1) {
          llev[[ll]] <- llv <- seq(ll%in%2:3, lnlev[ll]) ## m and d start at 1
          lv1 <- c(outer(llv, lv1*100, "+"))
        }
      }
      if (llvlg<llunit)
        ltatu <- c(outer(ltatu, lv1*100, "+")) ## information for getting label
      ## ---
      if (tickunit=="m")  ltat <- mstart[seq(1, length(mstart), tickint)]
      else { ## mstart contains start of month for whole year(s). select those needed
        ltat <- c(outer(llev[["d"]]-1, mstart[llev[["m"]]], "+")) ## day starts with 1
        lmd <- c(outer(llev[["d"]], 100*(llev[["m"]]), "+"))
        ltatu[lmd%in%c(limpossible, limpossible+1)] <- NA 
        liat <- lmd%nin% (limpossible+1) ## +1: keep tick at end of last day month
        ltat <- ltat[liat]
        ltatu <- ltatu[liat]
        if (tickunit!="d") { ## !!! mit llev arbeiten! oder 1:length(lv1)
          ltat <- c(outer(llev[["h"]]/24, ltat, "+"))
          if (tickunit!="h") {
            ltat <- c(outer(llev[["M"]]/(24*60), ltat, "+"))
            if (tickunit=="M") ltat <- ltat[seq(0, length(ltat), tickint)]
            else {
              ltat <- c(outer(llev[["s"]]/(24*60), ltat, "+"))
              if (tickunit=="s") ltat <- ltat[seq(0, length(ltat), tickint)]
              else  stop("!gendateaxis/lf.tickat! bug")
            }}}
      }
      if (length(ltat)!=length(ltatu)) 
        warning("debug lf.ticksat: length(ltat)!=length(ltatu)")
      ## avoid duplicated ticks. select second label
      li <- length(ltat)+1-which(duplicated(rev(ltat)))
      if (length(li)) {
        ltat <- ltat[-li]
        ltatu <- ltatu[-li]
      }
      structure(ltat, at.inunit=ltatu)
    } ## end lf.tickat
  lnlev <- c(y=100, m=12, d=31, h=24, M=60, s=60)
  limpossible <- c(229, 230, 231, 431, 631, 931, 1131) ## non-existing days
  ## ---------------------------------------
  ## --- prepare
  lcall <- match.call()
  ltickint <- i.getploption("tickintervals")
  ldtk <- i.getploption("date.ticks")
  ## gendate
  lcall$ploptions <- NULL
  lcall[[1]] <- quote(gendate)
  mode(lcall) <- "call"
  date <- eval(lcall, sys.parent())
  ## ---------
  lattr <- attributes(date)
  lvlab <- lattr$varlabel ## for later use
  lorigin <- if (length(lor <- lattr$origin))
    julian(as.Date(paste(lor[c("year","month","day")], collapse="-"))) else 0  
  if (length(lattr)) { ## avoid piling up of attributes
    li <-
      setdiff(names(lattr),
              c("numvalues", "ticksat", "ticklabelsat", "ticklabels", "units") )  
    attributes(date) <- lattr[li]
  }
  lndays <- diff(range(date))
  lidtk <- which(ldtk$limit<lndays)[1]
  ## -----------------------------------------------------------------
  ldate <- data.frame(chron::month.day.year(date+lorigin)[c(3,1,2)],
    hour = chron::hours(date), min = chron::minutes(date),
    sec = chron::seconds(date))
  ## avoid new day or month or year caused by a midnight obs.
  if (diff(range(ldate$hour))!=0)
    ldate <- ldate[!(ldate$day==1&ldate$hour==0&ldate$min==0),]
  llev <- lapply(ldate, lf.seq)
  names(llev) <- names(lnlev)
  ## which units vary?
  lvr <- sapply(llev, length) >1
  llvlg <- which(lvr)[1]
  llsm <- length(lvr)+1-which(rev(lvr))[1]
  ##
  lIym <- length(c(unique(ldate$year), unique(ldate$mon)))>2
  lyr <- lf.seq(ldate$year)
  lystart <- julian(as.Date(paste(lyr,"1-1", sep="-")))-lorigin
  lmn <- 1:12 ## if (llvlg>=3) lf.seq(ldate$mon)+1 else 1:12
  lyy <- if (length(lyr)==1) lyr else rep(lyr, each=12)
  lmstart <- julian(as.Date(paste(lyy,lmn,"1", sep="-")))-lorigin
  ##
  ltatsmall <-
    c(lf.tickat(ldtk$smallunit[lidtk], ldtk$smallint[lidtk],
                llev=llev, llvlg=llvlg, ystart=lystart, mstart=lmstart,
                lnlev=lnlev) ) ## ( attr "at.inunit" may be too short. not used here)
  ltatbig <-
    c(lf.tickat(ldtk$bigunit[lidtk], ldtk$bigint[lidtk],
                llev=llev, llvlg=llvlg, ystart=lystart, mstart=lmstart,
                lnlev=lnlev) )
  llunit <- ldtk$labelunit[lidtk]
  llint <- ldtk$labelint[lidtk]
  ltatlabel <-
    lf.tickat(llunit, llint,
              llev=llev, llvlg=llvlg, ystart=lystart, mstart=lmstart,
              lnlev=lnlev)
  ## --- ticklabels
  ltatu <- attr(ltatlabel, "at.inunit")
  lls <- any(ltatu>=100, na.rm=TRUE)
  ll <- list(u1 = ltatu%%100, u2 = if (lls) ltatu%/%100, sep=lls)
  lf.nabl <- function(x) {x[x=="NA"] <- "";x}
  llab <- 
    switch(llunit,
           y = as.character(lyr),
           m = paste(c.mon[ll$u1], ll$u2, sep = if(ll$sep) "." else ""),
           d = lf.nabl(paste(ll$u1, c.mon[ll$u2], sep=if(ll$sep) " " else "")),
           h = paste(ll$u2, ll$u1, sep = if(ll$sep) "|" else ""),
           M = paste(ll$u2, ll$u1, sep = if(ll$sep) ":" else ""),
           s = paste(ll$u2, ll$u1, sep = if(ll$sep) ":" else ""),
           {
             warning(":gendateaxis: labels went wrong. check ploptions(\"date.ticks\"")
             NULL
           }
           )
  llab[is.na(ltatu)] <- ""
##-   if (llunit=="d") {
##-   }
  if (length(lvlab)==0 || (length(l <- attr(lvlab, "setbyuser"))&&!l) ) {
    llv1 <- c("A",names(llev))[llvlg] ## level that does not vary
    lvlab <- structure(
      switch(llv1,
             A = "year",
             y = as.character(llev[[1]]),
             m = paste(c.months[llev[[2]]], llev[[1]]), 
             d = paste(llev[[3]], c.months[llev[[2]]]),
             h = "minute",
             M = "second",
             "time"),
      setbyuser = FALSE)
  }
  ## drop labels outside range
  lina <- is.na(clipat(ltatlabel, i.extendrange(unclass(range(date))), clipped=NA))
  if (any(lina)) {
    ltatlabel <- ltatlabel[!lina]
    llab <- llab[!lina]
    ltatu <- ltatu[!lina]
  }
  ## ticklabelsat
  llu <- match(llunit, c("y","m","d","h","M","s"), nomatch=0)
  if (llu<llsm) ## turn into interval
    ltatlabel <- outer(ltatlabel, c(0, c(365, 30, 1, 0,0)[llu]), "+")
  ## 1/24, 1/(24*60)
  ##      at.inunit = attr(ltatlabel, "at.inunit") )
  if (llunit=="d" & ldtk[lidtk,"labelint"]!=1) { ## drop mark at day 31
    li30 <- ltatu%%100==30
    li30[length(li30)] <- FALSE ## do not drop at end of scale
    llab[li30] <- ""
  }
  attr(ltatlabel, "at.inunit") <- NULL
  structure(date, numvalues=c(unclass(date)),
            ticksat=structure(ltatbig, small=ltatsmall),
            ticklabelsat=ltatlabel, ticklabels=llab, varlabel=lvlab)
}
## ===========================================================================
i.pchcens <-
  function(plargs, condquant)
    ##  Delta, nabla, >, <, quadrat : pch= c(24, 25, 62, 60, 32)
{
  if (is.null(condquant) | !is.null(lpc <- plargs$pldata$".pch."))
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
      warning(":plregr/pllimits: unsuitable argument  pllim ")
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
  function(x, subset = NULL, omit = NULL, select = NULL, drop = FALSE,
           keeprange = FALSE)
{
  ## subset, adjusting attributes
  ## adaptation of  subset.data.frame 
  if (!inherits(x, "data.frame"))
    stop("!plsubset! 'x' must inherit from data.frame")
  e <- if (lIsubset <- length(match.call()$subset))
         substitute(subset) else substitute(omit)
  r <- if (lIsubset & length(e)==0) rep_len(TRUE, nrow(x))
       else  dropNA(eval(e, x, parent.frame()))
  if (is.logical(r)) r <- which(r)
  if (is.character(r)) {
    ir <- match(r, row.names(x), nomatch=0)
    if (all(ir==0))
      stop("!plsubset! 'subset' is character, but contains no row.names of 'x'")
    if (any(ir==0))
      warning(":plsubset: 'subset' contains strings that are no row.names:\n  ",
              if (sum(ir==0)<=10) paste(r[ir==0], collapse=", ")
              else  paste(paste(r[ir==0][1:10], collapse=", "), "...")
              )
    r <- ir[ir!=0]
  }
  if (!lIsubset) r <- setdiff(1:nrow(x), r) ## omit has been used
  rr <- transferAttributes( x[r,, drop=FALSE], x)
  if (length(select)) {
    nl <- as.list(seq_along(x))
    names(nl) <- names(x)
    vars <- eval(substitute(select), nl, parent.frame())
    rr <- rr[vars,, drop=FALSE]
  }
  if (nrow(rr)==0) {
    warning(":plsubset: Empty subset")
    return(rr)
  }
  ## ----- attributes
  lf.rgratio <- function(x, r) {
    lrx <- diff(range(x, na.rm=TRUE))
    if (lrx>0) diff(range(x[r], na.rm=TRUE))/lrx else 1
  }
  for (lj in seq_len(ncol(rr))) {
    lxj <- x[,lj]
    if (length(lattr <- attributes(lxj))) {
      lrgratio <-
        if (is.numeric(lxj)) lf.rgratio(lxj, r) else 1
      if ("numvalues"%in%names(lattr)) {
        lnv <- lattr$numvalues
        lnvr <- lnv[r]
        lrgratio <- lf.rgratio(lnv, r)
        lattr$numvalues <- lnvr
      }
      if ("plcoord"%in%names(lattr)) {
        lpc <- lattr$plcoord
        lpcr <- lpc[r]
        lrgratio <- lf.rgratio(lpc, r)
        lattr$plcoord <- lpcr
      }
      if (!is.finite(lrgratio)) {
        warning(":plsubset: something wrong with attr(...,'numvalues') ",
                "or attr(...,'plcoord')")
        lrgratio <- 1
      }
      ## !!! ticksat if(...)
      if ((!keeprange) & lrgratio<i.getploption("subset.rgratio")) {
        lats <- c("ticksat","ticklabelsat","ticklabels") ## drop these attributes
        lattr <- lattr[setdiff(names(lattr), lats)]
        ## !!! getvarattributes(...)
        ## date variable
        if (inherits(lxj, c("Date", "times"))) {
          attr(rr[,lj], "numvalues") <- lattr$numvalues
          lattr <- attributes(gendateaxis(rr[,lj]))
          lattr[names(lattr)] <- lattr
        }
        ## scaled variable
        if (length(lpls <- lattr$plscale)) {
          lplsc <- plscale(lxj[r], lpls)
          lattr[lats] <- attributes(lplsc)[lats]
        }
      }
      attributes(rr[,lj]) <- lattr
      ## resp. genvarattributes
    }
  }
  rr[,,drop=drop]
}

## ========================================================================
plyx <-
  function(x=NULL, y=NULL, group=NULL, data=NULL, type="p", panel=NULL,
           xlab=NULL, ylab=NULL, 
           markextremes=0, rescale=TRUE, mar=NULL,
           plargs = NULL, ploptions = NULL, ...)
{
  ## --- intro
  lcall <- match.call() ## match.call modifies argument names -> sys.call
##-   lcnm <- i.def(names(lcall), names(lcl))
##-   names(lcall) <- ifelse(lcnm=="", names(lcl), lcnm)
  if (length(plargs)==0) {
    lcall$markextremes <- markextremes
    ldtnm <- substitute(data)
    ldtnm <- if (is.name(ldtnm))
               as.character(ldtnm) else if (length(ldtnm)) format(ldtnm)
    lcall$.subdefault <- if (length(ldtnm)<30) ldtnm else ""
    lcall[1] <- list(quote(pl.control))
    lcall <- as.call(lcall)
    plargs <- eval(lcall, envir=parent.frame())
  }
  ##
  if (length(ploptions)==0) ploptions <- plargs$ploptions
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
        if(inherits(lyj, "Surv"))
          survival::Surv(lplc, lyj[,2], type=attr(lyj, "type")) else lplc ,
        lyj )
    }
  }
  lny <- ncol(ly)
  ## style elements
  lIsmooth <- i.getploption("smooth")
  lIfirst <- TRUE
  lplab <- pldata[[".plab."]]
  lpch <- pldata[[".pch."]]
  lIpch <- length(lpch)>0
  lpcol <- pldata[[".pcol."]]
  lpexgrp <- i.getploption("cex.pch") ## par("cex")
  lcex <- i.getploption("cex") ## *par("cex")
  lmar <- i.getplopt(mar, ploptions)
  loma <- i.def(i.getploption("oma"), par("oma"), valuefalse=rep(0,4)) ## ??? 
  lmgp <- i.getploption("mgp")
  panel <- i.getplopt(panel)
  if (is.character(panel)) panel <- get(panel, globalenv())
  if (!is.function(panel)) {
    warning(":plyx: 'panel' not found. Using 'plgraphics::plpanel'")
    panel <- plgraphics::plpanel
  }
  lpsep <- i.getploption("panelsep")
  lyaxcol <- 1
  ## group
  lgroup <- pldata[[".group."]]
  lIgrp <- length(lgroup)>0
  if(lIgrp)   {
  ##  lgrpname <- i.def(attr(lgroup, "varname"), "group")
    lgrplab <- as.character(unique(lgroup))
##-     if (length(lpcol)==0) {
##-       lcl <- i.getploption("group.col")
##-       lpcol <- lcl[1+(as.numeric(lgroup)-1)%/%length(lcl)]
##-     }
  } else  lgroup <- rep(1, nrow(ly))
  if (is.factor(lgroup)) lgroup <- as.numeric(lgroup)
  lgrp <- unique(lgroup)
  lngrp <- length(lgrp)
  ## ranges
  lf.rg <- function(y)
    i.def(i.def(attr(y, "innerrange"), attr(y,"plrange")),
          range(y, na.rm=TRUE))
  lrgy <-
    sapply(ly, lf.rg)
  ## inner plotting range
  lnoutery <- as.matrix(sapply(ly, function(x) c(attr(x,"nouter"),0,0)[1:2]>0 ))
  lIinner <- apply(lnoutery, 1, any)
  if (lny>1) { ## need to examine all y's and possibly reset plrange for ly[,1]
    if (rescale==0) {
      lyy <- genvarattributes(as.data.frame(c(as.matrix(ly))),
                              ploptions=ploptions)
      lattr <- attributes(lyy[,1])[c("innerrange", "plrange", "ticksat",
                                     "ticklabelsat")]
      ly <- setvarattributes(ly, setNames(rep(list(lattr), lny), lynm))
      lrgy <- matrix(lf.rg(ly[,1]), 2, lny)
      for (lj in 1:lny) 
        attr(ly[,lj], "plcoord") <- plcoord(ly[,lj], range=lrgy)
    } else {
      if (rescale<0) { ## do not adjust tick marks etc
        lrgyy <- c(min(lrgy[1,]),max(lrgy[2,]))
        for (lj in 1:lny) 
          attr(ly[,lj], "plcoord") <- plcoord(ly[,lj], range=lrgyy)
      } 
      else  ## rescale >0
        lrgyy <- lrgy[,1]
      attr(ly[,1],"plrange") <-     ## extend the plrange of ly[,1] if needed
        lrgyy + diff(lrgyy)*c(-1,1)*
          ifelse(lIinner, i.getploption("innerrange.ext"), 0) ## ploptions$plext
    } ## !!! welche attr sollen wirklich gesetzt werden?
    attr(ly[,1], "nouter") <- lIinner ## plframe is called for ly[,1], needs attr
  } ## end lny>1
  ly1 <- ly1g <- ly[,1]
  lrgy1 <- lrgy[,1]
  ## mark extremes
  lmark <- i.getplopt(markextremes, ploptions)
  if (is.function(lmark)) lmark <- lmark(lnobs)
  lmk <- unlist(lmark)
  lImark <- length(lmk)>0 && any(ifelse(is.na(lmk),TRUE,lmk>0))
  lImark <- is.na(lImark)||lImark
  lymark <- if (lImark & lny==1) ly  ## cannot mark extremes if  lny>1
  ##
  if (lny>1 & is.null(mar)) lmar[4] <- lmar[2]
  ploptions$mar <- lmar
  lsmcol <- i.getploption("smooth.col")
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
      lmfig <- plmframes(lnx, lngrp, oma=loma)
      loldp <- attr(lmfig, "oldpar")
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
  lattrc <- c("innerrange","nouter")
  ## --- plot
  for (ipgr in 1:lnpgr) {
    lr <- (ipgr-1)*lnr ## start at row index  lr
    for (lj in 1:min(lnr, lnx-lr)) {
      lxj <- lx[,lr+lj]
      ## lxjattr <- 
      lxjg <- lxjv <-
        i.def(attr(lxj, "plcoord"), i.def(attr(lxj, "numvalues"), c(lxj)))
      attributes(lxjg)[lattrc] <- attributes(lxj)[lattrc] ## lxjattr
    if (lImark) {
      lplab <- plmark(lxj, y=lymark, markextremes=lmark, plabel=plargs$plabel)
      pldata[[".plab."]] <- lplab
      plargs$pldata <- pldata
    }
    lpchg <- lpch
    ## groups
    for (ipgc in 1:lnpgc) { ## columns on page
      lc <- (ipgc-1)*lnc
      for (lig in (lc+1):lnc) { ## groups
        lg <- lgrp[lig]
        if(lny>1)  {
          lyaxcol <- lsmcol <- lpcol <- attr(ly1, "vcol")
          plargs$ploptions$smooth.col <- lsmcol
          pldata[".pcol."] <- lpcol
          if(!lIpch) 
            lpchg <-
              if (lny>1) attr(ly1, "pch") else i.getploption("pch")[1]
        }
        ## frame
        ## attr(ly1, "nouter") <- lIinner 
        lop <- plframe(lxj, ly1, xlab=xlab, ylab=ylab,
                       axcol=c(NA,lyaxcol,NA,NA), ploptions=ploptions,
                       setpar=lIfirst)
        if (lIfirst) {
          loldp <- c(lop, loldp)
          on.exit(par(loldp[!duplicated(names(loldp))]))
        }
        lIfirst <- FALSE
        if (1==prod(par("mfg")[1:2])) {
          pltitle(plargs=plargs)
          stamp(sure=FALSE, ploptions=ploptions)
        }
        lrgold <- lrgy1
        lyg <- ly
        ## group
        if (lIgrp) {
          if (lmar[3]>=1 | par("mfg")[1]==1) 
            pltitle(main=NULL, sub=lgrplab[lg], ##paste(lgrpname, lgrplab[lg], sep=": "),
                    outer.margin=FALSE, xpd=TRUE,
                    line=min(1.2*lcex,par("mar")[3]-0.5*lcex), adj=0.03, 
                    cex=rep(lcex,2))  
          li <- which(lgroup==lgrp[lg])
          lxjg <- lxjv[li]
          lyg <- plsubset(ly, li, keeprange=TRUE) 
          ly1g <- lyg[,1] ## transferAttributes(lyg[,1], lyg)
          plargs$pldata <- pldata[li,]
          if (length(lpch)==lnr) lpchg <- lpch[li]
          if (length(lplab)==lnr) lplabg <- lplab[li]
        }
        if (lny>1) 
          ploptions$col <-ploptions$lcol <-
            plargs$pldata$".pcol." <- attr(ly1g,"vcol")
        panel(lxjg, ly1g, type=type, plargs=plargs, ploptions=ploptions)
        ## multiple y
        lusr <- par("usr")
        if (lny>1) {
          for (lj in 2:lny) {
            lyjg <- lyg[,lj]
            lrgj <- lrgy[,lj]
            lpcol <- attr(lyjg, "vcol") ##  the color must reflect the variable
            if (!lIpch) lpchg <- attr(lyjg, "pch")
            if (rescale>0) {
              lusr[3:4] <- lrgj[1] + diff(lrgj)/diff(lrgold)*(lusr[3:4]-lrgold[1])
              par(usr=lusr)
            }
            if (lIsmooth) {
              plargs$ploptions$smooth.col <- lpcol
              plsmooth(lxjg, lyjg, plargs=plargs)
            }
            plpoints(lxjg, lyjg, type=type, plab=lplab, pch=lpchg, col=lpcol,
                      plargs=plargs, setpar=FALSE)
            ##!!!linecolor
            lrgold <- lrgj
            if (lj==2) {
              lmfg <- par("mfg")
              plaxis(4, lyjg, lab=lmar[4]>=lmgp[2]+1 | lmfg[2]==lmfg[4],
                     range=lrgj, col=lpcol,
                     tickintervals=i.getploption("tickintervals"), setpar=FALSE)
            }
          }
        }
      }
    }
  }
  }
  invisible(plargs)
} ## end plyx
## ==========================================================================
setvarattributes <-
  function(data, attributes = NULL, list = NULL, ...)
{
  data <- as.data.frame(data)
  lnmdata <- names(data)
  list <- c(list, list(...))
  if (length(list) && !is.list(list))
      warning(":setvarattributes: argument 'list' must be a list")
  list <- c(list, list(...))
  if (length(list)) {
    lnames <- names(list)
    for (lnm in lnames) {
      lls <- as.list(list[[lnm]])
      llnm <- names(lls)
      if (any(linm <- llnm%nin%lnmdata))
        stop("!setvarattributes! names of  ",lnm,"  not in ",
             "names of 'data': ", paste(llnm[linm], collapse=", "))
      for (lnmv in llnm)
        attr(data[[lnmv]], lnm) <- lls[[lnmv]]
    }
  } else if (length(attributes)==0) warning(":setvarattributes: no attributes")
  if (length(attributes)) {
    if (!is.list(attributes))
      stop("!setvarattributes! argument 'attributes' must be a list")
    if (is.null(names(attributes))) {
      if (length(attributes)!=ncol(data))
        stop("!setvarattributes! argument 'attributes' must have names ",
             "or be of appropriate length")
      names(attributes) <- names(data)
    }
    lnames <- names(attributes)
    if (any(linm <- lnames%nin%names(data)))
      stop("!setvarattributes! names of argument 'attributes' not in ",
           "names of 'data': ", paste(lnames[linm], collapse=", "))
    for (lnm in lnames) {
      lattr <- attributes(data[[lnm]])
      lattr[names(attributes[[lnm]])] <- attributes[[lnm]]
      attributes(data[[lnm]]) <- lattr
      if (lnm=="innerrange") {
        ## call plcoord !!!
      }
    }
  } 
  invisible(data)
}
## ==========================================================================
plregr.control <-
  function(x, data = NULL, xvar = TRUE, transformed = FALSE,
           ## generate variables used for plotting
           weights = NULL, stdresid = TRUE, mar = NULL, 
           ## specify some specific aspects / contents of plots
           glm.restype = "working", condquant = TRUE, smresid = TRUE, 
           partial.resid = TRUE, cookdistlines = NULL,
           leveragelimit = NULL, condprobRange = NULL,
           testlevel = 0.05,
           ## smooth and refline
           refline = TRUE, 
           smooth = 2, ## smoothPar=NA, smoothIter=NULL,
           smooth.sim=NULL,
           xlabs = NULL, reslabs = NULL, markextremes = NULL, 
           ## multiple frames
           mf = TRUE, mfcol = FALSE, multnrow = 0, multncol = 0, ## multmar = NULL,
           oma = NULL, assign = TRUE, ... )
{ ## get data for plotting, collect and check arguments
  ## do preparations that are common to plregr and plresx
  ## --------------------------------------------------------------
  lcall <- match.call()
  lnaaction <- x$na.action
  if (!is.null(lnaaction)) class(lnaaction) <- "exclude"
  x$na.action <- lnaaction
  ## --- family
  lfam <- c(x$distrname, x$family$family, "")[1]
  if (lfam=="" & inherits(x, "lm")) lfam <- "gaussian"
  lfamgauss <- lfam%in%c("gaussian","Gaussian")
  lfamcount <- (lfam%in%c("binomial","poisson")&&length(unique(x$y))<=2) | 
    inherits(x,"polr") ## lfam=="multinomial" | 
  ## --- na.action: always get full data
  ## residuals first because they fix the number of observations
  lres <- residuals(x)
  lcq <- i.getplopt(condquant) 
  rtype <- i.def(i.def(if(lfamcount & lcq) "condquant", glm.restype), "working")
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
                      function(x) all(x==dropNA(x)[1], na.rm=TRUE ) ) )
  if (lres0)
    stop("!plregr/plresx! all residuals are equal -> no residual plots")
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
  llevlim <- i.getplopt(leveragelimit) ## should be , ploptions
  lstdres <- attr(lres, "stdresiduals")
  if (stdresid) { ## not needed for plresx
    if (length(lstdres)==0) { 
      lstdres <- stdresiduals(x, residuals=lres, leveragelimit=llevlim[2])
      llev <- attr(lstdres, "leverage")
      attr(lres, "stdresratio") <- attr(lstdres, "stdresratio")
      attributes(lstdres)[c("leverage","stdresratio","stddev")] <- NULL
      attr(lres, "stdresiduals") <- lstdres
    } ## else lstdres <- lres
    else llev <- naresid(lnaaction, llev)
    if (length(llev)==0) llev <- leverage(x)
  }
  ## --- xvar
  ## formula(x) inappropriate for merMod
  lform <- if (inherits(x, "regrMer")) x$formula else formula(x)  
    ## ... if it inherits from  lm
  lmodvdupl <- u.varsin2terms(lform[-2])
  if (u.notfalse(xvar)) {
    lxf <- if (transformed) lform else {
      lx <- all.vars(lform[-2])
      if (inherits(x, "nls")) lx <- setdiff(lx, names(coef(x)))
      u.asformula(lx)
    }
    if (!(is.null(xvar)|u.true(xvar))) {
      if (is.character(xvar)) {
        lxvarf <- u.asformula(setdiff(xvar,"."))
        lxf <- if ("."%in%xvar) update(lxf, lxvarf) else lxvarf
      }
      else {
        if (is.formula(xvar)) lxf <- update(lxf, xvar)
        else
          stop("!plregr.control! Inadequate argument 'xvar'")
      }
    }
    lxvar <- getvarnames(lxf, data=data, transformed=transformed)$xvar ## was TRUE
    lxvraw <- u.allvars(lxvar)
  } else lxvar <- NULL
  ## residual names
  lyexpr <- deparse(lform[[2]])
  lynm <- if (nchar(lyexpr)>10) "Y" else lyexpr ##
  lrn <- paste("res_", if (lmres>1) colnames(lres) else lynm, sep="")
  lresname <-
    gsub("\\(",".", gsub("\\)",".", gsub("\\*",".", gsub("/",".", gsub("\\^",".",
                                                                     lrn)))))
  names(lres) <- lresname
  ## --- prepare  pl.control
  lcall <- as.list(match.call())[-1]
  ladrop <- c("xvar", "glm.restype", "smresid", "partial.resid",
              "cookdistlines", "leveragelimit", "smooth", "smooth.sim",
              "refline", "testlevel", "xlabs", "reslabs",
              "mf", "mfcol", "multnrow", "multncol", "multmar", "oma")
  lcall$y <- lres
  lcall <- c(as.list(quote(pl.control)),
             as.list(lcall[setdiff(names(lcall),ladrop)]))
  ## extract  y  from model
  ldata <- x$model
  if (length(ldata)&&length(lnaaction))
    ldata <- i.naresid.exclude(lnaaction, ldata)
  ly <- i.def(x[["y"]], if (length(ldata)) ldata[,1, drop=FALSE]) 
  if (length(lxvar)&u.notfalse(lxvar)) {
    lcall$x <- lxvar
    lcall$transformed <- transformed
    ## lcall$.subdefault <- i.form2char(lform) ## transfer to main
    ## --- data argument
    if (length(lxvar)||any(names(lcall)%in%i.argPldata)) {
      if (length(ldata)==0||!(transformed & all(lxvar%in%names(ldata)))) {
        ldata <- data
        if (length(ldata)==0) {
          if (length(x$allvars))
            ldata <- i.naresid.exclude(lnaaction, x$allvars)
          else ldata <- getvariables(x$call$formula, eval(x$call$data), transformed=FALSE) ##eval(x$call$data)
          ##-     } else {
##-       if (length(lav <- x$allvars)&&NROW(data)==NROW(lav)) 
##-        ldata <- cbind(data, lav[colnames(lav)%nin%colnames(data)])
##-       else stop("!plregr.control! incompatible data and x$allvars")
    ##-     }
        }
        if (length(ldata)==0) {
          ldata <- i.naresid.exclude(lnaaction, x$model)
        } else
          if ("subset"%in%names(x$call)) ldata <- ldata[row.names(lres),]
        if (length(ldata)==0)  stop("!plregr.control! No data found")
      }
  ##
      if (lnres!=nrow(ldata)) {
        if (class(lnaaction)=="omit") ldata <- ldata[-lnaaction,]
        ##    if (lnr!=nrow(ldata)) ldata <- x$model ## needs at least a warning!
        if (lnres!=nrow(ldata))
          stop("!plregr.control! nrow of residuals and data do not agree.")
      }
    }
    if (any(lxvj <- lxvar%nin%names(ldata)))
      stop("!plregr/plresx! variable(s) ",
           paste(lxvar[lxvj], collapse=", "), " not found")
  } else {
    ldata <- lres ## needed to get number of observations
    lcall$x <- FALSE
  } ## needed to get  nobs  in pl.control
  lcall$data <- ldata
  lftext <- i.form2char(lform)
  lcall$.subdefault <- lftext
  lcall$assign <- FALSE
  mode(lcall) <- "call"
  ## -------------------------------------
  plargs <- eval(lcall, parent.frame())  ## pl.control
  ## -------------------------------------
  plargs$ploptions$smooth <- i.getplopt(smooth)
  ploptions <- plargs$ploptions
  plargs$main <- i.def(ploptions$main, i.form2char(lform))
  lpldata <- plargs$pldata 
  ## margins for multivariate regression
  lmmar <- rep(i.getploption("panelsep"), length=4)
  lmmar[2] <- i.getplopt(mar)[2] +0.5
##-   if (length(mar)) ## !!! marmult
##-     lmmar <- ifelse(is.na(lmr <- rep(mar, length=4)), lmmar, lmr)
##-     lmmar[1] <- i.def(i.getploption("mar")[1], 3)
##-     lmmar[3] <- i.def(i.getploption("mar")[3], 0.5) ##!!!
  xvar <- attr(lpldata, "xvar")
  ## -------------------------------------------------------
  ## attributes for residuals
  lres <- genvarattributes(lres, varlabels = lresname, ploptions=ploptions)
  if (lmres>1) { ## multivariate
    if (is.null(lcn <- colnames(lres))) lcn <- 1:ncol(lres)
    colnames(lres) <- lcn # paste("res", lcn, sep=".")
  }
  ## mark extreme  stdres
  lmxdef <- markextremes(lnobs)
  if (is.atomic(markextremes)) {
    markextremes <- i.def(markextremes, NA)
    if (anyNA(markextremes)) markextremes <- lmxdef
    markextremes <- 
      list("(res)"=markextremes, "(fit)"=0, "(lev)"=c(0,max(markextremes)) )
  }
  lmxres <- i.def(markextremes$"(res)", lmxdef)
  lresplab <-
    if (lmxres>0 & stdresid)
      apply(lstdres, 2,
            function(x) plmark(x, markextremes = lmxres, plabel=plargs$plabel) )
    else NULL
  ## --- smoothWeights, used for smooth calculation: get from  x  if needed
  lsmwgt <- lpldata[["(smoothWeights)"]]  ## possibly only logical
  lIsmweights <- is.logical(lsmwgt)&&all(lsmwgt) ## weights explicitly required
  lInosmweights <- is.logical(lsmwgt)&&!any(lsmwgt) ## weights explicitly denied
  if (lIsmweights | (length(lsmwgt)&&all(is.na(lsmwgt))))
      lsmwgt <- naresid(lnaaction, x$weights)
  lIsmwgt <-
    length(lsmwgt)>1 && any(lsmwgt!=dropNA(lsmwgt)[1],na.rm=TRUE) 
  if (lIsmweights&!lIsmwgt)
    warning(":plregr/plresx: no weights found for smooth calculation.")
  lpldata[["(smoothWeights)"]] <- lsmweights <-
    if (lIsmwgt)  lsmwgt / mean(lsmwgt, na.rm=TRUE) else NULL
  ## --- psize, used as sizes of plotting characters: same as weights
  lpsize <- lpldata[[".psize."]]  ## possibly only logical
  lIpsize <- is.logical(lpsize)&&all(lpsize) ## psize expl. required
  lInopsize <- is.logical(lpsize)&&!any(lpsize) ## psize expl. denied
  if (is.null(lpsize) | lIpsize | (length(lpsize)&&all(is.na(lpsize))))
      lpsize <- if (length(lsmwgt)) lsmwgt else naresid(lnaaction, x$weights)
  lIpsz <-
    length(lpsize)>1 && any(lpsize!=dropNA(lpsize)[1],na.rm=TRUE)
  if (lIpsize&!lIpsz)
    warning(":plregr/plresx: no plot sizes found.")
  lpldata[[".psize."]] <- 
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
  if (lnsims>0 & !inherits(x, c("lm","glm","lmrob","nls"))) {
    warning(":plregr/simresiduals: ",
            "I can simulate only for 'lm', 'lmrob', 'nls' and 'glm' objects")
    lnsims <- 0
  }
  lsimres <- NULL
  if (lnsims>0) {
    lsimres <- if(ploptions("debug"))
               simresiduals(x, lnsims, glm.restype=glm.restype)  else
      try(simresiduals(x, lnsims, glm.restype=glm.restype), silent=TRUE)
    if (class(lsimres)=="try-error") {
      warning(":plregr/simresiduals: simresiduals did not work. ",
              "No simulated smooths")
      lsimres <- NULL
      lnsims <- 0
    }
  }
  ## --- multiple frames
  mf <- i.def(mf, NULL, TRUE, NULL)
  oma <- i.def(oma, NULL, valuefalse=0)
  if (length(oma)==2) oma <- c(0,0,oma)
  ## --- more arguments
  reslabs <- i.def(reslabs, NULL, NULL, NULL)
  smresid <- i.def(smresid, TRUE)
  if (lmres>1) smresid <- FALSE ## !!! muss noch gemacht werden
  ## -----------------
  partial.resid <- i.def(partial.resid, TRUE)
  cookdistlines <- i.getplopt(cookdistlines)
  plargs$ploptions$refline <- i.getplopt(refline)
  testlevel <- i.getplopt(testlevel)
  if (testlevel<=0 | testlevel>=1)
    stop("!plregr.control! invalid test level")
  ## ------------------------------------------------------------
  ## result of plregr.control
  plargs$pldata <- lpldata
  plargs$smooth <- i.getplopt(smooth)
  rr <-
    c(plargs,
      list(      
      residuals = lres, rescol = lmres,
      response = ly, 
      ## weights = x$weights,
      leverage = llev,
      ## linear.predictor = x$linear.predictor,
      resmahal = x$resmahal, 
      ## resmahal = lresmahal,
      simres = lsimres, ## simstdres <- lsimstdres,
      yexpr = lyexpr, resname = lresname, ## stdresname = lstdresname,
      ##    absresname = labsresname, fitname = lfitname,
      family = lfam, famgauss = lfamgauss, famcount = lfamcount,
      formula = lform, na.action = lnaaction, 
      sigma=lsigma, df.residual = ldfres, df.model = ldfmod, 
      ## -- return arguments
      glm.restype = glm.restype, smresid = smresid,
      partial.resid = partial.resid, cookdistlines = cookdistlines,
      leveragelimit = llevlim, ## condprobRange = condprobRange,
      testlevel=testlevel,
      ## smooth and refline
      smooth.sim = lnsims,
      refline = refline,
      resplab = lresplab,
      mf=mf, multnrow = multnrow, multncol = multncol, multmar = lmmar,
      oma=oma #,
      ) )
  if (u.notfalse(assign)) assign(".plargs", rr, pos=1)
  rr
} ## end  plregr.control
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
                 "group", "smooth.group", "smooth.weights", "smooth.weight")
i.argPlcontr <-
  c("x", "y", "data", "transformed", "subset", ## "cex", "markextremes",
    "vcol", "vlty", "vpch", "plscale", ##  "smooth",
    "main", "sub", ".subdefault", ## "mar",
    "xlab", "ylab", "varlabels",
    "ploptions", ".environment.")

i.argPlControl <- ##!!!
  c("x", "y", "data", "xvar", "transformed", i.argPldata,
    "refline", "ploptions", 
    "main", "sub", "cex.main", "varlabels",
    ## plregr.control
    "xvar", "weights", "glm.restype", "smresid",
    "partial.resid", "cookdistlines", "leveragelimit", "condprobRange", 
    "testlevel", "refline", "smooth.sim",
    "xlabs", "reslabs", 
    "mf", "mfcol", "multnrow", "multncol", "multmar", "oma"
    )
## ====================================================================
i.argPlregr <- c("plotselect", "sequence", "addcomp", "smooth.legend")

plregr <-
  function(x, data=NULL, plotselect = NULL, xvar = TRUE,
           transformed = NULL, sequence=FALSE, weights=NULL, 
           addcomp = FALSE, smooth = 2, smooth.legend = FALSE, 
           markextremes = NA, mar = NULL, byrow = NULL, 
           plargs = NULL, ploptions = NULL,...)
{
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date:  7 May 93 / 2002
  argPlregr <- c("x", "data", "plotselect", ## "xvar", "transformed",
                 "sequence", "weights", "addcomp", "smooth.legend")
  lImer <- inherits(x, "merMod")
  if (lImer) x <- i.merprep(x)
  ## ----------------
  if (is.null(plargs)) {
    lcall <- match.call()
    lcall$x <- x
    lcall[1] <- list(quote(plgraphics::plregr.control))  ## need plgr:: because
    ## regr also uses this function
    mode(lcall) <- "call"
    plargs <- eval(lcall, parent.frame())
  }
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  ## -------------------------------------------------------------------
  ## all these results from  plregr.control  include the  na.action  observations
  lres <-  plargs$residuals
  lmres <- plargs$rescol
  lpldata <- plargs$pldata
  ##  ploptions <- plargs$ploptions
  lsmgrp <- lpldata$"(smooth.group)"
  lsmooth <- i.getploption("smooth")
  ## lsmpar <- i.getploption("smooth.par")
  lsigma <- plargs$sigma
  llev <- plargs$leverage
  lresname <- plargs$resname
  lsimres <- plargs$simres
  lnsims <- if (length(lsimres)==0) 0 else ncol(lsimres)
  llevlim <- plargs$leveragelimit
  lrefline <- i.def(plargs$ploptions$refline, TRUE)
  lsimstdres <- plargs$simstdres
  x$na.action <- lnaaction <- plargs$na.action
  ## from x
  lform <- plargs$formula ## formula(x)
  lweights <- naresid(lnaaction, x$weights)
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
    stop("!plregr! I do not know how to plot residuals of a multinomial regression")
  lglm <- inherits(x, "glm")
  lbinary <- lglm && length(unique(plargs$y))==2 ## binary binomial
  lcensored <- inherits(x[["y"]], "Surv")
  lnnls <- !inherits(x, "nls")
  if (u.true(plargs$famcount)) plargs$ploptions$smooth.iter <- 0
  condquant <- NULL
  lIcq <- u.true(i.getploption("condquant")) & (lbinary|lcensored)
  ## -----------------------------------
  ## plot selection
  lplsel <- unlist(
    i.plotselect(plotselect, smooth=plargs$smooth, Iwgt=lIwgt,
                 mult=lmres>1,
                 famgauss=plargs$famgauss, famglm=inherits(x,"glm"),
                 famcount=plargs$famcount)
    )
  lplsel <- lplsel[is.na(lplsel)|lplsel>0]
  lnplsel <- length(lplsel)
  ## --- fit
  lfit <- x$linear.predictor
  lfitname <- "linear predictor"
  if (inherits(x, "polr")) lfit <- fitted.regrpolr(x, type="link")
  else {
    if(is.null(lfit)) {
      lfit <- fitted(x)
      lfitname <- "fitted value"
    } else lfit <- naresid(lnaaction, lfit)
  }
  lfitname <- rep(lfitname, length=lmres)
  lir <- i.getploption("innerrange.fit") ## extra element of ploptions!
  lfit <- as.data.frame(lfit)
  lfit <- genvarattributes(as.data.frame(lfit), varlabels = lfitname,
                           innerrange.limits=lir)
  ## standardized residuals
  lstdres <- as.data.frame(attr(lres, "stdresiduals"))
  lstrratio <- as.data.frame(attr(lres, "stdresratio"))
  lresplab <- plargs$resplab
  ## cat("lresplab ", str(lresplab))
  lIrpl <- length(lresplab)>0
  if (any(lii <- c("qq","leverage","absresfit","absresweights")%in%names(lplsel))) {
    lIabs <- any(lii[3:4])
    labsres <- if (lIabs) lstdres
    lstdresname <- paste("st.", lresname, sep = "")
    labsresname <- paste("|st.",lresname,"|", sep="")
    ##
    for (lj in seq_len(lmres)) {
      lrsj <- lres[,lj]
      ## residuals from smooth
      if (plargs$smresid) { ## !!! transfer to plregr.control
        lfsm <- gensmooth(lfit[,lj], lrsj, plargs=plargs, ploptions=ploptions) 
        lfsmr <- residuals(lfsm)
        if (lna <- sum(is.na(lfsmr)&!is.na(lres[,lj])))
          warning(":plregr: residuals from smooth have ",
                  round(100*lna/lnobs,1), " % additional NAs")
        ## lstdres[,lj] <-
        lrsj <- lfsmr * lstrratio[,lj] ## !!! * f(leverage of smooth)
        if (lIcq) {
          lcq <- attr(lrsj, "condquant")
          if (length(lcq)) {
            li <- lcq[,"index"]
            lcq[,1:4] <- ( lcq[,1:4]-lcq[,1]+lfsmr[li] ) *lstrratio[li,lj]
            ## attr(lfsmr, "condquant") <- lcq
##            attr(lstdres[,lj], "condquant") <- lcq
            attr(lrsj, "condquant") <- lcq
          }
        } else {
          if (lIcq) {
            lcq <- attr(lrsj, "condquant")
            if (length(lcq)) {
              lcq[,1:4] <- lcq[,1:4]*lstrratio[lcq[,"index"],lj]
              attr(lrsj, "condquant") <- lcq
            }
          }
        lstdres[,lj] <- lrsj
        }
        if (lnsims) lsimstdres <- lsimres * lstrratio[,lj] ## index needed for multiv
      } ## fi plargs$smresid
      if (lIabs) labsres[,lj] <- abs(lrsj)
    }
    if (plargs$smresid) {
      lstdresname <- paste("st.sm.", lresname, sep = "")
      labsresname <- paste("|st.sm.",lresname,"|", sep="")
    }
    names(lstdres) <- lstdresname
    lstdres <- genvarattributes(lstdres, ploptions=ploptions)
    if (lIabs) {
      labsres <-
        genvarattributes(labsres, ploptions=ploptions, varlabels=labsresname)
      for (lj in seq_len(lmres))
        attr(labsres[,lj], "plrange")[1] <- 0
    }
  }
  ## --- multiple frames xxx
  lmf <- i.def(plargs$mf, NULL, TRUE, NULL)
  if (length(lmf)) {
    if (u.true(lmf)) ## is.logical(lmf)&&lmf)
      lmf <- if (lmres>1) {
               if (lmres<=4) c(lmres, NA) else lmres
             } else lnplsel
##               if (lnplsel<=2) lnplsel else c(2)
  }
  if (length(lmf)==2) {
    if (is.na(lmf[1])) {
      lmf1 <- ceiling(lnplsel/lmf[2])
      if(lmf1>lmf[2]+1) lmf1 <- lmf[2]
      lmf[1] <- lmf1
    }
    if (is.na(lmf[2])) {
      lmf1 <- ceiling(lnplsel/lmf[1])
      if(lmf1>lmf[1]+1) lmf1 <- lmf[1]
      lmf[2] <- lmf1
    }
  }
  lbyrow <- i.def(ploptions$byrow, FALSE)
  loma <- i.def(plargs$oma, c(3,3,2,1)*(length(lmf)>0), valuefalse=NULL)
  ##              c(0,0,2,1)*(length(lmf)>0)+c(3,3,0,1)*(lmres>1),
  if (length(loma)<4) loma <- c(0,0,loma,0)[1:4]
  loldpar <-
    c(par(c("cex","mar","mgp")),
      if (length(lmf)&(!is.logical(lmf))) {
        lop <-  
          if (length(lmf)==1)
            attr(plmframes(mft=lmf, oma=loma, byrow=lbyrow),"oldpar") else
        attr(plmframes(lmf[1], lmf[2], oma=loma, byrow=lbyrow),"oldpar")
        c(par("mfrow"), lop[setdiff(names(lop), c("mfig","mrow","mcol"))])
      } ## else par(oma=loma) ## , ask=plargs$ask
      )
  loldpar <- loldpar[!duplicated(names(loldpar))]
  on.exit(par(loldpar), add=TRUE)
  ## reduce  mar[3]
  ploptions$mar <- pmax(i.getploption("mar")-c(0,0,1.5,0), 0)
  lcex <- par("cex") ## *i.getploption("cex") is done in plframe
  ## par(cex=lcex)
  lmar <- if (lmres>1 && lmres==par("mfg")[3])
            plargs$multmar else  i.getploption("mar") ## multmar set by .control
  plargs$ploptions$mar <- lmar
  lnewplot <- TRUE ## !!!
  ## --------------------------------------------------------------------------
  ## start plots
  if (length(lplsel))
    for (liplot in 1:length(lplsel)) {
      lpllevel <- lplsel[liplot]
      lpls <- names(lpllevel)
      ## --- y on fit
      if(lpls=="yfit") {
##-         lsml <- if (length(lsmlegend["yfit"]))
##-                   setNames(rep(lsmlegend["yfit"],length(plargs$yname)),
##-                            plargs$yname)  else lsmlegend
        ly <- plargs$response ##attr(lpldata, "yvar")
        if (length(ly)==0)
          warning(":plregr: response not found")
        else {
          plargs$smooth <- lpllevel-1
##          plargs$mar <- lmar
          plargs$reflinecoord <- c(x=median(lfit[,lj], na.rm=TRUE),y=1)
          for (lj in seq_len(ncol(ly))) {
            lyj <- ly[,lj]
            lfj <- lfit[,lj]
            if(length(lsimres))
              lyj <- structure(data.frame(lyj, lfj+lsimres), primary=1)
            plpanel(x=lfj, y=lyj, frame=TRUE, ploptions=ploptions)
          }
        }
        plargs$reflinecoord <- NULL
      }
      ## --- Tukey Anscombe plot
      if(lpls=="resfit") {
        for (lj in seq_len(lmres)) {
          plargs$smooth <- lpllevel-1
          plargs$reflinecoord <- c(x=median(lfit[,lj], na.rm=TRUE),y=-1)
          lrsj <- lres[,lj]
          if(length(lsimres))
            lrsj <- structure(data.frame(lres[,lj], lsimres), primary=1)
          plpanel(lfit[,lj], lrsj, plargs=plargs, title=NA, frame=TRUE)
        }
        plargs$reflinecoord <- NULL
        par(cex=lcex)
      }
      ## --- scale plot
      if(lpls=="absresfit")
        if(length(labsres)==0) 
          warning(":plregr: No standardized residuals found")
        else {
          for (lj in seq_len(lmres)) {
            labsrj <- labsres[,lj, drop=FALSE]
            if (lnsims)
              labsrj <- structure(data.frame(labsrj, abs(lsimstdres)), primary=1)
            plpanel(lfit[,lj], labsrj, frame=TRUE, title=NA, 
                    plargs=c(plargs, list(smooth.power=0.5)) ) ## plsmooth needs 'power'
          }
          par(cex=lcex)
        }
      ## --- plot abs. res vs. weights
      if(lpls=="absresweights") {
        if (length(lweights)!=lnr)
          warning(":plregr: no suitable weights found.",
                  "cannot plot absres on weights")
        else { ## copy from absresfit
          lwg <- lweights
          lwg[lwg<=0] <- NA
          lwg <- genvarattributes(data.frame(lwg))
          if (attr(lwg[,1], "plrange")[1]<0.3) attr(lwg[,1], "plrange")[1] <- 0
          for (lj in seq_len(lmres)) {
            labsrj <- labsres[,lj, drop=FALSE]
            if (lnsims)
              labsrj <- structure(data.frame(labsrj, abs(lsimstdres)), primary=1)
            plpanel(lwg, labsrj, frame=TRUE, title=NA,
                    ylab=paste("|",lresname,"| * sqrt(w)", sep=""),
                    plargs=c(plargs, list(smooth.power=0.5))) ## plsmooth needs 'power'
          }
        par(cex=lcex)
        }
      }  
      ## --- normal plot qq plot
      if(lpls=="qq") {
        lnsims <- if (length(lsimstdres)) ncol(lsimstdres) else 0 ##plargs$smooth.sim
        if (lnsims)
          lsimstdr <-
            if (i.def(attr(lsimstdres, "type"), "resampled")=="resampled")
              simresiduals.default(x, nrep=lnsims, simfunction=rnorm,
                                   sigma=apply(lstdres, 2, mad, na.rm=TRUE) )
            else lsimstdres
        for (lj in seq_len(lmres)) {
          llr <- lstdres[,lj]
          lio <- order(llr)[seq_len(lnobs)]
          llr <- transferAttributes(llr[lio], llr)
          if (length(lat <- attr(llr, "numvalues")))
            attr(llr, "numvalues") <- lat[lio]
          if (length(lat <- attr(llr, "plcoord")))
            attr(llr, "plcoord") <- lat[lio]
          lIcqj <- length(lcq <- attr(llr, "condquant"))>0 ##!!! falsch ?
          lIcqu <- lIcq | lIcqj
          lpch <-
            if (lIcqu) i.pchcens(plargs, lcq)[order(llr)] else ploptions$pch[1]
          lxx <- qnorm(ppoints(lnobs))
          attr(lxx, "zeroline") <- 0
          plframe(lxx, llr, xlab = "theoretical quantiles",
                  ylab = lstdresname[lj], mar=lmar, plargs=plargs)
          ##-  lxy <- qqnorm(llr, ylab = lstdresname[lj], main="", type="n", )
          if (lnsims>0) {
            for (lr in 1:lnsims) {
              llty <- last(i.getploption("smooth.lty"))
              lines(lxx,sort(lsimstdr[,lr]), lty=llty,
                    lwd=i.getploption("linewidth")[llty],
                    col=last(i.getploption("smooth.col")) )
            }
          }
          plpoints(lxx, llr, plargs=list(pldata=plargs$pldata[lio,]),
                   ploptions=ploptions, setpar=FALSE)
          lquart <- quantile(llr, c(0.25,0.75), na.rm=TRUE)
          ## qq line
          plrefline(c(0, diff(lquart)/(2*qnorm(0.75))), plargs=plargs)
          if(lIcq & lj==lmres)
            legend("bottomright",
                   pch=c(rep(ploptions$censored.pch,length=2)),
                   legend=c("uncensored","censored"))
          pltitle(plargs=plargs, show=FALSE)
        }
        par(cex=lcex)
      }
      
      ## --- leverage plot. If weight are present, use "almost unweighted h"
      if(lpls=="leverage")
      if ((!is.na(lpllevel))&&lpllevel>0) {
        if (lIwgt) llev <- llev/lweights
        if (diff(range(llev,na.rm=TRUE))<0.001)
          notice("plregr: all leverage elements equal, no leverage plot")
        else {
          llevpl <- genvarattributes(
            data.frame(leverage = plcoord(llev, c(0,llevlim[2]), ploptions=ploptions)) )[,1]
          lstdres <- genvarattributes(lstdres, varlabels = lstdresname)
          ## mark extremes
          lmx <- i.getploption("markextremes")
          if (is.list(lmx)) lmx <- lmx[["(lev)"]]
          if (is.function(lmx)) lmx <- lmx(lnobs)
          lmx <- last(i.def(lmx, 0))
          if (lImxlev <- lmx>0)
            lpllev <- plmark(llev, markextremes=c(0,lmx), plabel=plargs$plabel)
          ##
          llevtit <- paste("leverage", if(lIwgt) "(de-weighted)")
          ldfmod <- plargs$df.model
          lcookl <- plargs$cookdistlines
          if (i.def(ldfmod, 0)<=1) {
            warning(":plregr: model degrees of freedom <=1. No Cook-distance lines")
            lcookl <- NULL
          }
          if (lIcook <- length(lcookl)>0) {
            llx <- seq(0, max(llev,na.rm=TRUE), length=50)
            ## see formula for curves of constant Cook's distance in terms of
            ##   standardized residuals
            llrcd <- outer(sqrt((1-llx)/((ldfmod-1)*llx)), c(lcookl,-lcookl)) 
          }
          for (lj in 1:lmres) {
            lstrj <- lstdres[,lj]
            lplj <- if (lIrpl)  lresplab[,lj]  else  rep("", lnr)
            if (lImxlev) {
              lplj <- ifelse(lpllev=="", lplj, lpllev)
              if (any(lplj!=""))  plargs$pldata$.plab. <- lplj
            }
            if (lIcook) plargs$reflinecoord <- list(x=llx, y=llrcd)
            lplopt <- ploptions(smooth=FALSE, assign=FALSE)
            plpanel(llevpl, lstrj, xlab=llevtit, plargs=plargs, ploptions = lplopt,
                    frame=TRUE)
##-             plframe(llevpl, lstrj, xlab=llevtit, mar=lmar, plargs=plargs)
##-             if (lIcook) plrefline(list(x=llx, y=llrcd), x=llevpl, plargs=plargs)
##-             lplj <- if (lIrpl)  lresplab[,lj]  else  rep("", lnr)
##-             if (lImxlev) lplj <- ifelse(lpllev=="", lplj, lpllev)
##-             plpoints(llevpl, lstrj, psize=if(lIwgt) lweights, ## condquant=0,
##-                      plab=if(lIrpl|lImxlev) lplj, plargs=plargs, setpar=FALSE)
##-             pltitle(plargs=plargs, show=FALSE)
          }
        }
        par(cex=lcex)
      }
      
      ## -----------------------------------------------------------------
      ## --- multivariate:
      ## residual matrix for multivariate regr
      if(lpls=="resmatrix") {
        lpanel <- function(x, y, indx, indy, pch, col, plab, panelargs=plargs, ...) {
          plpoints(x,y, plargs=panelargs, setpar=FALSE)
        }
        if (plargs$rescol>1) {
          lpa <- plargs
          lpa$pldata <- lpa$residuals
          lpa$mar <- NULL
          plmatrix(lres, panel=lpanel, plargs=lpa) #main=plargs$main, pch=plargs$pch
        }
      }
      ## --- mahalanobis residuals
      if(lpls=="qqmult")  ## qq plot of Mahalanobis lenghts for multivariate regr
        if ((!is.na(lpllevel))&&lpllevel>0) {
          lresmahal <- plargs$resmahal
          if (is.null(lresmahal)) lresmahal <- mahalanobis(lres,0,var(lres,na.rm=TRUE))
          lxx <- sqrt(qchisq(ppoints(lresmahal),ncol(lres)))
          lio <- order(llr)[seq_len(lnobs)]
          llr <- transferAttributes(llr[lio], llr)
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
  ## ----------------------------------------------------------
  plargs$mf <- FALSE ## avoid a new page
  ## plargs$ylim <- lylim  ## no need to calculate again
  plargs$ploptions$mar <- lmar + c(0,0,1.5,0)
  lxvar <- if (u.notfalse(xvar)) attr(lpldata, "xvar")
  if (!i.def(sequence, FALSE)) lxvar <- setdiff(lxvar, ".sequence.")
  ## avoid plot against x in simple regression
  if (length(lxx <- setdiff(lxvar, c(".sequence.", ".weights.")))==1 &&
      "resfit"%in%names(lplsel)) {
   ## ltr <- i.def(transformed, FALSE) || lxx%in%
    if (lxx %in% names(lpldata)) {
      message("plregr: plot of residuals on  ", lxx,
              "  not shown because it is equivalent to 'resfit'")
      lxvar <- setdiff(lxvar, lxx)
    }
  }
  if (length(i.def(lxvar, NULL, valuefalse=NULL)))
    plresx(x, data=data, resid=lres, xvar=lxvar, 
           transformed = transformed, sequence=sequence,
           weights= if ("weights"%in%names(lplsel)) FALSE else NULL,
           addcomp = addcomp, smooth.legend=lsmlegend,
           plargs = plargs)
  ## --- end plregr
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
        warning(":plregr: Inadequate argument plotselect")
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
        warning(":plregr: Inadequate elements in plotselect: ",
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
} ## end i.plotselect
## ==========================================================================
plresx <-
  function(x, data = NULL, xvar = TRUE, transformed = NULL,
           sequence = FALSE, weights = NULL, 
           addcomp = FALSE, smooth = 2, smooth.legend = FALSE,
           markextremes = NA,
           plargs = NULL, ploptions = NULL, ...)
## ------------------------------------------------------------
{ ## plresx
  lcall <- match.call()
  ## lIplargs <- is.null(lcall[["plargs"]])
  if (length(plargs)==0) {
  ## plargs <-
  ##  if (lIplargs) {
    lac <- as.list(lcall)[-1]
    lac$stdresid <- FALSE
    ladrop <- c("sequence", "weights", "addcomp", "smooth.legend")
    lcall <-
      c(list(quote(plregr.control)), lac[setdiff(names(lac), ladrop)])
    mode(lcall) <- "call"
    plargs <-eval(lcall, parent.frame())
  } ## else  eval(lcall[["plargs"]]) ## , sys.frame(1))
  ##  plargs$ploptions$smooth <- i.getplopt(smooth)
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  plargs$ploptions$mar <- ploptions$mar <-
    pmax(i.getploption("mar")-c(0,0,1.5,0), 0)
  ## -------------------------------------------------------------------
  if (inherits(x,"mulltinom"))
    stop("!plresx! I do not know how to plot residuals of a mulitnomial regression")
  if (is.null(transformed))
    transformed <- i.def(plargs$transformed, length(xvar)>0)
  ##  reference lines
  lshrefl <- i.getploption("refline") & !inherits(x,"coxph")
  ## data
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
    pldata$".sequence." <- structure(1:lnr, varlabel="sequence", zeroline=FALSE)
    lvars <- c(lvars,".sequence.")
  }
  ##  --- weights  as x variable
  lIweights <- i.def(weights, lIwgt, TRUE, FALSE)
  if (lIweights)
    if (!lIwgt) 
      warning(":plresx; No weights found.",
              " Cannot plot residuals against weights")   else {
    pldata[,".weights."] <- naresid(lnaaction, plargs$weights)
    lvars <- c(lvars, ".weights.")
  }
  ## ------------------
  lnvars <- length(lvars)
  if (lnvars==0) {
    warning(":plresx: I did not find any x variables")
    return() }
  ## terminouterel
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
  terminouterel <- lvars%in%lvmod
  lvcomp <- intersect(lvars, lvmod)
  ## type
  addcomp <- as.logical(i.def(addcomp, FALSE, TRUE, FALSE))
  ## lpa <- plargs
  lInnls <- !inherits(x, "nls")
## -----------------------------------
  ## fit components for refline
  lIcomp <- addcomp|lshrefl
  if (length(lvcomp)) {
    if (any(terminouterel) && lInnls) {
      if (length(x$call$data)==0) x$call$data <- pldata
      lcmp <- fitcomp(x, vars=lvcomp, transformed=transformed,
                      xfromdata=FALSE, se=lshrefl>1)
      lcompx <- lcmp$x
      lcompy <- if (addcomp) lcmp$comp else -lcmp$comp
      lcompse <- lcmp$se
      if (addcomp) {
        lcompdt <-
          fitcomp(x, pldata, vars=lvcomp, transformed=transformed,
                  xfromdata=TRUE)$comp
        ## !!! add to lres, careful for condq
      }
    } else lshrefl <- FALSE
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
    lpa <- eval(i.getploption("panel"))
    if (is.character(lpa)) lpa <- get(lpa, envir=globalenv())
    lpanel <- function(x, y, indx, indy, pch, col, panelargs=NULL, ...) {
      lcmpx <- lcmpy <- NULL
      ltin <- terminouterel[indx]
      lvx <- lvars[indx]
      lcnt <- !is.factor(lvx)
      if (ltin&lshrefl) {
        lcmpy <- lcompy[,lvx,indy]
        if (lcnt) lcmpx <- lcompx[,lvx]
      }
      lpa(x, y, plargs=plargs)
    }
    plmatrix(pldata[,lvars,drop=FALSE], lres, panel=lpanel,
             ##pch=plargs$plab, plcol=plargs$pldata$plcol,
             nrow = plargs$multnrow, ncol = plargs$multncol,
             mar = plargs$multmar,
             plargs=plargs) 
    return()
  }
  ## ------------------------------------------------------------------
  lmf <- i.def(plargs$mf, NULL, TRUE, NULL)
  if (length(lmf)) {
    if (u.true(lmf))
      lmf <- if (lnvars<=6) lnvars else
        min(lnvars,ceiling(lnvars/((lnvars-1)%/%6+1)))
  }
  lbyrow <- i.def(ploptions$byrow, FALSE)
  loma <- i.def(plargs$oma, c(2,1)*(length(lmf)>0), valuefalse=NULL)
  if (length(loma)<4) loma <- c(0,0,loma,0)[1:4]
  loldpar <-
    c(par(c("cex","mar","mgp")),
      if (length(lmf)&(!is.logical(lmf))) {
        lop <- 
          if (length(lmf)==1)
            attr(plmframes(mft=lmf, oma=loma, byrow=lbyrow),"oldpar") else
        attr(plmframes(lmf[1], lmf[2], oma=loma, byrow=lbyrow),"oldpar")
        c(par("mfrow"), lop[setdiff(names(lop), c("mfig","mrow","mcol"))])
      } ## else par(oma=loma) ## , ask=plargs$ask
      )
  on.exit(par(loldpar), add=TRUE)
  lcex <- par("cex") ## *i.getploption("cex")
  ## par(cex=lcex)
  ##
  lmbox <- ploptions$factor.show=="mbox"
  lIjitter <- !lmbox ## ploptions$factor.show=="jitter"
  ljitfac <- i.getploption("jitter.factor")
  lrpl <- plargs$resplab
  if (lIrpl <- length(lrpl)>0) lrpl <- lrpl[,1]
  ## lpla <- plargs
  lr <- lres
  if (lIsmooth && lnsims>0) lr <- cbind(lr, lsimres)
  lrs <- lr ## need a copy for the case  addcomp  is true
##-         if (inherits(lr, "condquant"))
  ##-           lr <- lr[,1]
  ## ------------------------------------------------------------------
  ## --- loop --- plresx
  for (lj in 1:lnvars) {
    lvr <- lvars[lj] ## if (transformed) lvars[lj] else lrawv[lj]
    lv <- unname(lvr)
    lcmpj <- terminouterel[lj] && lshrefl && lInnls
    lci <- if (lcmpj&lshrefl) lcompy[, lvr] else 0 ## !!!
    if (plargs$partial.resid) 
      if (lshrefl && addcomp && lcmpj)  lrs <- lr+lcompdt[, lvr]
    lvv <- pldata[, lv]
    ##    lpa$pldata <- cbind(lr, lvv, pldata)
    mar <- i.def(plargs[["mar"]], c(NA, par("mar")[-1]))
    ## ---
    if (is.null(attr(lvv, "varlabel")))  attr(lvv, "varlabel") <- lv
    lrs1 <- lrs[,1,drop=FALSE]
    if (is.factor(lvv)) {
      ## factors
      ## lmar <- par("mar") ## might be adapted to nchar(levels)
##      lvv <- structure(lvv, varlabel=attr(lvv, "varlabel"))
      ll <- levels(lvv)
      lnl <- length(ll)
      if (lmbox) {
        ##       if (lIcq) lrs <- lres[,"random"]
##        attr(lpla$pldata,"xvar") <- lv
        plmboxes.default(lvv, lrs1, data=pldata, mar=mar,
                 plargs=structure(plargs, pldata=structure(pldata, xvar=lv)))
      } else {
##        attr(lvv, "plcoord") <- jitter(as.numeric(lvv), factor=ljitfac)
        plframe(lvv,lrs1, ploptions=ploptions)
      }
      ## reference values !!! simplify, defer to plpanel ?
      if (lcmpj) {
        lx <- seq_along(ll)
        lcil <- lci[1:lnl]
        lrsrg <- attr(lrs1,"innerrange")
        if ((!is.null(lrsrg))&&diff(lrsrg)) {
          lcilp <- plcoord(lcil, lrsrg, innerrange.ext=plargs$innerrange.ext)
          if (any(attr(lcilp,"nouter"))) {
            liout <- lcilp!=lcil
            lrsout <- pmax(pmin(lcilp[liout], lrsrg[2]), lrsrg[1])
            segments(lx[liout]-0.45, lrsout, lx[liout]+0.45, lrsout,
                     lty=ploptions$refline.lty, lwd = ploptions$refline.lwd,
                     col=ploptions$refline.col[1])
            lcil[liout] <- NA
          }
        }
        segments(lx-0.45, lcil, lx+0.45, lcil,
                 lty=ploptions$refline.lty[1], lwd = ploptions$refline.lwd[1],
                 col=ploptions$refline.col[1])
        if (lshrefl>1) {
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
      if (!lmbox) plpoints(lvv,lrs1, plargs=plargs, setpar=FALSE)
    } else { # ---
      ## --- continuous explanatory variable
      ##lrs <- lres[,lj]
      if (lnsims>0) lrs1 <- structure(cbind(lrs1, lsimres), primary=1)
      if (lshrefl && lcmpj) {
        lrefx <- lcompx[,lvr]
        lrefyb <- if (lshrefl>1) outer(lqnt*lcompse[,lvr], c(-1,1))
        plargs$reflinecoord <- list(x=lrefx, y=lci, band=lrefyb)
      }
      plpanel(lvv, lrs1, plargs=plargs, title=NA, frame=TRUE)
    } ## ends  if factor else
    par(cex=lcex)
  }
  invisible(plargs)
} ## end plresx
## ==========================================================================
smoothRegr <-
  function(x, y, weights=NULL, par=NULL, iterations=50, minobs=NULL, ...)
{
  minobs <- i.def(minobs, i.getploption("smooth.minobs"), valuefalse=1)
  lnx <- sum(is.finite(x))
  if (lnx<minobs ) return(NULL)
  liter <- max(i.def(iterations, i.getploption("smooth.iter")),0)
  lform <- if (is.formula(x)) x else y~x
  lpar <- i.def(par,i.getploption("smooth.par"))
  if (is.function(lpar)) lpar <- lpar(lnx)
  ## ----------------------------------------------------------------
  lcall <- call("loess", formula=lform, data=data.frame(x=x, y=y),
              weights=weights, span=lpar, iterations=max(liter,1),
              family=if (liter>0) "symmetric" else "gaussian",
              na.action=na.exclude)
  if (is.null(weights)) lcall$weights <- NULL
  lsm <- if (ploptions("debug")) eval(lcall, parent.frame())
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
           plargs = NULL, ploptions = NULL, ...)
{
  ## Purpose:   smooth for multiple y : one column from data, the other sim
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  9 Feb 2016, 14:57
  if (length(plargs)==0) plargs <- get(".plargs", globalenv())
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  ##
  lsmfunc <- i.getploption("smooth.function")
  if (is.character(lsmfunc)) lsmfunc <- get(lsmfunc)
  if (is.null(lsmfunc)) lsmfunc <- smoothRegr
  power <- i.def(power, 1,1,1)
  ## ---
  lnx <- NROW(x)
  if (inherits(y, "Surv")) y <- y[,1] ## !!! needs improvement!
  ly <- as.matrix(y)
  if (nrow(ly)!=lnx) stop("!gensmooth! Incompatible dimensions of 'x' and 'y'")
  ## if (length(weights)<=1) weights <- rep(1, lnx)
  lweights <- plargs$pldata$"(smoothWeights)"
  lIwgt <- length(lweights)>0
  if (lIwgt&&length(lweights)!=lnx)
    stop("!gensmooth! Incompatible dimensions of 'x' and 'weights'")
  lgroup <- plargs$pldata$"(smooth.group)"
  if (is.null(lgroup)) lgroup <- plargs$pldata$".group."
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
  if (is.function(lpar)) lpar <- lpar(lnobs)
  lparband <- lpar[1]* i.def(lpar[2], 1.5, 1.5, 1)
  liter <- i.getploption("smooth.iter")
  lminobs <- i.getploption("smooth.minobs")
  ## data: look for numvalues
  lx <- i.def(attr(x,"numvalues"),x)
##  if (inherits(lx, "POSIXt")) lx <- as.numeric(lx) 
  ly <- apply(ly,2, function(y) i.def(attr(y, "numvalues"), y))
  lnna <- apply(cbind(lx,ly), 1, sumNA)
  lx[lnna>0] <- NA
  lio <- order(as.numeric(lgrp), lx) ## order by group
  lio <- lio[!is.na(x[lio])]
  lxo <- lx[lio] # sorted without NA
  lyo <- ly[lio,,drop=F]
  lgrpo <- lgrp[lio]
  lngrp <- length(levels(lgrpo))
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
    lxg <- lxo[lig]
    if (sum(!is.na(lxg))<lminobs) {
      notice("gensmooth: too few non-missing observations",
             if(lngrp>1) paste(" for group ",lgr) )
      next
    }
    for (j in ncol(lyo):1) {
      lsm <- lsmfunc(lxg, lyo[lig,j]^power,
                     weights=if(lIwgt) lwgto[lig] else NULL,
                     par=lpar[1], iterations=ploptions$smooth.iter, ...)
      if (length(lsm)==0) {
        notice("gensmooth: too few observations for a smooth")
      } else  lysm[lig,j] <- lsm^(1/power)
    }
    ## band
    if (band & length(lsm)) {
      lysmb <- rep(NA, length(lsm))
      lsmr <- lyo[lig,1]^power-lsm^power ## residual
      lsmrpos[lig] <- lip <- lsmr>=0
      ## high end
      lii <- which(lip)
      if(length(lii)) {
        lsmrh <- lsmr[lii]
      ##  ligi <- lig[lii]
        lsmh <- lsmfunc(lxg[lii], sqrt(lsmrh),
                        weights=if (lIwgt) lwgto[lig[lii]] else NULL,
                        par=lparband, iterations=liter)
        if (length(lsmh)==0) {
          notice("gensmooth: too few observations for a 'high' smooth",
                 if(lngrp>1) paste(" for group ",lgr) )
        } else lysmb[lii] <- lsmh^2
      }
      ## low end
      if(length(lii)) {
        lii <- which(!lip)
        lsmrl <- - lsmr[lii]
      ##  ligi <- lig[lii]
        lsml <- lsmfunc(lxg[lii], sqrt(lsmrl),
                        weights=if (lIwgt) lwgto[lig[lii]] else NULL,
                        par=lparband, iterations=liter)
        if (length(lsml)==0) {
          notice("gensmooth: too few observations for a 'low' smooth",
          if(lngrp>1) paste(" for group ",lgr) )
        } else lysmb[lii] <- - lsml^2
      }
      ## resulting band
      lysmband[lig] <- (lysmb + lsm^power)^(1/power)
    }
  }
  lysmin <- matrix(NA, lnx, ncol(lyo), dimnames=list(names(x),colnames(lyo)))
  lysmin[lio,] <- lysm
  lres <- if (resid==2) ly/lysmin else ly-lysmin
  rr <- list(x = lxo, y = lysm, group = if(!lInogrp) factor(lgrpo),
             index = lio, xorig = x, ysmorig = lysmin, residuals = lres,
             xtrim = attr(lsm, "xtrim") )
  if (band) rr <- c(rr, yband = list(lysmband), ybandindex = list(lsmrpos) )
  rr
} ## end gensmooth
## ==========================================================================
smoothLm <- function(x, y, weights = NULL, ...) {
  rr <- if (is.null(weights)) lm.fit(cbind(1,x), y, ...)$fitted
        else  lm.wfit(cbind(1,x), y, weights, ...)$fitted
  structure(rr, xtrim=0)
}
## smoothRegrrob <- function(x,y,weights,par=3*length(x)^log10(1/2),iter=50)
## =======================================================================
plmatrix <-
  function(x, y=NULL, data=NULL, panel=NULL, panelargs = plargs, 
           nrow=0, ncol=nrow, reduce=TRUE, keeppar=FALSE,
           xaxmar=NULL, yaxmar=NULL, xlabmar=NULL, ylabmar=NULL,
           xlab=NULL, ylab=NULL, ## partial match!?!
           oma=NULL, ## mar=NULL,
           cex=NULL, diaglabel.cex = NULL, ## cex needed because otherwise
           ## cex will be translated into diaglabel.cex
           plargs = NULL, ploptions = NULL, ...) 
{
## Purpose:    pairs  with different plotting characters, marks and/or colors
##             showing submatrices of the full scatterplot matrix
##             possibly on seeral pages
## -------------------------------------------------------------------------
  lcall <- match.call()
##  lIplargs <- is.null(lcall[["plargs"]])
  if (is.null(plargs)) {
    lcall[[1]] <- quote(pl.control)
    lcall$x <- x  ## needs evaluation
    lcall$.subdefault <- format(as.expression(substitute(x))) ##as.character(as.expression(substitute(x)))
    lcall$y <- y
    lcall$data <-
      if(length(data)) {
        if (is.name(substitute(data)))
          lcall$.subdefault <- format(as.expression(substitute(data)))
        data
      }
    lcall$gensequence <- FALSE
    lcall$.environment. <- parent.frame()
    ##
    plargs <- eval(lcall, envir=parent.frame())
##    plargs$main <- lmain
  }
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  ## --------------
  lf.axis <- function(k, x, axm, labm, txt, ...) {
    if (k %in% axm) plaxis(k, x, varlabel="", ploptions=ploptions, setpar=FALSE)
    if (k %in% labm) 
      mtext(txt, side=k, line=(1.5+(k %in% axm)), cex=lcex, ...)
  }
  lf.eq <- function(v1,v2) {
    if (is.factor(v1)) is.factor(v2)&& all(dropNA(as.numeric(v1)==as.numeric(v2)))
    else all(dropNA(v1==v2))
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
  ## margins, will be used by  plframe
  panel <- i.getplopt(panel)
  if (is.character(panel)) panel <- get(panel, globalenv())
  if (!is.function(panel)) {
    warning(":plmatrix: 'panel' not found. Using 'plgraphics::plpanel'")
    panel <- plgraphics::plpanel
  }
  if (length(xlab)|length(ylab))
    warning(":plmatrix: Arguments 'xlab' and 'ylab' not used. ",
            "Set 'varlabels' instead!")
  ## -----------------------------------------------
  ## data and bookkeeping
  pldata <- plargs$pldata
  lgrp <- c(".group.",".pcol.") %in% names(pldata)
  ## if  group , generate  pcol
  if (".group."%in%names(pldata)) {    
    if (!".pcol."%in%names(pldata))
      plargs$pldata[,".pcol."] <- pldata[,".group."]
##-     if (!"(smooth.group)"%in%names(pldata))
##-       plargs$pldata[,"(smooth.group)"] <- pldata[,".group."]
  }
  if (is.null(i.getploption("main"))) plargs$main <- plargs$dataLabel
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
  lnobs <- lnr-mean(sumNA(ldata))
  lvsurv <- sapply(ldata, function(x) inherits(x, "Surv") )
  if (any(lvsurv)) { ## survival vars
    lf.surv <- function(dt) structure(dt[,1], pch=dt[,2]+1)
    ldata[lvsurv] <- lapply(ldata[lvsurv], lf.surv)
  }
  ## title !!!
  lsub <- plargs$sub
  lmain <- plargs$main
  lcexmain <- i.getploption("title.cex")
  ## title: how many lines? 
  lltit <- (length(lmain)>0&&lmain!="") * sum(lcexmain[1:2])
  if (lltit==0) {
    if(length(lsub)>0&&as.character(lsub)!="") lltit <- lcexmain[2] ## + 0.2
  } else lltit <-  lltit ## +0.2
  ## --- position of tick marks and axis labels, oma
  xaxmar <- i.def(xaxmar, 1+(nv1*nv2>1))
  xaxmar <- ifelse(xaxmar>1,3,1)
  yaxmar <- i.def(yaxmar, 2+(nv1*nv2>1))
  yaxmar <- ifelse(yaxmar>2,4,2)
  xlabmar <- i.def(xlabmar, if (nv1*nv2==1) xaxmar else 4-xaxmar )
  ylabmar <- i.def(ylabmar, if (nv1*nv2==1) yaxmar else 6-yaxmar )
  if (length(oma)!=4)
    oma <- c(1+(xaxmar==1)+(xlabmar==1), 1+(yaxmar==2)+(ylabmar==2),
             1+(xaxmar==3)+(xlabmar==3) + lltit,  1+(yaxmar==4)+(ylabmar==4))
  ## set par
  ## if (!keeppar)
  loldp <- plmframes(nrow, ncol, nrow=nv2, ncol=nv1, oma=oma) # , mgp=c(1,0.5,0)
  on.exit(par(loldp), add=TRUE, after=FALSE)
  ## else par(mfrow=par("mfrow")) ## new page! if mfcol was set, it will not work
  ## on.exit(par(attr(loldp, "oldp"))) ## !!!
  lcex <- i.getploption("cex")*par("cex")
  lcexdiag <- i.getplopt(diaglabel.cex)*par("cex")
  ## ploptions$mar <-
  lmar <- rep(i.getploption("panelsep"), length=4) 
  lmfig <- par("mfg")
  lnr <- lmfig[3]
  lnc <- lmfig[4]
  lnpgr <- ceiling(nv2/lnr)
  lnpgc <- ceiling(nv1/lnc)
  ## cex.pch
  lcex.pch <- ploptions$cex.pch
  if (is.function(lcex.pch)) lcex.pch <- lcex.pch(lnobs)
  plargs$ploptions$cex.pch <- lcex.pch
  plargs$ploptions$axes <- FALSE
  lipanelargs <-
    intersect(names(as.list(args(panel))), c("indx","indy","panelargs"))
##
##-   ## log
##-   if (length(grep("x",log))>0) ldata[ldata[,1:nv1]<=0,1:nv1] <- NA
##-   if (length(grep("y",log))>0) ldata[ldata[,lv2+1:nv2]<=0,lv2+1:nv2] <- NA
  loldp <- c(par(cex=lcex), par(mar=lmar, mgp=i.getploption("mgp")))
  ## on.exit(loldp) ## already requested
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
          plframe(v1, v2, xlab="", ylab="", ticklabels = FALSE, mar=lmar,
                  ploptions=ploptions, setpar=FALSE) # plargs=plargs
          do.call(panel,
                  c(list(x=v1, y=v2), ## panel must have arguments x and y
                    list(indx=jd1, indy=jd2, panelargs=panelargs)[lipanelargs]) )
##          panel(v1,v2, indx=jd1, indy=jd2, plargs=plargs)
        }
        else {
          lv0 <- as.numeric(v1)
          plot(lv0,lv0, type="n", axes=FALSE, xlab="",ylab="")
          uu <- par("usr") # diagonal: print variable name
          text(mean(uu[1:2]),mean(uu[3:4]), lylab, cex=lcexdiag)
        }
##        usr <- par("usr")
        ##       lat=c(mean(usr[1:2]),mean(usr[3:4]))
        if (jr==lnr||jd2==nv2) lf.axis(1, v1, xaxmar, xlabmar, lxlab)
        if (jc==1) lf.axis(2, v2, yaxmar, ylabmar, lylab)
        if (jr==1) lf.axis(3, v1, xaxmar, xlabmar, lxlab)
        if (jc==lnc||jd1==nv1) lf.axis(4, v2, yaxmar, ylabmar, lylab)
        if (lltit>0)   pltitle(plargs=plargs, show=NA)
      } else frame()
    }
  }}
    ## stamp(sure=FALSE, outer.margin=TRUE) 
  }
  invisible(plargs)
} ## end plmatrix

## ====================================================================
plpanel <-
  function(x = NULL, y = NULL, indx=NULL, indy=NULL, type="p", frame = FALSE,
           title = FALSE, plargs = NULL, ploptions = NULL,
           ...) ## data=plargs$pldata
{
  if (length(plargs)==0) plargs <- get(".plargs", globalenv())
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  pldata <- plargs$pldata
  plargs$ploptions$stamp <- FALSE
  lshrefl <- i.getploption("refline")
  ##
  mbox <- i.getploption("factor.show")=="mbox"
  lIsm <- i.getploption("smooth")
  ## intro, needed if formulas are used or data is given or ...
##-   lcl <- match.call()
##-   lcall <- sys.call() ## match.call modifies argument names
##-   lcnm <- i.def(names(lcall), names(lcl))
  ##-   names(lcall) <- ifelse(lcnm=="", names(lcl), lcnm)
  lcall <- match.call()
  lxnm <- "x"
  lynm <- "y"
  if (is.formula(x)|is.formula(y)|
      any(c("data","pcol","psize")%in%names(lcall))) {
    lcall$assign <- FALSE
    lcall$ploptions <- ploptions
    if (is.null(lcall$data)) lcall$data <- pldata
    lplargs <- do.call(pl.control, as.list(lcall[-1]), envir=parent.frame())
    ploptions <- lplargs$ploptions
    plargs$pldata <- pldata <- lplargs$pldata
    x <- pldata[,lxnm <- attr(pldata, "xvar")[1], drop=FALSE]
    y <- pldata[,lynm <- attr(pldata, "yvar"), drop=FALSE]
  } else {
  ## ---
    if (length(x)==0) {
      x <- pldata[,2, drop=FALSE]
      lxnm <- names(pldata)[2]
    }
    if (length(y)==0) {
      y <- pldata[,1, drop=FALSE]
      lxnm <- names(pldata)[1]
    }
  }
  if (is.data.frame(x)) x <- x[,1]
  if (sum(!is.na(x))==0) {
    warning(":plpanel: no finite values of variable ", lxnm,
            ". Nothing to plot")
    return()
  }
  if (sum(!is.na(as.matrix(y)))==0) {
    warning(":plpanel: no finite values of variable(s) ", paste(lynm, collapse=", "),
            ". Nothing to plot")
    return()
  }
  if (is.character(x)) x <- factor(x) ## !!! attributes!
  if (is.character(y)) {
    if (NCOL(y)>1) {
      warning("!plpanel! multiple  y  cannot be of type character or factor",
              "Using the first column")
      y <- y[,1]
    }
    y <- factor(y)
  }
  if (length(dropNA(unique(x)))==1) x <- factor(x)
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
  lIyp <- attr(y, "primary")
  lIys <- attr(y, "secondary")
  if (length(lIyp) | length(lIys)) {
    lIyp <- i.def(lIyp, 1:NCOL(y))
    lIys <- i.def(lIys, setdiff(1:NCOL(y), lIyp))
  }
  lyp <- if (length(lIyp)) y[,lIyp, drop=FALSE] else y
  lys <- if (length(lIys)) y[,lIys] else NULL
  if (frame) plframe(x,lyp, ploptions=ploptions) ## !!! , setpar=FALSE
  ## secondary smooths
  if (lIsm & length(lys))
    plsmooth(x, ysec=lys, band=FALSE, power=plargs$smooth.power, plargs=plargs)
  ## refline
  if (lshrefl && length(lrfl <- plargs$reflinecoord))
    plrefline(lrfl, x=x, y=y, plargs=plargs)
  ## points
  plpoints(x, lyp, type=type, plargs=plargs, ploptions=ploptions, setpar=FALSE, ...)
  ## primary smooth
  if (lIsm) plsmooth(x, y=lyp, plargs=plargs)
  ## title
  if (u.notfalse(title)) pltitle(plargs=plargs, show=title)
} ## end plpanel
## ====================================================================
panelSmooth <-
  function(x, y, indx, indy, plargs = NULL, ploptions = NULL, ...) {
    if (length(plargs)==0) plargs <- get(".plargs", globalenv())
    graphics::panel.smooth(x, y, pch=plargs$pch, col=plargs$pcol,
                           cex=i.getploption("cex.pch"), ...)
  }
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
} ## end plmbox
## ====================================================================
plmboxes <- function(x, ...)
  UseMethod("plmboxes")

plmboxes.formula <-
  function(x, data, ...)
{
  ldt <- genvarattributes(getvariables(x, data))
##-   l1backback <- length(x[[3]])>1 && as.character(x[[3]][[2]])=="1"
##-   if (l1backback)
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
  function(x, y, data, width = 1, at = NULL, horizontal = FALSE,
    probs = NULL, outliers = TRUE, na = FALSE, backback = NULL,
    refline = NULL, add = FALSE, 
    xlim = NULL, ylim = NULL, axes = TRUE, xlab = NULL, ylab = NULL, 
    labelsperp = FALSE, xmar = NULL, 
    widthfac = NULL, minheight = NULL, colors = NULL, lwd = NULL,
    .subdefault = NULL, plargs = NULL, ploptions = NULL, 
    ...)
{
  ## Purpose:    multibox plot
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 14 Dec 2013, 23:38
  f.ylim <- function(ylm, ext)
    c((1+ext)*ylm[1]-ext*ylm[2], (1+ext)*ylm[2]-ext*ylm[1])
  ##
  lcall <- match.call()
  if (is.null(plargs)) {
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
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  ##
  title <- i.def(ploptions$title, TRUE)
##-   lzl <- i.getploption("zeroline")
##-   lzl <- list(x=FALSE, y=if(is.list(lzl)) lzl[[2]] else lzl)
  ## ----------------------------------------------------------------------
  pldata <- plargs$pldata
  if (is.null(x))
    x <- pldata[,i.def(attr(pldata,"xvar"),1), drop=FALSE] ## may have two columns
  else x <- as.data.frame(x)
  if (is.null(y))
    y <- pldata[,i.def(attr(pldata,"yvar")[1],2), drop=FALSE]
  ly <- if (is.data.frame(y)) y[,1] else y
##  ly <- i.def(i.def(attr(y, "numvalues"), attr(y, "plcoord")), y)
  if (!is.numeric(ly))
    stop("!plmboxes.default! 'y' must be numeric")
  lhoriz <- as.logical(i.def(horizontal, FALSE, valuetrue=TRUE))
##  if (lhoriz) lzl <- lzl[2:1]
##  plargs$ploptions$zeroline <- lzl
  ## widths
  lwfac <- modarg(widthfac, c(max=2, med=1.3, medmin=0.3, outl=NA, sep=0.003))
  ## colors, line widths
  lcol <- modarg(colors,
                 c(box="lightblue",med="blue",na="gray90",refline="magenta") )
  llwd <- modarg(lwd, c(med=3, range=2))
  ## data
  ## preliminary 
  llr <- ncol(x)>=2 ## asymmetric mboxes required for binary (second) factor
  if (llr && length(unique(x[,2]))!=2) {
    warning(":plmboxes: second x-variable must be binary. I ignore it.")
    x <- x[,1, drop=FALSE]
    llr <- FALSE
  }
  x[,1] <- lx <- transferAttributes(factor(x[,1]),x[,1], except="levels")
                                        # unused levels are dropped
  llist <- split(ly,x)
  llev <- levels(lx)
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
  if (lnat <- length(at)) {
    if (lnat!=lng) {
      warning(":plmboxes: 'x' has wrong length")
      at <- if (lnat>lng) at[1:lng] else NULL
    }
  }
  if (is.null(at)) at <- 1:lng
  backback <-
    i.getplopt(backback, ploptions) &&
      length(dropNA(unique(x[,1])))==2 && NCOL(x)==1
  if (backback) {
    llr <- TRUE
    lng <- 1
  }
##-     llev <- ""
##-     llev2 <- c(levels(x[,1]),"","")[1:2]
  ##-   } else
  lpos <- if(backback) rep(mean(at),2) else at
  ## probabilities
  if (is.null(probs))
      probs <- if (sum(!is.na(y))/(lng*(1+llr))<20) c(0.1,0.5,1)/2 else
               c(0.05,0.1,0.25,0.50,0.75,1)/2
  ## box for NA's?
  ##if (is.null(na)||is.na(na)||(is.logical(na)&&!na)) na.pos <- NULL else 
  ##    if (is.logical(na))
  na <- i.def(na, NA)
  na.pos <- i.def(na, c(min(ly, na.rm=TRUE)*(1-0.3)-0.3*max(ly, na.rm=TRUE)),
                  valuefalse=NULL)
  if (length(na.pos)==1)
    na.pos <- na.pos+ 0.03*diff(range(ly, na.rm=TRUE))*c(-1,1)
  lusr <- par("usr")
  ## plot range
  ##  lrg <- if (add) lusr[3:4] else attr(ly, "plrange")
  lirg <- attr(y, "innerrange")
  ljlim <- any(c(attr(y, "nouter"),0)>0)
  lyat <- i.def(attr(y,"ticksat"),
                pretty(ly, n=i.def(ploptions$tickintervals, 7)) )
  if (add) {
    xlim <- lusr[1:2]
    ylim <- lusr[3:4]
  } else {
    lxldf <- i.extendrange(
      range(at, na.rm=TRUE) + ## max(diff(at)) * 
        (max(width[c(1,length(width))]) - backback)*c(-1,1)*0.4) ## !!! is 0.4 ok???
      ## *i.getploption("plext")
    xlim <- i.def(xlim, lxldf, valuefalse=lxldf)
    if (any(li <- is.na(xlim))) xlim[li] <- lxldf[li]
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
    llev <- levels(x[,1])
    attr(lx, "ticklabelsat") <- attr(lx, "ticksat") <- seq_along(llev)
    attr(lx, "ticklabels") <- llev
    ## ---------------------------------
    if (lhoriz)
      plframe(y, xlim, plargs=plargs)
      else plframe(xlim, y, plargs=plargs)
    if (axes) {
      plaxis(1+lhoriz, lx, las=2*as.logical(labelsperp)) ## at=at, labels=llev, 
      if(lIna && anyNA(y)) 
        mtext("NA", 2-lhoriz,1, at=mean(na.pos), las=2) 
##-       if (backback) {
##-         mtext(llev2[1], 1+lhoriz,1, at=0.75)
##-         mtext(llev2[2], 1+lhoriz,1, at=1.25)
##-       }
      ##      mtext(xlab, side=1+lhoriz, line=ploptions$mar[1+lhoriz]-1)
      plaxis(2-lhoriz, x=ly)
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
  lattr <- attributes(ly)
  for (li in 1:lng) {
    if (is.na(at[li])) next
    lli <- llist[[li]]
    if (length(lli)) {
      attributes(lli) <- lattr
      plmbox(lli,lpos[li]-lsep, probs=probs, outliers=outliers,
             horizontal=lhoriz,
             wfac=lfac[li], adj=1-0.5*(1-llr), na.pos=na.pos, extquant=TRUE,
             widthfac=lwfac, colors=lcol, lwd=llwd, warn=-1)
    }
    if (llr) { ## second half of asymmetrix  mbox
      llir <- llist[[li+lng]]
      if (length(llir)) {
        attributes(llir) <- lattr
        plmbox(llir,lpos[li]+lsep,probs=probs, outliers=outliers,
               horizontal=lhoriz,
               wfac=lfac[li], adj=0, na.pos=na.pos, extquant=TRUE,
               widthfac=lwfac, colors=lcol, warn=-1)
      }
    }
  }
  pltitle(plargs=plargs, show=title)
  stamp(sure=FALSE)
  invisible(at) ##!!? return  plargs 
} ## end plmboxes
## ====================================================================
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
  function(formula=NULL, reg=NULL, data=NULL, restrict=NULL, size = 1,
           xlab = NULL, ylab= NULL, plextext = NULL, pale = 0.2,
           plargs = NULL, ploptions = NULL, ...)
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
  if (is.null(plargs)) {
    lac <- as.list(lcall)[-1]
    lac$x <- lac$reg
    lac$xvar <- lac$formula
    ladrop <- c("formula", "restrict", "size")
    lcall <- c(list(quote(plregr.control)),
               lac[setdiff(names(lac), ladrop)])
    mode(lcall) <- "call"
    plargs <- eval(lcall, parent.frame())
  } 
  if (length(ploptions)==0) ploptions <- plargs$ploptions
  ## ------------------------------------------------------------------
  lz <- plargs$residuals
  pldata <- plargs$pldata
  lxvar <- attr(pldata, "xvar")
  if (length(lz)==0) lz <- attr(pldata, "yvar")
  ##--- restrict z values: --- !!! use innerrange of z
  lrestr <-
    if (any(attr(lz, "nouter")>0)) max(abs(attr(lz, "innerrange"))) else NULL
  restrict <- i.def(i.def(restrict, lrestr), NULL, valuefalse=NULL)
  if(length(restrict)==0)   restr <- FALSE else {
    restr <- abs(lz) > restrict
    lz <- pmin( pmax( lz, -restrict), restrict) }
## size
  size <- i.def(size, 1)
  lppin <- par("pin")
  lratio <- 2*size * par("cin")[1] / lppin
  llwd <- i.getploption("lwd")
  lpcol <- pldata[[".pcol."]]
  if (length(lpcol)) {
    lgrpcol <- i.getploption("group.col")
    if (is.factor(lpcol)) lpcol <- as.numeric(lpcol)
    if (is.numeric(lpcol)) lpcol <- rep(lgrpcol, length=max(lpcol))[lpcol]
  } else lpcol <- i.getploption("col")[1]
  pldata$".pcol." <- colorpale(lpcol, pale=pale)
  ## for ...
  lzj <- lz[,1]
  lzj <- lzj/max(abs(lzj), na.rm = TRUE)
  lpanel <-
    function(xx, yy, indx, indy, pch, col, plab, lwd, zz, ...) {
      lusr <- par("usr")
      lfx <- lratio[1] * diff(lusr[1:2])
      lfy <- lratio[2] * diff(lusr[3:4])
      lsxz <- c(lfx * abs(zz))
      lsyz <- c(lfy * zz)
      lx <- attr(xx, "plcoord")
      if (is.null(lx)) lx <- xx
      ly <- attr(yy, "plcoord")
      if (is.null(ly)) ly <- yy
      plpoints(lx, ly, plargs=plargs, setpar=FALSE)
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
  lplxx <- i.getploption("plextext") * 2 * size
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
} ## end plres2x
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
## ============================================================================
plmframes <-
  function(mfrow=NULL, mfcol=NULL, mft=NULL, nrow=NULL, ncol=NULL, byrow=TRUE,
           oma=NULL, mar=NULL, mgp=NULL, ...)
{
## Purpose:    par(mfrow...)
  ## Author: Werner Stahel, 1994 / 2001
  lf.mft2mf <- function(mft, mfrow, din) {
    if (mfrow==0)
      mfrow <- max(1, ceiling(sqrt(mft*0.8*ldin[2]/ldin[1])) )
    ## 0.8: landscape preferred
    lmcol <- ceiling(mft/mfrow)
    c(ceiling(mft/lmcol), lmcol)
  }
  ## number of rows and cols
  lmfg <- if (length(mfrow)==2) mfrow else c(i.def(mfrow, 0), i.def(mfcol,0))
  ldin <- par("din")
  if (length(mft)) lmfg <- lf.mft2mf(mft, lmfg[1], ldin)
  ## nrow, ncol
  lnfig <- lf.mft2mf(i.getploption("mfgtotal"), 0, ldin) 
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
  mar <- rep(i.getplopt(mar, ploptions), length=4)
  if (anyNA(mar)) mar <- ifelse(is.na(mar), par("mar"), mar)
  mgp <- i.getplopt(mgp, ploptions)
  oma <- if (prod(lmfg)>1) i.getplopt(oma, ploptions) else i.def(oma, rep(0,4))
  oldpar <- if(byrow)
              par(mfrow=lmfg, oma=oma, mar=mar, mgp=mgp, ...)
            else par(mfcol=lmfg, oma=oma, mar=mar, mgp=mgp, ...)
  loldo <- ploptions(c("mar","mgp"))
  if (length(mar)) {
    ploptions(mar=mar)
    oldpar <- c(oldpar, list(mar=par("mar")))
  }
  if (length(mgp)) {
    ploptions(mgp=mgp)
    oldpar <- c(oldpar, list(mgp=par("mgp")))
  }
  invisible(
    structure(list(mfig = lmfg, mrow = if (byrow) lmfg, mcol = if(!byrow) lmfg,
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
  function (x=NULL, list=NULL, default=NULL, ploptions = NULL,
            assign=TRUE, ...)
{ ##
  lpldef <- get("ploptionsDefault", pos=1)
  lnewo <- loldo <-
    if (is.null(ploptions)) {
    if (exists(".ploptions", where=1)) get(".ploptions", pos=1)
    else  lpldef
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
  if (length(default) && u.notfalse(default)) { ## get default values
    if (u.true(default)) default <- "all"
    if (!is.character(default))
      stop("!ploptions! Argument 'default' must be of mode character or logical")
    if (default[1]=="all") return(ploptions(list=lpldef, assign=assign))
    ## resets all available components
    if (default[1]=="unset")
      return(ploptions(list=lpldef[names(lpldef)%nin%names(loldo)],
                       assign=assign) )
    if (any(default!="")) {
      llopt <- lpldef[default[default%in%names(lpldef)]]
      return( ploptions(list=llopt) )
    }
  }
  ## set options
  ## check
  largs <- check.ploption(list=largs)
  if (length(largs))  lnewo[names(largs)] <- largs
  lo <- intersect(names(largs),names(loldo))
  if (length(lo)) attr(lnewo, "old") <- loldo[lo]
  if (assign) assign(".ploptions", lnewo, pos=1)
  ## assignInMyNamespace does not work
  invisible(lnewo)
}
## ====================================================================
i.getploption <- function(opt, plo=NULL) {
  ## opt is character, plo list or NULL
  lpldef <- get("ploptionsDefault", pos=1)
  if (is.null(plo))
    plo <- get("ploptions", envir=parent.frame()) ## list in calling fn
  if (is.function(plo)) plo <- NULL
  lopt <- plo[[opt]]  
  if (is.null(lopt)||(!is.function(lopt)&&all(is.na(lopt)))) ## NULL or NA
      lopt <- ploptions(opt)
  else {lopt <- check.ploption(opt, lopt)
    if (length(lopt)) lopt <- lopt[[1]]
  }
  if (length(lopt)==0) lopt <- lpldef[[opt]]
##  names(lopt) <- opt
  lopt
}
i.getplopt <- function(opt, plo = NULL) {
  lpldef <- get("ploptionsDefault", pos=1)
  if (is.null(plo))
    plo <- get("ploptions", envir=parent.frame()) ## list in calling fn
  lnam <- as.character(substitute(opt))
  lopt <- opt 
  if (is.function(plo)) plo <- NULL
  if (is.null(lopt)||(is.atomic(lopt)&&all(is.na(lopt))))
    lopt <- plo[[lnam]]
  if (is.null(lopt)||(is.atomic(lopt)&&all(is.na(lopt)))) 
    lopt <- ploptions(lnam)
  else unlist(check.ploption(lnam, lopt))   ## check
  if (is.null(lopt)) lopt <- lpldef[[lnam]]
##  names(lopt) <- opt
  lopt
}
## -----------------------------------------------------------------------
cexSize <- function(n)  min(1.5/log10(n),2)
markextremes <- function(n) ceiling(sqrt(n)/2)/n
smoothpar <- function(n) c(min(1.8/n^0.2, 0.98), 1.5)
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
c.pchvalues <- c(0:25,33:120)
c.ltyvalues <- 1:6
## ==========================================================================
c.colors <- c("black","red","blue","darkgreen",  ##"deepskyblue3",
              "orange","purple","deepskyblue2","green3","brown",
              "pink3", "aquamarine3", "brown2") ## "darkgoldenrod3", "burlywood",
## visualize the colors, possibly to improve the definition
##- palette(c.colors)
##- plot(0:1, c(1,length(ly <-1:length(c.colors))), type="n")
##- for (lj in 1:length(ly)) lines(0:1, rep(lj,2), col=lj, lty=1)
##- c.colors <- c("black","red","blue","darkgreen","brown","orange","purple",
##-               "olivedrab", "burlywood", "violet")
c.dateticks <- data.frame(
  limit =c(30*365,3800,3.1*365,370,190, 65, 35, 10,  5,  3,  1,.33,.16,.09),
  smallunit = I(c( "y","y","m","m","m","d","d","d","h","h","h","h","M","M")),
  smallint =    c(  1 , 1 , 3 , 1 , 1 , 10, 5 , 1 , 12, 3 , 1 , 1 , 30, 15),
  bigunit =   I(c( "y","y","y","m","m","m","d","d","d","h","h","h","h","M")),
  bigint  =     c(  5 , 2 , 1 , 3 , 3 , 1 , 10, 5 , 1 ,12 , 3 , 3 , 1 , 30),
  labelunit = I(c( "y","y","y","m","m","m","d","d","d","d","h","h","h","h")),
  labelint =    c(  5 , 5 , 1 , 6 , 3 , 1 , 10, 5 , 1 , 1 , 6, 6 , 1 , 1 )
  )
## ----------------------------------------------------------------------
.ploptions <- ploptionsDefault <-
  list(
    colors = c.colors,
    linewidth = c(1,1.3,1.7,1.3,1.2,1.15), cex = 1,
    ticklength = c(-0.5, 0, 0.2, -0.2),
    ## basic
    pch = 1, cex.pch=cexSize, cex.plab=1, 
    lty=1, lwd=1, col=c.colors, lcol=c.colors,
    ## group
    group.pch=2:18, group.col=c.colors[-1], group.lty=2:6,
    group.lcol=c.colors[-1],
    ## variables
    variables.pch=1:18, variables.col=c.colors, variables.lty=1:6,
    variables.lcol=c.colors,
    ## censored
    censored.pch =  c(62, 60, 2, 23, 23, 6, 23, 23),
    ##                 >,  <, Delta, q,q, nabla, q,quadrat 
    censored.size=1.3, censored.pale = 0.3,
    ## frame
    axes = 1:2, mar=c(3.1,3.1,2.1,1.1), oma=c(2.1,2.1,2.1,2.1), mgp=c(2,0.8,0),
    panel = "plpanel", panelsep = 0.5, 
    tickintervals = c(7,3),
    date.ticks = c.dateticks, date.origin = 2000, date.format=c("y-m-d", "h:m:s"), 
    xlab = "", ylab = "",
    stamp=1, doc=TRUE, 
    mfgtotal = 30, 
    innerrange = TRUE, innerrange.factor=4, innerrange.ext=0.1,
    innerrange.function = "robrange", 
    plext=0.05, plextext=0.03,
    markextremes = markextremes, ## is a function...
    ## title (mtext)
    title.cex=c(1.2,1,1), title.cexmin=0.6, title.adj = c(0.5,0.97,0.03),
    ## grid
    grid = TRUE, grid.lty = 1, grid.lwd = 1,
    grid.col = "gray85",
    zeroline = TRUE, zeroline.lty = 1, zeroline.lwd = 1,
    zeroline.col = "gray50",
    ## refline
    refline = TRUE, 
    refline.lty = c(4,6), refline.lwd = c(1,0.7), refline.col = "darkgreen",
    ## smoothline
    smooth.lty = 2, smooth.lwd = c(2, 0.7),
    smooth.col = "blue", smooth.pale = 0.3,
    ## smooth
    smooth = TRUE, 
    smooth.function = "smoothRegr", smooth.par = smoothpar, smooth.iter = 50,
    smooth.minobs = 8, smooth.band = TRUE,
    ## bars
    bar.midpointwidth = 1, bar.lty = 1, bar.lwd = c(2,1), bar.col = "burlywood4",
    ## factors
    factor.show = "mbox", backback = TRUE, jitter = TRUE, jitter.factor = 2,
    ## time axes
    timerangelim = list(year=c(4,20), month=c(4,6), day=c(4,10), hour=4, min=4),
    ## subset
    subset.rgratio = 0.9,
    ## condquant
    condquant = TRUE, condprob.range = c(0,1), condquant.pale = c(0.5, 0.5),
    condquant.pch = c(3,4),
    ## plmatrix
    diaglabel.cex = 1.5,
    ## plregr
    functionxvalues = 51, smooth.xtrim = smoothxtrim, leveragelimit = c(0.5, 0.99),
    cookdistlines = 1:2,
    printnotices = TRUE, debug = FALSE )
.plargs <- list(ploptions=.ploptions)
  ## makes sure that  .plargs  extists when starting
## -----------------------------------------------------------------------
ploptionsCheck <-
  list(
    colors=ccl(),
    linewidth = cnr(c(0.1,5)), cex = cnr(c(0.1,5)), 
    ticklength = cnr(c(-2,2)),
    ## basic
    pch = cnv(c.pchvalues),
    cex.pch = list(cfn(),cnr(c(0.1,5))), cex.plab=cnr(c(0.1,5)),
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
    tickintervals = cnr(c(2,20)), ## date.ticks = cdf(names=...)
    date.origin = cnr(c(1900,2050)), date.format = ccv(),
    stamp=list(clg(),cnr(c(-1,2))),
    mfgtotal = cnr(c(4,100)), 
    innerrange = list(clg(),cnr()), innerrange.factor=cnr(c(0.5,10)),
    innerrange.ext=cnr(c(0,0.5)),
    plext=cnr(c(0,0.5)), plextext=cnr(c(0,0.5)),
    ## title (mtext)
    title.cex=cnr(c(0.1,5)), title.cexmin=cnr(c(0.1,5)),
    title.adj=cnr(c(-0.1,1.1)),
    ## grid
    grid = list(clg(),cnr(),cln(NA)),
    grid.lty = cnv(c.ltyvalues), grid.lwd = cnr(c(0.1,5)),
    grid.col = ccl(),
    ## zeroline
    zeroline = list(clg(),cnr()),
    zeroline.lty = cnv(c.ltyvalues), zeroline.lwd = cnr(c(0.1,5)),
    zeroline.col = ccl(),
    ## refline
    refline = list(clg(),cnr(0,2)),
    refline.lty = cnv(c.ltyvalues), refline.lwd = cnr(c(0.1,5)),
    refline.col = ccl(), 
    ## smoothline
    smooth.lty = cnv(c.ltyvalues), smooth.lwd = cnr(c(0.1,5)),
    smooth.col = ccl(), smooth.pale = cnr(c(0,1)),
    smooth = clg(), 
    smooth.function = cfn(), smooth.par = list(cfn(), cnr(c(0,2))),
    smooth.minobs = cnr(c(3,20)), smooth.band = clg(),
    ## bars
    bar.lty = cnv(c.ltyvalues), bar.lwd = cnr(c(0.1,5)), bar.col = ccl(),
    bar.midpointwidth = cnr(c(0.1,5)),
    ## factors
    factor.show = ccv(c("mbox","jitter","asis","")), backback = clg(),
    jitter = clg(),
    jitter.factor = cnr(c(0.1,5)),
    ## time axes
    ## timerangelim = c(year=365*4, month=30*12, day=4),
    ## subset
    subset.rgratio = cnr(c(0.1,1)),
    ## condquant
    condquant = clg(), condquant.pale = cnr(c(0,1)),
    condquant.pch = cnv(c.pchvalues),
    condprob.range = cnr(c(0,1)),
    ## plmatrix
    diaglabel.cex = cnr(c(0.2,10)),
    ## plregr
    functionxvalues = cnr(c(5,500)),
    smooth.xtrim = list(cfn(), cnr(c(0,0.4), na.ok=FALSE)),
    leveragelimit = cnr(c(0.1,0.99999)), cookdistlines = cnr(c(0.05, 5)),
    printnotices = clg(), debug = clg()
  )
## ==========================================================================
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
  if (length(clipped)==1) x[-li] <- clipped
  else {
    x[which(x<lrg[1])] <- clipped[1]
    x[which(x>lrg[2])] <- clipped[2]
  }
  x
}

