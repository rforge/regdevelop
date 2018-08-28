## Plotting functions, low and high-level
## -------------------------------------------------------------------------
pl.control <-
  function(x=NULL, y=NULL, data = NULL, transformed = TRUE,
           ## generate variables used for plotting
           psize = NULL, plab = FALSE, pch = NULL, pcol = NULL,
           cex = NULL, markextremes = NULL,
           ycol = NULL, ylty = NULL, ypch = NULL,
##           smooth = NULL,
           main = NULL, sub = ":", .subdefault = NULL,
           mar = NULL, ## needed because  markextremes  hides it
           varlabels = NULL,
           ploptions = .ploptions, .environment. = parent.frame(), ...
         )
    ## get data for plotting, collect and check arguments
    ## do preparations that are common to all plots
    ## --------------------------------------------------------------
{
  ## ---
  lcall <- sys.call()
  ## ploptions
  ploptions <- i.def(ploptions, .ploptions)
  lnmd <- setdiff(names(ploptionsDefault), names(ploptions) )
  if (length(lnmd)) ploptions <- c(ploptions, ploptionsDefault[lnmd])
  lnmdots <- setdiff(names(lcall), c("",i.argPlcontr,i.argPldata))
  if (length(lnmdots)) {
    lls <- do.call("list",as.list(lcall[lnmdots]))
    ploptions <- ploptions(list=lls, ploptions=ploptions, assign=FALSE)
  }
  ##
  if (length(ploptions$smoothlines.col)==1)
    ploptions$smoothlines.col[2] <-
      colorpale(ploptions$smoothlines.col, ploptions$smoothlines.pale)
  ## x
  lIdf <- FALSE
  lform <- NULL
  lformarg <- lynames <- lxnames <- NULL
  ldtl <- NULL
  if (length(x)) {
    if (is.atomic(x)&&is.character(x))  ## names of variables
      lxnames <- x
    else {  ## matrix or data.frame
      if (is.matrix(x)|is.atomic(x)) x <- as.data.frame(x)
      if (lIdf <- is.data.frame(x)) {
        lxnames <- names(x)
        ldtl <- if (is.null(ltit <- tit(x)))
                  as.character(attr(x, "dname")) else ltit
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
      if (is.matrix(y)|is.atomic(y)) y <- as.data.frame(y)
      if (lIdf <- is.data.frame(y)) {
        lynames <- names(y)
        x <- if (lIdf) cbind(x, y) else y
        lIdf <- TRUE
        ldtl <- NULL
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
  lformarg <- list(lxnames,lynames)
  ldtl <- NULL
  if (length(data)) {
    ldtl <- if (length(ltit <- tit(data))) ltit
    if (is.null(ldtl) & is.name(substitute(data))) ldtl <- substitute(data)
  }
  if (lIdf)  data <- if (length(data)) cbind(x, data) else x
  lvarnames <- c(lxnames, lynames)
  ## --- data
  ## ltransformed <- i.def(lcall$transformed, TRUE)
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
  if (length(lvarnames)||any(names(lcall)%in%i.argPldata)) {
    lcall <- c(list(quote(getvariables), formula=lformarg, data=data),
               as.list(lcall[intersect(largs, names(lcall))]),
               envir=.environment.)
    lcall <- as.call(lcall)
    ##        ----
    lpldata <- eval(lcall, envir=environment(lform))
    ##        ----
    if (inherits(lpldata, "pl-error"))
      stop("!pl.control! ",attr(lpldata, "message"))
    lvarnames <- attr(lpldata,"variables")
    ##
    if (length(lpcol <- lpldata$"(pcol)")) {
      lgrpcol <- i.getploption("group.col")
      if (is.factor(lpcol)) lpcol <- as.numeric(lpcol)
      if (is.numeric(lpcol)) lpcol <- lgrpcol[lpcol%%length(lgrpcol)]
      lpldata$"(pcol)" <- lpcol
    }
    if (length(lgroup <- lpldata[["(group)"]]))
      if (length(llb <- as.character(lcall$group))<=20)
        attr(lpldata[["(group)"]], "varname") <- llb
    ## --- attributes of variables
    if (length(lvarnames)) {
      lpldata[,lvarnames] <- 
        genvarattributes(lpldata[,lvarnames, drop=FALSE], ynames = lynames,
                         ycol = ycol, ylty = ylty, ypch = ypch,
                         labels = varlabels, ploptions=ploptions)
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
  if (length(lpldata) && NCOL(lpldata)==1) { ## NCOL(NULL) is 1
    lpldata <-
      cbind("(seq)"= structure(1:NROW(lpldata), label="sequence"), lpldata)
    attr(lpldata,"yvar") <- lvarnames
    attr(lpldata,"xvar") <- "(seq)"
    lvarnames <- c("(seq)", lvarnames)
  }
  ## -------------------------
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
  ## cex
##-   if (length(cex)) ploptions$cex <- cex
##-   if (length(cex.pch)) ploptions$default$cex.pch <- cex.pch
##-   if (length(cex.plab)) ploptions$default$cex.plab <- cex.plab
  ## lcp <- max(min(1,log(50)/log(lnobs))^2,0.3)
  ## select default proportion between  symbol size end text size
  ## lcpp <- 0.5 * if (all(lpch==".")) 4 else 0.7 
##  cex.pch <- i.def(cex.pch, ploptions$default$cex, valuefalse = 1 ) ## lcpp*
  ##  if (length(cex.plab)==1) cex.plab <- c(1, lcpp)*cex.plab
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
  ploptions$markextremes <- i.getplopt(markextremes)
  ploptions$mar <- i.getplopt(mar)
  ## --- condprobRange
  ploptions$condprobRange <-
    if (length(ploptions$condprobRange)==0) {
      if (lnobs>50) c(0,0) else c(0.05,0.8) }
    else c(ploptions$condprobRange,1)[1:2]
  ## ----------------------------------------------------
  ## --- smooth
  lsmgrp <- lpldata$"(smooth.group)"
  lsmgrplab <- levels(lsmgrp)
  ## smooth
  lnsm <- lnobs
  lnsmgrp <- length(unique(lsmgrp))
  if (length(lsmgrp)) lnsm <- lnobs/lnsmgrp
  ploptions$smooth <- i.def(ploptions$smooth, 0, 2, 0)
  ploptions$smooth.par <-
    i.def( ploptions$smooth.par, 5*lnsm^log10(1/2)*(1+inherits(x,"glm")) )
  ## --- main
  main <- i.def(main, "", "", "")
  sub <- i.def(sub, NULL, ":", NULL)
  ## --- more arguments
  ## reflines <- i.def(reflines, TRUE)
  ## ------------------------------------------------------------
  ## result of pl.control
  rr <- list(
    pldata = lpldata, formula = lform,
    xvar = attr(lpldata,"xvar"), yvar = attr(lpldata,"yvar"),
    nobs = lnobs, transformed = transformed,
    pch = lpch, plabel = lplabel, plab = lIplab, ##plabna = lplabna, ???
    smooth.ngroups = lnsmgrp, smooth.grouplab = lsmgrplab,
    dataLabel = ldtl, main = main, sub = sub, .subdefault = .subdefault,
    ploptions = ploptions, datetime = date()
  )
  assign(".plargs", rr, pos=1, immediate=FALSE)
  rr
} ## end of  pl.control

## ===================================================================
getvarnames <-
  function(formula, transformed=FALSE)
{ ## get  varnames 
  if (is.character(formula))
    return(list(varnames=formula, xvar=formula, yvar=NULL))
  ##    formula <- as.formula(paste("~",paste(formula,collapse="+")))
  if (is.list(formula)) formula <- formula(formula)
  if (!is.formula(formula)) stop("!getvarnames! invalid argument 'formula'")
  lyv <- NULL
  lxv <- lvnm <-
    if (transformed) rownames(attr(terms(formula[1:2]), "factors"))
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
    lvnm <- getvarnames(formula, transformed=transformed)
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
    lvmode <- sapply(variables, mode)
    if (any(li <- lvmode%nin%c("numeric","character","logical")))
      stop("!getvariables! variable(s)  ",paste(lvarnames[li],collapse=", "),
           "  have wrong mode")  ## !!! convert into 'delayed error' as below
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
      rr <- cbind(rr,data.frame(extras, check.names=FALSE))
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
i.getPlattributes <- function(prop, arg, vars, data, ploptionscomp, drop=0)
{ ## get plot property  prop  either from a direct argument  arg ,
  ## or a stored property  ploptcomp
  ## The function avoid a duplicated use of the property
  pla <- setNames(rep(NA, length(vars)), vars)
  lIvars <- length(vars)>0
  if(length(arg)) {
    lnm <- names(arg)
    if(length(lnm)&length(vars)) {
      lnm <- intersect(lnm, vars)
      if(length(lnm)) {
        pla[lnm] <- arg[lnm]
        vars <- setdiff(vars, lnm)
      }
    } else {
      if(length(arg)%in%c(1,length(pla))) {
        pla[] <- arg
        return(pla)
      }
      else 
        warning(":getPlprop: unsuitable length of (unnamed) properties\n",
                paste(arg, collapse=", "))
    }
  }
  if (length(vars)&length(data)) {
    lpr <- sapply(data, function(x) attr(x,prop))
    lpr <- unlist(lpr[sapply(lpr, length)>0])
    lnm <- intersect(names(lpr), vars)
    if(length(lnm)) {
      pla[lnm] <- lpr[lnm]
      vars <- setdiff(vars, lnm)
    }
  }
  ldef <- ploptionscomp
  if (drop>0) ldef <- ldef[-seq_len(drop)]
  if(lIvars) {
    if(length(vars)==0) return(pla)
    ldef <- rep(setdiff(ldef, pla[vars]), length=length(vars))
    pla[vars] <- ldef
  } else
    pla <- rep(ldef, length=length(pla))
  pla
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
           labels = NULL, innerrange.limits = NULL,
           ploptions = NULL, replace = FALSE)
{
  if (!is.data.frame(data))
    stop("!genvarattributes! 'data' must be a data frame")
  ## ---
  tickintervals <- i.getploption("tickintervals")
  jitter <- i.getploption("jitter")
  ljitter <- length(jitter)
  if (ljitter)
    jitter <- if (is.list(jitter)) jitter
            else setNames(rep(jitter, ncol(data)),
                          colnames(data))
  jitter.factor <- i.getploption("jitter.factor")
  ## innerrange
  lir <- i.def(i.getploption("innerrange"))
  lirglim <- innerrange.limits
  if (length(lirglim)) lir <- TRUE ## else lirglim <- lir
  if (lir)
    lirg <- if (is.list(lirglim)) lirglim
            else setNames(rep(list(lirglim), ncol(data)),
                          colnames(data))
  lirfactor <- i.getploption("innerrange.factor")
  ## labels
  llb <- NULL
  if (length(labels))
    if(!is.character(labels))
      warning("!genvarattributes! 'labels' must be of mode 'character'")
    else {
      if (length(lnm <- names(labels))==0) {
        if (length(labels)==ncol(data)) llb <- setNames(labels, names(data))
        else warning("!genvarattributes! 'labels' must be of length  ",
                     ncol(data), "  or have names")
      }  else llb <- labels
    }
  ##
  lnm <- colnames(data)
    if (anyNA(lnm)) colnames(data) <- lnm <-
                    ifelse(is.na(lnm), paste("V",1:NCOL(data), sep=""), lnm)
  ## line color and type
  if (is.null(ynames))
    ynames <- union(union(names(ycol),names(ylty)),names(ypch))
  lny <- length(ynames)
  if (lImulty <- lny>0) {
    lypch <- 
      i.getPlattributes("pch", ypch, ynames, data[,ynames, drop=FALSE],
                        ploptions$variables.pch, drop=lny>1)
    lylty <- 
      i.getPlattributes("lty", ylty, ynames, data[,ynames, drop=FALSE],
                        ploptions$variables.lty, drop=lny>1)
    lycol <-
      i.getPlattributes("col", ycol, ynames, data[,ynames, drop=FALSE],
                        ploptions$variables.col, drop=lny>1)
  }
  for (lv in lnm) {
    lvv <- data[,lv]
    lcls <- class(lvv)[1]
    attr(lvv, "varname") <- lv
    lnv <- sum(!duplicated(lvv),na.rm=TRUE)
    if (replace || is.null(attr(lvv, "nvalues")) )
      attr(lvv, "nvalues") <- lnv
    ## turn character into factor
    if (lcls=="character") lvv <- factor(lvv)
    ##    if (lv %in% lfacgen)  class(lvv) <- c(class(lvv, "usedAsFactor"))
    if (lImulty) {
      if(lv%in% names(lycol)) attr(lvv, "col") <- lycol[lv]
      if(lv%in% names(lylty)) attr(lvv, "lty") <- lylty[lv]
      if(lv%in% names(lypch)) attr(lvv, "pch") <- lypch[lv]
    }
    if (inherits(lvv, c("factor", "usedAsFactor"))) {
      ## factor
      lat <- seq_along(levels(lvv))
      if (replace || is.null(attr(lvv, "plrange")))
        attr(lvv, "plrange") <- c(0.5, max(lat)+0.5)
      if (replace || is.null(attr(lvv, "axisat")))
        attr(lvv, "axisat") <- lat
      if (replace || is.null(attr(lvv, "axislabels")))
        attr(lvv, "axislabels") <- levels(lvv)
      ## jitter
      if(replace || is.null(attr(lvv, "plcoord")) && (lij <- jitter[lv])) {
      attr(lvv, "plcoord") <-
        jitter(as.numeric(lvv), factor=jitter.factor,
               amount=if(is.numeric(lij)) lij else NULL)
      }
    } else {
      lvvv <- if(inherits(lvv,"Surv")) lvv[,1] else lvv
      lnmod <- c(0,0)
      if (is.numeric(lvvv)) {
        if (lir && u.notfalse(lirgv <- lirg[[lv]]) ) {
          lrg <- attr(lvv, "innerrange")
          lIrg <- is.null(lrg) ## new range
          lrg <- if (length(lirgv)==2 && is.numeric(lirgv)) lirgv
                 else plinnerrange(lirgv, lvvv, factor=lirfactor)
          if (replace | lIrg) {
            attr(lvv, "innerrange") <- lrg
            lpc <- plcoord(lvvv, lrg, ploptions=ploptions)
            ## attributes: avoid a level...
            lpca <- attributes(lpc)
            attributes(lvv)[names(lpca)] <- lpca
            attributes(lpc) <- NULL
            attr(lvv, "plcoord") <- lpc
            lnmod <- lpca$nmod
            if (is.null(lnmod)) lnmod <- c(0,0)
          }
        } else lrg <- range(lvvv)
        ltint <- min(max(tickintervals,3)+sum(lnmod>0),20)
        lat <- pretty(lrg, n=ltint, min.n=ltint-2)
        if (lnmod[1]) lat <- lat[lat>=lrg[1]]
        if (lnmod[2]) lat <- lat[lat<=lrg[2]]
        if (replace  || is.null(attr(lvv, "axisat")))
          attr(lvv, "axisat") <- lat ## lat[lat>lrg[1]&lat<lrg[2]]
      }
    }
    attr(lvv,"label") <- if (lv%in%names(llb)) llb[lv] else lv
    attributes(data[[lv]]) <- attributes(lvv)
  }
  data
}
## =====================================================================
plframe <-
  function(x, y, plargs=NULL, axes=1:2, xlab=NULL, ylab=NULL,
           plextext=NULL, axcol=rep(1,4), ploptions=NULL)
{
  ## -------------------
  if (is.null(ploptions)) ploptions <- plargs$ploptions
  lext <- rep(i.getploption("plext"),length=4)
  plextext <- rep(i.def(plextext, 0, i.getploption("plextext"), 0), length=4)
  axcol <- rep(i.def(axcol,1), length=4)
  lx <- if (is.data.frame(x)) x[,1] else x ## x[,2] ## lpd[,y]
  lrgx <- attr(lx,"plrange")
  if (is.null(lrgx))
    lrgx <- {
      lr <- range(as.numeric(lx), na.rm=TRUE)
      lr+c(-1,1)*diff(lr)* lext[1:2]
    }
  ly <- if (is.data.frame(y)) y[,1] else y ## x[,2] ## lpd[,y]
  lrgy <- attr(ly,"plrange")
  if (is.null(lrgy))
    lrgy <- {
      lr <- range(as.numeric(y), na.rm=TRUE)
      lr+c(-1,1)*diff(lr)*lext[3:4]
    }
  if (length(lx)*length(ly)==0)
    stop("!plframe! unsuitable argument(s) 'x' and/or 'y'")
  ## 
  lmar <- rep(i.getploption("mar"), length=4)
  lmgp <- i.getploption("mgp")
  lop <- par(mar=lmar, mgp=lmgp)
  ## on.exit(par(lop))  ## produces artifact!
  plot(lrgx + diff(lrgx)*c(-1,1)*plextext,
       lrgy + diff(lrgy)*c(-1,1)*plextext,
       xlab = "", ylab = "", type="n", axes=FALSE, xaxs="i", yaxs="i")
  ## axis labels
  lmfg <- par("mfg")
##-   if(lmar[1]>lmgp[1]+1 | lmfg[1]==lmfg[3])
##-     mtext(c(xlab, attr(lx, "label"), "")[1], side=1, line=lmgp[1], xpd=TRUE)
##-   if(lmar[2]>lmgp[1]+1 | lmfg[2]==1)
##-     mtext(c(ylab, attr(ly, "label"), "")[1], side=2, line=lmgp[1], xpd=TRUE)
  ## inner range
  lnmodx <- c(attr(lx, "nmod"),0,0)[1:2]
  lirgx <- c(attr(lx, "innerrange"), lrgx)[1:2]
  abline(v=unique(c(lirgx[lnmodx>0], lrgx)),lty=3)
  lnmody <- c(attr(ly, "nmod"),0,0)[1:2]
  lirgy <- c(attr(ly, "innerrange"), lrgy)[1:2]
  ## gridlines
  lgrl <- i.getploption("gridlines")
  ## !!! need attr("axisat") -> otherwise, generate it!
  if (u.notfalse(lgrl)) {
    if (is.atomic(lgrl)) lgrlx <- lgrly <- lgrl
    else {
      lgrlx <- lgrl[[1]]
      lgrly <- lgrl[[2]]
    }
    if (is.logical(lgrlx)) lgrlx <- attr(lx,"axisat")
    if (is.logical(lgrly)) lgrly <- attr(ly,"axisat")
    abline(v=lgrlx[lgrlx>lirgx[1]&lgrlx<=lirgx[2]],
           h=lgrly[lgrly>lirgy[1]&lgrly<=lirgy[2]],
           lty=ploptions$gridlines.lty, lwd=ploptions$gridlines.lwd,
           col=ploptions$gridlines.col)
  }
  ## zeroline
  if (length(lzl <- ploptions$zeroline)) {
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
           lty=ploptions$zeroline.lty, lwd=ploptions$zeroline.lwd,
           col=ploptions$zeroline.col)
  }
  ## bounding boxes
  abline(h=unique(c(lirgy, lrgy)), lty=3)
  lxrg <- ifelse(lnmodx>0, lirgx, lrgx)
  lyrg <- ifelse(lnmody>0, lirgy, lrgy)
  lines(lxrg[c(1,2,2,1,1)], lyrg[c(1,1,2,2,1)])
  ## axes
  if (1%in%axes) plaxis(1, lx, lmar[1]>=lmgp[2]+1 | lmfg[1]==lmfg[3], lxrg,
                        col=axcol[1])
  if (2%in%axes) plaxis(2, ly, lmar[2]>=lmgp[2]+1 | lmfg[2]==1, lyrg,
                        col=axcol[2])
  if (3%in%axes) plaxis(3, lx, lmar[3]>=lmgp[2]+1 | lmfg[1]==1, lxrg,
                        col=axcol[3])
  if (4%in%axes) plaxis(4, ly, lmar[4]>=lmgp[2]+1 | lmfg[2]==lmfg[4], lyrg,
                        col=axcol[4])
  invisible(lop)
}
## --------------------------------------------------------------------
plaxis <-
  function(axis, x, lab=TRUE, range=NULL, label=NULL, col=1, ...)
{
  range <- i.def(range, i.def(attr(x, "innerrange"), attr(x, "plrange")),
              valuefalse = range(x,na.rm=TRUE))
  lat <- attr(x,"axisat")
  label <- i.def(label, attr(x,"label"), valuefalse="")
  lmar <- par("mar")
  lmgp <- par("mgp")
  lmfg <- par("mfg")
  lIouter <- switch(axis, lmfg[1]==lmfg[3], lmfg[2]==1,
                    lmfg[1]==1, lmfg[2]==lmfg[4])
  if(lab & (lmar[axis]>lmgp[1]+1 | lIouter) )
    mtext(label, side=axis, line=lmgp[1], xpd=TRUE, col=col)
  ## tick labels
  llab <- NULL
  if (lnat <- length(lat)) {
    li <- lat>=range[1] & lat<=range[2]
    lat <- lat[li]
    if (lab && length(lat)>1) {
      llab <- attr(x,"axislabels")
      if (length(llab)) {
        if (is.numeric(llab)) {
          llab <- ifelse(lat%in%llab, as.character(lat), "")
        } else {
          if (length(llab)!=lnat) {
            warning(":plframe: incorrect number of axis labels")
            llab <- NULL
          } else llab <- llab[li]
        }
      }
    } else {
      lat <- pretty(range, 7)
      lat <- lat[lat>=range[1] & lat<=range[2]]
      llab <- if (lab) format(lat) else rep("",length(lat))
    }
  }
  axis(axis, at=lat, labels=llab, col=col, ...)
}
## ----------------------------------------------------------------------
pltitle <-
  function(main=NULL, sub=NULL, plargs=NULL, cex=NULL, cexmin=NULL, 
           side=3, line=NULL, adj=NULL, outer.margin=NULL, col="black",
           doc=NULL, sure=TRUE)
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
  ##
  lf.text <- function(text, cex, cexdef) {
    lcex <- i.def(cex, max(cexmin, min(cexdef, lfac/max(nchar(text)))),
                  valuefalse = 0 )      
    lmaxchar <- lfac/lcex
    if (any(li <- nchar(text)>lmaxchar))
      text[li] <- paste(substr(text[li], 1, lmaxchar-3),"...")
    ladj <- i.def(adj, max(minadj,0.5*(lcex>cexmin)), 0.5, minadj)
    mtext(text, side, line-lcex, cex = lcex, adj=ladj, outer = outer.margin,
          col=col)
    lcex
  }
  outer.margin <- i.def(outer.margin, par("oma")[3]>0,
                        valuefalse = FALSE)
  if ((!sure) && outer.margin && 1!=prod(par("mfg")[1:2])) return()
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
  if (length(main)) {
    lcex <- lf.text(main, cex=cex[1], cexdef=lcexdef[1])
    line <- line-lcex
  }
  if (length(lsub) && as.character(lsub)!=":") {
    lcex <- lf.text(lsub, cex=cex[2], cexdef=lcexdef[2])
    line <- line-lcex
  }
  if (line>=0 && (!is.null(doc)) && doc && length(tit(main)))
    lf.text(tit(main), cex=cex[3], cexdef=lcexdef[3])
}
## -----------------------------------------------------------------
plpoints <-
  function(x, y, type="p", plargs=NULL, plab=NULL,
           pch=NULL, col=NULL, lty=NULL, condquant = TRUE)
{
  ## lpd <- plargs$pldata
  if (is.data.frame(x)) {
    lattrx <- attributes(x[,1])
    x <- x[,1]
  } else lattrx <- attributes(x)
  x <- i.def(lattrx$plcoord, x)
  if (is.data.frame(y)) {
    lattry <- attributes(y[,1])
    y <- y[,1]
  } else lattry <- attributes(y)
  y <- i.def(lattry$plcoord, y)
  ## condquant
  condquant <- i.def(condquant, 1, 1, 0) 
  lIcq <- length(condquant)>0  ## condquant representation by bars
  lIcqx <- length(lcqx <- lattrx$condquant) >0
  lIcqy <- length(lcqy <- lattry$condquant) >0
  ##
  lnr <- length(x)
  ploptions <- i.def(plargs$ploptions, .ploptions)
  lpd <- plargs$pldata
  psize <- lpd[["(psize)"]]
  lpsize <-
    if (is.null(psize)) 1  else  sqrt(psize/median(psize, na.rm=TRUE))
  if (is.null(plab)) plab <- lpd[["(plab)"]]
  lpchdef <- if ((lIcqx|lIcqy) && condquant==0)
               ploptions$censored.pch[1] else ploptions$basic.pch
  lpchcq <-  rep(i.def(ploptions$censored.pch[-1], 3), length=3)
  if (is.null(pch))
    pch <- i.def(i.def(lpd[["(pch)"]], lattry$ypch), lpchdef[1], valuefalse="")
  pch <- rep(pch, length=lnr)
  lpcol <- i.def(col, lpd[["(pcol)"]])
  if (length(lpcol)) {
    lgrpcol <- i.getploption("group.col")
    if (is.factor(lpcol)) lpcol <- as.numeric(lpcol)
    if (is.numeric(lpcol)) lpcol <- lgrpcol[lpcol%%length(lgrpcol)]
  } else lpcol <- i.getploption("basic.col")[1]
  pcol <- rep( lpcol, length=lnr)
  lty <- i.getploption("basic.lty")
  lwd <- ploptions$linewidth[lty]
  lcol <- i.def(lpcol, i.getploption("basic.linecol"))[1]
##  lcol <- rep(lcol, length=lnr)
  lnobs <- sum(is.finite(x) & is.finite(y))
  cex <- i.getploption("cex")
  cex <- if (is.function(cex)) cex(lnobs) else i.def(cex, cexSize(lnobs))
##  cex.pch <- cex*i.getploption("default.cex")
  cex.plab <- cex*i.getploption("basic.cex.plab")
  lsplab <- rep(abs(cex.plab*lpsize), length=lnr) 
  lspch <- rep(cex*lpsize, length=lnr)
  ## condquant
  if (lIcq)
    if (condquant>0) {
      if (lIcqx) 
        x[lcqix <- lcqx[,"index"]] <- NA
      if (lIcqy) 
        y[lcqiy <- lcqy[,"index"]] <- NA
      if (lIcqx) plbars(lcqx[,1:3], y[lcqix], plargs=plargs)
      if (lIcqy) {
        if(length(lrgy <- lattry$innerrange))
          lcqy <- plcoord(lcqy[,1:3], range=lrgy, ploptions=ploptions)
        plbars(x[lcqiy], lcqy, plargs=plargs)
      }
    } else {
      pch[lcqx[,"index"]] <- lpchcq[1]
      pch[lcqy[,"index"]] <- lpchcq[2]
      pch[intersect(lcqx[,"index"],lcqy[,"index"])] <- lpchcq[3]
    }
  lIpl <- (length(plab)>0) && any(!is.na(plab))
  plab <- as.character(plab)
  ## --- plot!
  if (lIpl) {
    text(x, y, plab, cex=lsplab*lpsize, col=pcol)
    lipch <- ifelse(is.na(plab), TRUE, plab=="")
  } else lipch <- rep(TRUE,length(x))
  if (any(lipch)) {
    if (type=="l") lines(x, y, col=lcol, lty=lty, lwd=lwd)
    else {
      if (type=="b") 
        points(x, y, pch=NA, col=lcol, cex=lspch, type="b", lty=lty, lwd=lwd)
      points(x[lipch], y[lipch], pch=pch[lipch], cex=lspch[lipch],
             col=pcol[lipch])
    }
  }
}
## -----------------------------------------------------------------
plbars <-
  function(x, y, ploptions = NULL, plargs=NULL)
{
  if (is.null(ploptions)) ploptions <- plargs$ploptions
  lcol <- i.getploption("bars.col")
  llty <- rep(i.getploption("bars.lty"), length=2)
  llwd <- rep(i.getploption("bars.lwd"), length=2)
  lmpw <- i.getploption("bars.midpointwidth")*
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
  function(x, y=NULL, markextremes=NULL, plabel=NULL, plargs=.plargs)
{
  lf.pmark <- function(mprop, x) {
    lrk <- (rank(x, na.last="keep")-0.5)/sum(is.finite(x))
    lrk<mprop[1] | lrk>1-mprop[2]
  }
  ##
  plabel <- i.def(plabel, .plargs$plabel)
  if (is.null(plabel)) {
    warning(":plmark: no labels found")
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
  lmxx[is.na(lmxx)] <- 0
  lmxy[is.na(lmxy)] <- 0
  if (any(c(lmxx, lmxy)>0)) {
    li <- lf.pmark(rep(lmxx, length=2), lx)
    if (length(y)) li <- li | lf.pmark(rep(lmxy, length=2), ly)
    ifelse(li, plabel, "")
  } else plabel
}
## -----------------------------------------------------------------
plsmooth <-
  function(x, y, plargs=NULL, band=FALSE)
{
  band <- i.def(band, plargs$smooth>2, FALSE)
  lsm <- gensmooth(x, y, plargs=plargs, band=band)
  plsmoothlines(lsm, x, y, plargs=plargs)
}
## --------------------------------------------------------------
plsmoothlines <-
  function(smoothlines, x, y, plargs=NULL)
{
  if (!is.list(smoothlines)) {
    warning(":plsmoothlines: 'smoothlines' is not suitable. No smooth lines")
    return()
  }
  lrgx <- attr(x,"innerrange")
  lx <- smoothlines$x
  ploptions <- plargs$ploptions
  lgrp <- as.numeric(smoothlines$group)
  if (lInogrp <- length(lgrp)==0) {
    lgrp <- rep(1, NROW(lx))
    lcol <- i.getploption("smoothlines.col")
  } else {
    lcol <- i.getploption("group.col")
  }
  lx[lx<lrgx[1]|lx>lrgx[2]] <- NA
  lny <- NCOL(smoothlines$y)
  lngrp <- max(lgrp)
  llty <- rep(i.getploption("group.lty"), length=2*lngrp)
  ## may be a 2-vector or  a matrix of 2 rows
  llwd <- i.getploption("linewidth")
  lpale <- i.getploption("smoothlines.pale")
  for (lgr in seq_len(lngrp)) {
    lig <- which(lgrp==lgr)
    if (length(lig)) {
      lcl <- lcol[min(lgr,length(lcol))]
      if (lny>1) 
        matlines(lx[lig], smoothlines$y[lig,-1], lty=llty[2*lgr],
                 lwd=llwd[llty[2*lgr]],
                 col = colorpale(lcl, lpale))
      else
        lines(lx[lig], smoothlines$y[lig,1],
              lty=llty[2*lgr-1], lwd=llwd[llty[2*lgr-1]],
              col=lcl)  ## xxx
    if (length(lbd <- smoothlines$yband)) {
      li <- smoothlines$ybandindex[lig]
      if (any(li))
        lines(smoothlines$x[lig[li]], smoothlines$yband[lig[li]], lty=llty,
              lwd=llwd/2, col = lcl) 
      if (any(!li))
        lines(smoothlines$x[lig[!li]], smoothlines$yband[lig[!li]], lty=llty,
              lwd=llwd/2, col = lcl) 
    }
  }}
}
## -----------------------------------------------------------------
plreflines <-
  function(reflines, x=NULL, innerrange=NULL, y=NULL,
           cutrange = c(x=TRUE, y=FALSE), ploptions=NULL, plargs=NULL)
{
  ## draws a reference line (with extended range) and
  ##   band given by reflinesyw (only inner range) if requested
  lf.irna <- function(x, rg) {
    x[x<rg[1]|x>rg[2]] <- NA
    x
  }
  if (is.null(ploptions)) ploptions <- plargs$ploptions
  lrfyb <- NULL
  if (is.function(reflines))
    reflines <- reflines(y~x)
  if (is.list(reflines)&&length(reflines$coef)) reflines <- reflines$coef
  if (is.atomic(reflines)) {
    if (length(reflines)!=2) {
      warning(":plreflines: 'reflines' not suitable. No refliness")
      return()
    }
    lusr <- par("usr")
    lrfx <- seq(lusr[1],lusr[2],
                length=i.getploption("functionxvalues"))
    lrfy <- reflines[1]+reflines[2]*lrfx ## needs correction if lIxir
  } else {
    lrfx <- reflines$x
    lrfy <- reflines$y
    if (length(lrfx)==0|length(lrfx)!=NROW(lrfy)) {
      warning(":plreflines: 'reflines' not suitable. No refliness")
      return()
    }
    lrfyb <- reflines$band
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
  llty <- rep(i.def(ploptions$reflines.lty, c(4,6)), length=2)
  llwd <- rep(i.def(ploptions$reflines.lwd, c(2,1)), length=2)
  lcol <- rep(i.def(ploptions$reflines.col, "darkgreen"), length=2)
  matlines(lrfx, lrfy, lty=llty[1], col=lcol[1], lwd=llwd[1])
  if (lIrfyb)
    matlines(lrfx, as.matrix(lrfy+lrfyb),
             lty=ploptions$reflines.lty[2], lwd=ploptions$reflines.lwd[2],
             col=lcol[2])
}
## =========================================================================
plcoord <-
  function(x, range=NULL, innerrange.factor=NULL, innerrange.ext=NULL,
           plext=NULL, ploptions=NULL)
{
  ## Purpose:    values for plot with limited "inner" plot range
  ldtrg <- range(x, na.rm=TRUE)
  lirfunc <- i.getploption("innerrange.function")
  innerrange.factor <- i.getplopt(innerrange.factor)
  innerrange.ext <- i.getplopt(innerrange.ext)
  plext <- i.getplopt(plext)
  lrg <- if (length(notna(range))==0)
    lirfunc(x, fac=innerrange.factor)  else  range(range, na.rm=TRUE)
  if (length(lrg)==0) lrg <- ldtrg
  if (diff(lrg)==0) lrg <- c(-1,1)*lrg
  ## if data fits into extended inner range, then avoid inner range
  lrgext <- lrg + c(-1,1)*plext*diff(lrg) 
  if (ldtrg[1]>=lrgext[1]) lrg[1] <- ldtrg[1]
  if (ldtrg[2]<=lrgext[2]) lrg[2] <- ldtrg[2]
  ## --------- transformation
  rr <- pmax(pmin(x,lrg[2]),lrg[1])
  lxd <- x-rr
  lnmod <- c(sum(lxd<0,na.rm=TRUE),sum(lxd>0,na.rm=TRUE))
  if (sum(lnmod)>0) rr <- rr+lxd/(1+abs(lxd)/(diff(lrg)*innerrange.ext))
  ## ---------
  ## inner range must not extend beyond data
  lrg <- ifelse(lnmod, lrg, range(x, na.rm=TRUE))
  ## extend range to plotting range
  lplrg <- lrg + c(-1,1)*diff(lrg)*ifelse(lnmod>0, innerrange.ext, plext)
  ## extend inner range if there are no modified points
  lrg <- ifelse(lnmod, lrg, lplrg) 
  attr(rr,"innerrange") <- lrg
  attr(rr,"innerrange.ext") <- innerrange.ext
    ## needed for transforming further quantities
  attr(rr,"nmod") <- lnmod
  attr(rr,"plrange") <- lplrg
  class(rr) <- class(x)
  rr
}
## -------------------------------------------------------------------
i.pchcens <-
  function(plargs, condquant)
    ##  Delta, nabla, >, <, quadrat : pch= c(24, 25, 62, 60, 32)
{
  if (is.null(condquant) | !is.null(lpc <- plargs$pldata$"(pch)"))
  return(lpc)
  ##
  lpch <- plargs$ploptions$censored.pch
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
plyx <-
  function(x=NULL, y=NULL, group=NULL, data=NULL, type="p",
           xlab=NULL, ylab=NULL, 
           markextremes=0, rescale=TRUE, mar=NULL, ...)
{
  lcall <- sys.call() ## match.call modifies argument names
  lIplargs <- is.null(lcall[["plargs"]])
  plargs <-
    if (lIplargs) {
      lcall$markextremes <- markextremes
      ldtnm <- as.character(substitute(data))
      lcall$.subdefault <- if (length(ldtnm)<30) ldtnm else ""
      do.call(pl.control, as.list(lcall[-1]))
    } else  eval(lcall[["plargs"]], sys.frame(1))
  ploptions <- plargs$ploptions
  ##
  ldt <- plargs$pldata
  lnr <- NROW(ldt)
  lxnm <- attr(ldt,"xvar")
  if (is.null(lynm <- attr(ldt,"yvar"))) {
    lynm <- last(lxnm)
    lxnm <- last(lxnm, -1)
  }
  lx <- ldt[,lxnm, drop=FALSE]
  ly <- ldt[,lynm, drop=FALSE]
  ly <- lyg <- as.data.frame(transferAttributes(
    lapply(ly, function(y) i.def(attr(y,"plcoord"),y) ),
    lapply(ly, function(y) {attr(y,"plcoord") <- NULL; y}) ))
  ## why so complicated? I need this when  group  is active
  lny <- ncol(ly)
  ly1 <- ly1g <- ly[,1]
  lrgy1 <- i.def(attr(ly1, "innerrange"), attr(ly1, "plrange") ) 
  ## group
  lgroup <- ldt[["(group)"]]
  lIgrp <- length(lgroup)>0
  if(lIgrp)   {
    lgrpname <- i.def(attr(lgroup, "varname"), "group")
    lgrplab <- as.character(unique(lgroup))
  } else  lgroup <- rep(1, nrow(ly))
  if (is.factor(lgroup)) lgroup <- as.numeric(lgroup)
  lgrp <- unique(lgroup)
  ## ranges
  lrg <-
    sapply(ly, function(y) i.def(attr(y, "innerrange"), attr(y,"plrange")))
  ## inner plotting range
  lnmody <- as.matrix(sapply(ly, function(x) c(attr(x,"nmod"),0,0)[1:2]>0 ))
  lIinner <- apply(lnmody, 1, any)
  if (lny>1) {
    if (!rescale) {
      lrgy1 <- c(min(lrg[1,]),max(lrg[2,]))
      ly <- apply(ly, 2, plcoord, range=lrgy1)
        ##for (lj in 1:lny)  ly[,lj] <- plcoord(ly[,lj], lrgy1)
      ly1 <- ly[,1]
      attr(ly1,"innerrange") <- lrgy1
    }
    ## extend
    attr(ly1,"plrange") <-
      lrgy1 + diff(lrgy1)*c(-1,1)*
      ifelse(lIinner, ploptions$innerrange.ext, ploptions$plext)
  }
  ## mark extremes
  lmark <- i.getplopt(markextremes)
  lmk <- unlist(lmark)
  lImark <- length(lmk)>0 && (any(is.na(lmk)||lmk>0))
  lImark <- is.na(lImark)||lImark
  lymark <- if (lImark & lny==1) ly  ## cannot mark extremes if  lny>1
  ##
  lIsmooth <- i.def(ploptions$smooth, FALSE, TRUE, FALSE)
  lIfirst <- TRUE
  lplab <- ldt[["(plab)"]]
  lpch <- ldt[["(pch)"]]
  lIpch <- length(lpch)>0
  lpcol <- ldt[["(pcol)"]]
  lmar <- i.getploption("mar")
  loma <- i.def(ploptions$oma, par("oma"), valuefalse=0) ## ??? 
  lmgp <- i.getploption("mgp")
  lyaxcol <- 1
  if (lny>1) loma[4] <- max(loma[4],4-lmar[4])
  lmfg <- par("mfg")
  par(oma=loma)
  par(mfg=lmfg, new=FALSE)
  if (is.null(lpcol))
    lpcol <- if (lny==1) i.getploption("basic.col")[lgroup] ## !!! lny>1 NULL
  ## if (is.null(lpch)) lpch <- lgroup
  ## --- plot
  for (lj in 1:ncol(lx)) {
    lxjg <- lxj <- lx[,lj]
    if (lImark)
      lplab <- plmark(lxj, y=lymark, markextremes=lmark, plabel=plargs$plabel)
    lplabg <- lplab
    lpchg <- lpch
    for (lg in seq_along(lgrp)) {
      ## frame
      if(lny>1)  {
        lyaxcol <- lpcol <- attr(ly1, "col")
      if(!lIpch) 
        lpchg <-
          if (lny>1) attr(ly1, "pch") else i.getploption("basic.pch")[1]
      }
      lop2 <- plframe(lxj, ly1, ploptions=ploptions, xlab=xlab, ylab=ylab,
                      axcol=c(NA,lyaxcol,NA,NA))
      lop <- attr(lop2,"old")
      if (lIfirst) on.exit(par(lop))
      lIfirst <- FALSE
      lrgold <- lrgy1
      ## group
      if (lIgrp) {
        pltitle(main=NULL, sub=paste(lgrpname, lgrplab[lg], sep=": "),
                outer.margin=FALSE)  ## !!! sub!
        li <- which(lgroup==lgrp[lg])
        lxjg <- lxj[li]
        ly1g <- transferAttributes(ly1[li], ly1)
        lyg <- transferAttributes(ly[li,], ly)
        if (length(lpch)==lnr) lpchg <- lpch[li]
        if (length(lplab)==lnr) lplabg <- lplab[li]
      }
      ## smooth
      if (lIsmooth) { ## !!! lny>1
        if (lIgrp) plargs$pldata <- ldt[li,]
        plargs$ploptions$smoothlines.col <- lpcol
        plsmooth(lxjg, ly1g, plargs=plargs)
      }
      ## reflines
      if (length(lrfl <- ploptions$reflines))
        plreflines(lrfl, x=lxj, y=ly1, plargs=plargs)
      ## points
      if(!lIpch) 
        lpch <- if (lny>1) attr(ly1, "pch") else ploptions$basic.pch[1]
      plpoints(lxjg, ly1g, type=type, plargs=plargs, plab=lplab, pch=lpch,
               col=lpcol)
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
            plargs$ploptions$smoothlines.col <- lpcol
            plsmooth(lxjg, lyjg, plargs=plargs)
          }
          plpoints(lxjg, lyjg, type=type, plargs=plargs,
                   plab=lplab, pch=lpchg, col=lpcol)
          ##!!!linecolor
          lrgold <- lrgj
          if (rescale & lj==2) {
            lmfg <- par("mfg")
            plaxis(4, ly[,lj], lab=lmar[4]>=lmgp[2]+1 | lmfg[2]==lmfg[4],
                   range=lrgj, col=lpcol)
          }
        }
      }
    }
  }
  pltitle(plargs=plargs, sure=FALSE)
  stamp(sure=FALSE)
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
  }
  invisible(data)
}
## ==========================================================================
plotregr.control <-
  function(x, data = NULL, xvar = TRUE, transformed = FALSE,
           ## generate variables used for plotting
           weights = NULL, 
           ## specify some specific aspects / contents of plots
           glm.restype = "working", smresid = TRUE, 
           partial.resid = TRUE, cookdistlines = 1:2,
           leveragelim = c(0.99, 0.5), condprobRange = NULL,
           ## smooth and reflines
           reflines = TRUE,
           reflinesband = FALSE, testlevel = 0.05,
           smooth = 2, ## smoothPar=NA, smoothIter=NULL,
           smooth.sim=NULL,
           xlabs = NULL, reslabs = NULL, markextremes = NULL, 
           ## multiple frames
           mf = TRUE, mfcol = FALSE, multnrow = 0, multncol = 0,
           oma = NULL, ... )
{ ## get data for plotting, collect and check arguments
  ## do preparations that are common to plot.regr_ and plresx
  ## --------------------------------------------------------------
  lcall <- match.call()
  ## --- family
  lfam <- c(x$distrname, x$family$family)[1]
  if (is.null(lfam) || lfam=="" || is.na(lfam)) lfam <- "gaussian"
  lfamgauss <- lfam%in%c("gaussian","Gaussian")
  lfamcount <- u.true(lfam%in%c("binomial","poisson","multinom")) |
    inherits(x,"polr")
  ## --- na.action: always get full data
  ## residuals first because they fix the number of observations
  rtype <- NULL
  if (inherits(x, "glm")) rtype <- glm.restype
  lnaaction <- x$na.action
  lres <- as.data.frame(residuals(x, type=rtype, standardized=TRUE))
  lmres <- ncol(lres)
  lnobs <- sum(is.finite(lres[,1]))
  lres0 <- all( apply(lres[,1:lmres, drop=FALSE],2,
                      function(x) all(x==notna(x)[1], na.rm=TRUE ) ) )
  if (lres0)
    stop("!plot.regr/plresx! all residuals are equal -> no residual plots")
  lnr <- nrow(lres)
  ## --- ldfres 
  ldfres <- df.residual(x)
  if (is.null(ldfres)) {
    warning(":plot.regr: bug: no df of residuals. setting n-1")
    ldfres <- lnr-1
  }
  ## --- sigma
  lsigma <- x$sigma
  if (length(lsigma)==0) lsigma <- c(x$scale, summary(x)$sigma)[1]
  if (length(lsigma)==0)
    lsigma <- if (lfamcount) 0 else sqrt(apply(lres^2,2,sum)/ldfres)
  x$sigma <- lsigma
  ## --- standardized residuals
  lstres <- attr(lres, "stresiduals")
  if (length(lstres)==0) {
    llevlim <- lcall$leveragelim
    if (is.null(llevlim)) llevlim <- i.getplopt(leveragelim)
    lrs <- i.stres(x, leveragelim=llevlim)
    attributes(lres) <- c(attributes(lres), lrs)
    lstres <- lrs$stresiduals
  }
  ## --- xvar
  lform <- formula(x)
  if (length(xvar)>0) {
    if (is.character(xvar))
      lxvar <- xvar
    else {
      if (is.logical(xvar)) {
        if (is.na(xvar)||xvar) ## true
          xvar <- lform[-2]
        else  lxvar <- xvar <- NULL ## FALSE
      }
      if (length(xvar)) {
        if(!is.formula(xvar))
          stop("!plot.regr/plresx! argument 'xvar' not suitable")
        lxvar <-  ## formula
          getvarnames(update(lform, xvar), transformed=transformed)$xvar
      }
    }
  } else lxvar <- NULL
  ## residual names
  lyexpr <- deparse(lform[[2]])
  lynm <- if (nchar(lyexpr)>10) "Y" else lyexpr ## 
  lresname <- paste("res_", if (lmres>1) colnames(lres) else lynm, sep="")
  names(lres) <- lresname
  ## --- variables to be evaluated in data
  lcall <- as.list(match.call())[-1]
  ladrop <- c("xvar", "glm.restype", "smresid", "partial.resid",
              "cookdistlines", "leveragelim", "smooth", "smooth.sim",
              "reflines", "reflinesband", "testlevel", "xlabs", "reslabs",
              "mf", "mfcol", "multnrow", "multncol", "oma")
  lcall <- c(as.list(quote(pl.control)),
             as.list(lcall[setdiff(names(lcall),ladrop)]))
  if (length(xvar)) {
    lcall$x <- lxvar
    lcall$transformed <- transformed
    lcall$.subdefault <- i.form2char(lform)
    ## --- data argument
    if (length(lxvar)||any(names(lcall)%in%i.argPldata)) {
      ldata <- x$model
      if (length(ldata)==0||!(transformed & all(lxvar%in%names(ldata)))) {
        ldata <- data
        if (length(ldata)==0) {
          if (length(x$allvars))
            ldata <- x$allvars else ldata <- eval(x$call$data)
          ##-     } else {
##-       if (length(lav <- x$allvars)&&NROW(data)==NROW(lav)) 
##-        ldata <- cbind(data, lav[colnames(lav)%nin%colnames(data)])
##-       else stop("!plotregr.control! incompatible data and x$allvars")
    ##-     }
        }
        if (length(ldata)==0) {
          ldata <- x$model
          if (any(lxvj <- lxvar%nin%names(ldata)))
            stop("!plot.regr/plresx! variable(s) ",
                 paste(lxvar[lxvj], collapse=", "), " not found")
        }
        if (length(ldata)==0)  stop("!plotregr.control! No data found")
      }
  ##
      if (lnr!=nrow(ldata)) {
        if (class(lnaaction)=="omit") ldata <- ldata[-lnaaction,]
        ##    if (lnr!=nrow(ldata)) ldata <- x$model ## needs at least a warning!
        if (lnr!=nrow(ldata))
          stop("!plotregr.control! nrow of residuals and data do not agree.")
      }
    }
  } else {
    ldata <- lres ## needed to get number of observations
    lcall$x <- numeric(0) ## NULL is received as FALSE !
  } ## needed to get  nobs  in pl.control
  lcall$data <- ldata
  lftext <- i.form2char(lform)
  lcall$.subdefault <- lftext
  mode(lcall) <- "call"
  ## -------------------------------------
  plargs <- eval(lcall, parent.frame())  ## pl.control
  ## -------------------------------------
  lpldata <- plargs$pldata
  xvar <- attr(lpldata, "xvar")
  ## attributes for residuals
    ##
  lres <- genvarattributes(lres, labels = lresname, ploptions=ploptions)
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
      apply(lstres, 2, function(x) plmark(lmxres, x, plabel=plargs$plabel) )
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
  lresmahal <- x$resmahal
  if (lmres>1) 
    if (is.null(lresmahal)) lresmahal <- mahalanobis(lres,0,var(lres))
  ## --- simulated residuals
  ## when using  smooth.group , default is  0
  lnsims <- i.def(smooth.sim, 19, 19, 0)
##  if (inherits(x, c("nls", "nlm", "survreg", "polr", "coxph"))) lnsims <- 0
  if (lmres>1) lnsims <- 0 # not yet programmed for mlm
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
  reflines <- i.def(reflines, TRUE)
  reflinesband <- i.def(reflinesband, FALSE, TRUE, FALSE)
  testlevel <- i.def(i.def(testlevel, i.getploption("testlevel")), 0.05)
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
      ## resmahal = lresmahal,
      simres = lsimres, ## simstres <- lsimstres,
      yexpr = lyexpr, resname = lresname, ## stresname = lstresname,
      ##    absresname = labsresname, fitname = lfitname,
      family = lfam, famgauss = lfamgauss, famcount = lfamcount,
      na.action = lnaaction, 
      sigma=lsigma, df.residuals = ldfres, 
      ## -- return arguments
      glm.restype = glm.restype, smresid = smresid,
      partial.resid = partial.resid, cookdistlines = cookdistlines,
      leveragelim = leveragelim, ## condprobRange = condprobRange,
      ## smooth and reflines
      smooth.sim = lnsims,
      reflines = reflines,
      reflinesband = reflinesband, testlevel=testlevel,
      resplab = lresplab,
      mf=mf, oma=oma #,
    ) )
  assign(".plargs", rr, pos=1, immediate=FALSE)
  rr
}
## -----------------------------------------------------------------------
i.argPldata <- c("psize", "plab", "pch", "pcol",
                 "group", "smooth.group", "smooth.weights")
i.argPlcontr <-
  c("x", "y", "data", "transformed", "psize", "plab", "pch", "pcol", "cex",
  "markextremes", "ycol", "ylty", "ypch", "main", "sub", ".subdefault", 
  "mar", "varlabels", "ploptions", ".environment.")

i.argPlControl <- ##!!!
  c("x", "y", "data", "xvar", "transformed", i.argPldata,
    "reflines", "ploptions", 
    "main", "sub", "cex.main", "varlabels",
    ## plotregr.control
    "xvar", "weights", "glm.restype", "smresid",
    "partial.resid", "cookdistlines", "leveragelim", "condprobRange", 
    "reflinesband", "testlevel", "smooth.sim",
    "xlabs", "reslabs", 
    "mf", "mfcol", "multnrow", "multncol", "oma"
    )
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
                 "sequence", "weights", "addcomp", "smooth.legend")
  ## ----------------
  lac <- as.list(match.call())[-1]
  lextra <- setdiff(names(lac), c(argPlregr, i.argPlControl))
  if (length(lextra))
    warning(":plot.regr: argument(s)  ",paste(lextra, collapse=", "),
            "  not used")
  lcall <- c(list(quote(plotregr.control)), 
            lac[intersect(i.argPlControl, names(lac))])
##            list(fit=TRUE, hat=TRUE, stresiduals=TRUE)  )
  mode(lcall) <- "call"
  plargs <-eval(lcall, parent.frame())
  ## -------------------------------------------------------------------
  lres <-  plargs$residuals
  lmres <- plargs$rescol
  lpldata <- plargs$pldata
  ploptions <- plargs$ploptions
  lsmgrp <- lpldata$"(smooth.group)"
  lsmooth <- ploptions$smooth
  lsmpar <- ploptions$smooth.par
  lsigma <- plargs$sigma
  lresname <- plargs$resname
  lsimres <- plargs$simres
  lnsims <- if (is.null(lsimres)) 0 else ncol(lsimres)
  lreflines <- i.def(plargs$reflines, TRUE)
  ## lsimstres <- plargs$simstres
  ## from x
  lnaaction <- x$na.action
  lweights <- naresid(lnaaction, x$weights)
  lIwgt <- length(lweights)>0
  ## number of observations
  lnr <- nrow(lres)
  lnna <- !apply(is.na(lres),1,any)
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
    stop("!plot.regr! I do not know how to plot residuals of a mulitnomial regression")
  lglm <- inherits(x, "glm")
  lnnls <- !inherits(x, "nls")
##  lIcq <- length(attr(lres, "condquant"))>0
  lform <- formula(x)
  ## plot selection
  lplsel <-
    i.plotselect(plotselect, smooth=plargs$smooth, Iwgt=lIwgt,
                 mult=lmres>1,
                 famgauss=plargs$famgauss, famglm=inherits(x,"glm"),
                 famcount=plargs$famcount)
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
  if(is.null(lfit)) {
    lfit <- fitted(x)
    lfitname <- "fitted value"
  } else lfit <- naresid(lnaaction, lfit)
  lfitname <- rep(lfitname, length=lmres)
  lir <- ploptions$innerrange["fit"]
  lplopt <- if (length(lir)) ploptions(innerrange=lir) else ploptions
  lfit <-
    genvarattributes(as.data.frame(lfit), labels = lfitname, ploptions=lplopt)
  ## standardized residuals
  lstres <- attr(plargs$residuals, "stresiduals")
  if (length(lstres)==0) {
    warning(":plot.regr: I do not find standardized residuals.",
            " I will use raw residuals instead.")
    lstres <- lres
    lstrratio <- rep(1, lnr)
  } else {
    lstrratio <- i.def(attr(lres, "strratio"), rep(1,lnr))
  }
  lresplab <- plargs$resplab
  lIrpl <- length(lresplab)>0
##  lstrratio <- as.matrix(lstrratio)
  llev <- i.def(attr(lres, "leverage"), NULL, valuefalse = NULL)
  if (is.null(llev)) llev <- leverage(x) ## do not combine:
  ##           avoid evaluation of leverage(x) if not needed
##-   lsimabsres <- lsimstres <- attr(lsimres, "stres")
  ##-   if (!is.null(lsimabsres)) lsimabsres <- abs(lsimabsres)
  ## --- standardized residuals
  if (any(c(lplsel[c("absresfit","qq","leverage")],0)>0)) {
    lstres <- as.data.frame(as.matrix(lres)*lstrratio)
    lstresname <- paste("st.", lresname, sep = "")
    labsresname <- paste("|st.",lresname,"|", sep="")
    ##
    for (lj in seq_len(lmres)) {
      ## residuals from smooth
      if (plargs$smresid) {
        lrsj <- lres[,lj]
        lfsm <- gensmooth(lfit[,lj], lrsj, plargs=plargs) 
        lfsmr <- residuals(lfsm)
        if (lna <- sum(is.na(lfsmr)&!is.na(lres[,lj])))
          warning(":plot.regr: residuals from smooth have ",
                  round(100*lna/lnobs,1), " % additional NAs")
        lstres[,lj] <- lstrsj <-
          lfsmr * lstrratio[,lj] ## !!! * f(leverage of smooth)
      }
      if (length(lcqj <- attr(lres[,lj], "condquant")))
        attr(lstres[,lj], "condquant") <- lcqj[,"index",drop=FALSE] 
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
    if (is.logical(lmf)&&lmf)
      lmf <- if (lmres>1) {
               if (lmres<=4) c(NA, lmres) else lmres
             } else  c(NA, 2)
    }
  if (length(lmf)==2 && is.na(lmf[1])) {
    lmf1 <- sum(names(lplsel)%in%
                c("yfit","resfit","absresfit","qq","absresweight") )
    if(lmf1>lmf[2]+1) lmf1 <- lmf[2]
    lmf[1] <- lmf1
  }
  loma <- i.def(plargs$oma, c(2,1)*(length(lmf)>0), valuefalse=NULL)
  if (length(loma)<4) loma <- c(0,0,loma,0)[1:4]
  loldpar <- 
    if (length(lmf)&(!is.logical(lmf))) {
      if (length(lmf)==1) plmfg(mft=lmf, oma=loma) else
          plmfg(lmf[1], lmf[2], oma=loma)
    } else par(oma=loma) ## , ask=plargs$ask
  on.exit(par(loldpar), add=TRUE)
  lnewplot <- TRUE ## !!!
  ## -----------------------------------
  lplsel <- lplsel[is.na(lplsel)|lplsel>0]
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
        plframe(ly[,lj],lfit[,lj], plargs=plargs)
        plpoints(ly[,lj],lfit[,lj], plargs=plargs,
                 plab=if(lIrpl) lresplab[,lj])
        pltitle(plargs=plargs, sure=FALSE)
      }
    }
  }
  ## ---
  if(lpls=="resfit") {
    for (lj in seq_len(lmres)) {
      ## smooth
      plframe(lfit[,lj], lres[,lj], plargs=plargs)
      if (lpllevel>1)
        plsmooth(lfit[,lj], cbind(lres[,lj], lsimres),
                 plargs=plargs, band=lpllevel>2)
      if (lreflines)
        plreflines(c(x=median(lfit[,lj], na.rm=TRUE),y=-1),
                   x=lfit[,lj], y=lres[,lj], plargs=plargs)
      plpoints(lfit[,lj], lres[,lj], plargs=plargs,
               plab=if(lIrpl) lresplab[,lj])
      pltitle(plargs=plargs, sure=FALSE)
    }
  }
## ---
  if(lpls=="absresfit")
    if(length(lstres)) {
      plargsmod <- plargs
      lplext <- rep(ploptions("plext"), length=4)
      lplext[3] <- 0
      plargsmod$ploptions$plext <- lplext
      lir <- i.def(i.getploption("innerrange"))
      for (lj in seq_len(lmres)) {
        lfitj <- lfit[,lj]
        lstrj <- lstres[,lj]
        lattrj <- attributes(lstres[,lj])
        if (lir) {
          lirgj <- c(0, max(abs(as.numeric(lattrj$innerrange))))
          labssrj <-
            genvarattributes(as.data.frame(abs(lstrj)),
                             innerrange.limits=lirgj)
        }
        attr(labssrj, "condquant") <- attr(lstrj,"condquant")
        ##lattrj[["condquant"]]
        plframe(lfitj, labssrj, ylab=labsresname[lj], plargs=plargsmod)
        if (lpllevel>1) {
          lsimabsres <- if (lnsims) abs(lsimstres) else NULL
          plsmooth(lfitj, cbind(as.matrix(labssrj), lsimabsres), 
                   plargs=plargsmod, band=lpllevel>2)
          }
        plpoints(lfitj, labssrj, plargs=plargsmod, condquant=0,
               plab=if(lIrpl) lresplab[,lj]) 
        pltitle(plargs=plargs, sure=FALSE)
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
            gensmooth(lwg, cbind(labssrj, lsimabsres), plargs=plargs,
                      band=lsmooth>2)  else  NULL
        plframe(lwg, labssrj, ylab=lrlab[lj], plargs=plargsmod)
        plsmoothlines(lsm, plargs=plargs)
        plpoints(lwg, labssrj, plargs=plargsmod, condquant=0,
                 plab=if(lIrpl) lresplab[,lj])
        pltitle(plargs=plargs, sure=FALSE)
  } } }
## --- normal plot qq plot
  if(lpls=="qq") {
    ##    if (lpllevel>0)
    lnsims <- plargs$smooth.sim
    lIcq <- FALSE
    for (lj in seq_len(lmres)) {
      llr <- lstres[,lj]
      lIcqj <- length(lcq <- attr(llr, "condquant"))>0
      lIcq <- lIcq | lIcqj
      lpch <-
        if (lIcqj) i.pchcens(plargs, lcq) else ploptions$basic.pch[1]
      lxy <- qqnorm(llr, ylab = lstresname[lj], main="", type="n")
      lquart <- quantile(llr, c(0.25,0.75), na.rm=TRUE)
      plreflines(c(0, diff(lquart)/(2*qnorm(0.75))), plargs=plargs)
      if (lnsims>0) {
        lxx <- qnorm(ppoints(lnobs))
        for (lr in 1:lnsims) {
          llty <- last(ploptions$smoothlines.lty)
          lines(lxx,sort(lsimstres[,lr]), lty=llty,
                lwd=ploptions$linewidth[llty],
                col=last(ploptions$smoothlines.col))
        }
      }
      ## qq line
      li <- order(lxy$x)
      lxx <- lxy$x[li]
      lyy <- lxy$y[li]
##      lines(lxx,lyy, col=plargs$ploptions$basic.col[1])
      lpa <- plargs
      lpa$pldata <- lpldata[li,]
      plpoints(lxx, lyy, pch=lpch, plargs=lpa)  ## ??? simplify?
      if(lIcq & lj==lmres)
        legend("bottomright",
               pch=c(rep(ploptions$censored.pch,length=2)),
               legend=c("uncensored","censored"))
      pltitle(plargs=plargs, sure=FALSE)
    }
  }
## --- leverage plot. If weight are present, use "almost unweighted h"
  if(lpls=="leverage")
  if ((!is.na(lpllevel))&&lpllevel>0 && lnnls) {
    if (diff(range(llev,na.rm=TRUE))<0.001)
      warning(":plot.regr: all leverage elements equal, no leverage plot")
    else {
      llevpl <- plcoord(llev, c(0,0.5), ploptions=ploptions)
      lstres <- genvarattributes(lstres, labels = lstresname)
      ## mark extremes
      lmx <- ploptions$markextremes
      if (is.list(lmx)) lmx <- lmx[["(lev)"]]
      if (is.function(lmx)) lmx <- lmx(lnobs)
      lmx <- last(i.def(lmx, NA))
      if (lImxlev <- is.na(lmx)||lmx>0)
        lpllev <- plmark(llev, markextremes=lmx, plabel=plargs$plabel)
      ##
      llevtit <- paste("leverages", if(lIwgt) "(unweighted)")
      lcookl <- plargs$cookdistlines
      if (lIcook <- length(lcookl)>0) {
        ldfres <- df.residual(x)
        llx <- seq(min(c(lcookl,4),na.rm=TRUE)^2 *(1-ldfres/lnobs)/6,
                   max(llev,na.rm=TRUE), length=50)
        ## see formula for curves of constant Cook's distance in terms of
        ##   standardized residuals
        llrcd <- outer(sqrt((1-llx)*(lnobs-ldfres)/(llx*ldfres)),
                     c(lcookl,-lcookl)) 
      }
      for (lj in 1:lmres) {
        lstrj <- lstres[,lj]
        plframe(llevpl, lstrj, xlab=llevtit, plargs=plargs)
        if (lIcook) plreflines(list(x=llx, y=llrcd), x=llevpl, plargs=plargs)
        lplj <- if (lIrpl)  lresplab[,lj]  else  rep("", lnr)
        if (lImxlev) lplj <- ifelse(lpllev=="", lplj, lpllev)
        plpoints(llevpl, lstres[,lj], plargs=plargs, condquant=0,
                 plab=if(lIrpl|lImxlev) lplj)
        pltitle(plargs=plargs, sure=FALSE)
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
      plmatrix(lres, panel=lpanel, plargs=lpa) #main=plargs$main, pch=plargs$pch
    }
  }
  ## mahalanobis residuals
  if(lpls=="qqmult")  ## qq plot of Mahalanobis lenghts for multivariate regr
    if ((!is.na(lpllevel))&&lpllevel>0) {
      lresmahal <- x$resmahal
      if (is.null(lresmahal)) lresmahal <- mahalanobis(lres,0,var(lres))
      lxx <- sqrt(qchisq(ppoints(lresmahal),ncol(lres)))
      lor <- order(lresmahal)
      lyy <- sqrt(lresmahal[lor])
      lop <- par(mfrow=c(1,1), oma=c(0,0,2,0))
##-       plot(lxx,lyy, xlab="sqrt(Chisq.quantiles)",type="n",
##-            ylab = "Mahal.oulyingness")
      plframe(lxx,lyy, plargs=plargs,
              xlab="sqrt(Chisq.quantiles)", ylab = "Mahal.oulyingness")
      lines(lxx,lyy)
      ## !!! needs work!!!
      ##      if (ltxt) text(lxx,lyy, plab[lor]) # else points(lxx,lyy,pch=lplab[lor])
      abline(0,1,lty = ploptions$gridlines.lty, col=ploptions$gridlines.col)
      pltitle(plargs=plargs, sure=FALSE)
      stamp(sure=FALSE)
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
  plargs$mf <- FALSE ## avoid a new page
  ## plargs$ylim <- lylim  ## no need to calculate again
  if (length(plargs$xvar))
    plresx(x, data=data, resid=lres, xvar=plargs$xvar, 
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
           addcomp = FALSE, smooth.legend = FALSE, markextremes = NA,...)
## ------------------------------------------------------------
{ ## plresx
  lcall <- match.call()
  lIplargs <- is.null(lcall[["plargs"]])
  plargs <-
    if (lIplargs) {
      lac <- as.list(lcall)[-1]
      ladrop <- c("xvar", "sequence", "weights", "addcomp", "smooth.legend")
      lcall <- c(list(quote(plotregr.control)), lac[setdiff(names(lac), ladrop)])
      mode(lcall) <- "call"
      plargs <-eval(lcall, parent.frame())
    } else  eval(lcall[["plargs"]], sys.frame(1))
  ploptions <- plargs$ploptions
  ## -------------------------------------------------------------------
  if (inherits(x,"mulltinom"))
    stop("!plresx! I do not know how to plot residuals of a mulitnomial regression")
  lres <- plargs$residuals
  ldt <- plargs$pldata
##  lIcq <- inherits(lres, "condquant")
  lform <- plargs$formula
  lvars <- plargs$xvar
  lsimres <- plargs$simres
  lnsims <- if (is.null(lsimres)) 0 else ncol(lsimres)
  lInnls <- !inherits(x, "nls")
  lnr <- nrow(ldt)
  lIwgt <- length(x$weights)>0
  lnaaction <- x$na.action
  lweights <- naresid(lnaaction, x$weights)
  ##  --- sequence
  lIseq <- i.def(sequence, FALSE, TRUE, FALSE)
  if (lIseq) {
    if (length(lvars)) {
      ## is the seqence represented by any other variable?
      lseqvar <-
        if (length(lvars)>0)
          sapply(ldt[,lvars,drop=FALSE],function(x) {
            if (is.factor(x)||is.character(x)) FALSE else {
             ld <- diff(x)
             sum(ld==0)<0.1*length(x) && (all(ld<=0) | all(ld>=0)) }
           } ) else FALSE
      lIseq <- !any(lseqvar)
      if (!lIseq) warning(paste(":plresx: sequence represented by",
                                paste(lvars[lseqvar],collapse=", ")))
    }
    ldt$"(sequence)" <- structure(1:lnr, label="sequence")
    lvars <- c(lvars,"(sequence)")
  }
  ##  --- weights  as x variable
  lIweights <- i.def(weights, lIwgt, TRUE, FALSE)
  if (lIweights)
    if (!lIwgt) 
      warning(":plresx; No weights found.",
              " Cannot plot residuals against weights")   else {
    ldt[,"(weights)"] <- naresid(lnaaction, x$weights)
    lvars <- c(lvars, "(weights)")
  }
  ## ------------------
  lnvars <- length(lvars)
  if (lnvars==0) {
    warning(":plresx: I did not find any x variables")
    return() }
  ## terminmodel
  lrawv <- unlist(
    lapply(as.list(lvars),
           function(x) all.vars(as.formula(paste("~",x))) )
    )
  lvmod <- all.vars(formula(x))
  terminmodel <- lrawv%in%lvmod
  ##  reference lines
  reflines <- i.def(plargs$reflines, !inherits(x,"coxph"))
  ## type
  addcomp <- as.logical(i.def(addcomp, FALSE, TRUE, FALSE))
  ## 
  lpa <- plargs
## -----------------------------------
  ## data to be plotted
  lInnls <- !inherits(x, "nlls")
  lvcomp <- if (transformed)
              intersect(lvars, rownames(attr(terms(x),"factors")))
            else  unique(lrawv[terminmodel])
  if (any(terminmodel) && reflines && lInnls) {
    lcmp <- fitcomp(x, vars=lvcomp, transformed=transformed,
                    xfromdata=FALSE, se=i.def(plargs$reflinesband, TRUE))
    lcompx <- lcmp$x
    lcompy <- if (addcomp) lcmp$comp else -lcmp$comp
    lcompse <- lcmp$se
    if (addcomp) {
      lcompdt <-
        fitcomp(x, ldt, vars=lvcomp, transformed=transformed,
                xfromdata=TRUE)$comp
      ## !!! add to lres, careful for condq
    }
  }
  ## quantile 
  lqnt <-
    if(length(plargs$dfres)>0) qt(1-plargs$testlevel/2, plargs$dfres)
    else  qnorm(1-plargs$testlevel/2)
  ## --- smooth
  lIsmooth <- plargs$smooth
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
      if (ltin) {
        lcmpy <- lcompy[,lvx,indy]
        if (lcnt) lcmpx <- lcompx[,lvx]
      }
      lsm <- if(is.factor(xx)) FALSE else lIsmooth
      if (lsm) plsmooth(xx, yy, plargs=plargs)
      plpoints(xx, yy, plargs=plargs)
    }
    plmatrix(ldt[,lvars,drop=FALSE],lres, panel=lpanel,
             ##pch=plargs$plab, plcol=plargs$pldata$plcol,
             nrow = plargs$multnrow, ncol = plargs$multncol, plargs=plargs) 
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
    if (length(lmf)&(!is.logical(lmf))) {
      if (length(lmf)==1) plmfg(mft=lmf, oma=loma) else
          plmfg(lmf[1], lmf[2], oma=loma)
    } ## else par(oma=loma)  ## , ask=plargs$ask
  if (length(loldpar)) on.exit(par(loldpar), add=TRUE)
  ##
  lmbox <- ploptions$factor.show=="mbox"
  lrpl <- plargs$resplab
  if (lIrpl <- length(lrpl)>0) lrpl <- lrpl[,1]
  ## --- loop --- plresx
  for (lj in 1:lnvars) {
    lv <- lvars[lj]
    lvr <- lvars[lj] ## if (transformed) lvars[lj] else lrawv[lj]
    lcmpj <- terminmodel[lj] && reflines && lInnls
    lci <- if (lcmpj) lcompy[, lvr] else 0 ## !!!
    rr <- lres
    if (inherits(rr, "condquant")) rr <- rr[,1,drop=FALSE]
    if (plargs$partial.resid) 
      if (addcomp && lcmpj) {
        rr <- rr+lcompdt[, lvr]
      }
    lvv <- ldt[, lv]
    ##    lpa$pldata <- cbind(rr, lvv, ldt)
    mar <- i.def(plargs[["mar"]], c(NA, par("mar")[-1]))
    ## ---
    if (is.factor(lvv)) {
      ## factors
      ## lmar <- par("mar") ## might be adapted to nchar(levels)
      lvv <- structure(factor(lvv), label=attr(lvv, "label"))
      ll <- levels(lvv)
      lnl <- length(ll)
      lrs <- rr[[1]]
      if (lmbox) {
        ##       if (lIcq) rr <- lres[,"random"]
        plmboxes(lvv, lrs, data=ldt, mar=mar)
      } else {
        plframe(lvv,rr, plargs)
      }
      ## -
      if (lcmpj) {
        lx <- seq(along = ll)
        lcil <- lci[1:lnl]
        if ((!is.null(lrsrg <- attr(lrs,"innerrange")))&&diff(lrsrg)) {
          lcilp <- plcoord(lcil, lrsrg, innerrange.ext=plargs$innerrange.ext)
          if (any(attr(lcilp,"nmod"))) {
            liout <- lcilp!=lcil
            lrsout <- pmax(pmin(lcilp[liout], lrsrg[2]), lrsrg[1])
            segments(lx[liout]-0.4, lrsout, lx[liout]+0.4, lrsout,
                     lty=ploptions$reflines.lty, lwd = ploptions$reflines.lwd,
                     col=ploptions$reflines.col)
            lcil[liout] <- NA
          }
        }
        segments(lx-0.4, lcil, lx+0.4, lcil,
                 lty=ploptions$reflines.lty[1], lwd = ploptions$reflines.lwd[1], col=ploptions$reflines.col[1])
        if (plargs$reflinesband) {
          wid <- lqnt * lcompse[1:lnl, lv]
          lcill <- pmax(lcil-wid,lrsrg[1])
          lcilu <- pmin(lcil+wid,lrsrg[2])
          lines(c(rbind(lx-0.1,lx+0.1,lx+0.1,lx-0.1,lx-0.1,NA)),
                c(rbind(lcill,lcill,lcilu,lcilu,lcill,NA)),
                lty=ploptions$reflines.lty[2], lwd = ploptions$reflines.lwd[2], col=last(ploptions$reflines.col))
        }
      }
      if (!lmbox) plpoints(lvv,rr, plargs)
    } else { # ---
## --- continuous explanatory variable
      plframe(lvv,rr, plargs)
      if (lIsmooth) {
        if (lnsims>0)
          rr <- cbind(rr, if (addcomp && lcmpj)
                            lsimres+lcompdt[, lvr] else lsimres )
##-         if (inherits(rr, "condquant"))
##-           rr <- rr[,1]
        plsmooth(lvv, rr, plargs=lpa, band=plargs$smooth>=2) 
      }
      ## reflines
      if (lcmpj) {
        lrefx <- lcompx[,lvr]
        lrefyb <-
          if (plargs$reflinesband)
            outer(lqnt*lcompse[,lvr], c(-1,1))  else  NULL
        plreflines(list(x=lrefx, y=lci, band=lrefyb), x=lvv, y=rr,
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
  function(x, y, plargs=NULL, band=FALSE, power=1, resid="difference", ...)
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
  lpar <- ploptions$smooth.par[1]
  if (length(lpar)==0 || is.na(lpar)) lpar <- smoothpar(lnobs)
  lparband <- lpar*i.def(ploptions$smooth.par[2], 1.5, 1)
  liter <- i.def(ploptions$smooth.iter, 50)
  lnna <- apply(cbind(x,y), 1, sumna)
  x[lnna>0] <- NA
  lio <- order(as.numeric(lgrp),x)
  lio <- lio[!is.na(x[lio])]
  lxo <- x[lio] # sorted without NA
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
                     par=lpar, iterations=ploptions$smooth.iter, ...)
      if (length(lsm)) lysm[lig,j] <- lsm^(1/power)
    }
    if (band & length(lsm)) {
      lysmb <- lsm
      lsmr <- lyo[lig,1]-lsm^(1/power)
      lii <- lsmr>=0
      lsmrh <- lsmr[lii]
      ligi <- lig[lii]
      lsmh <- lsmfunc(lxo[ligi], sqrt(lsmrh),
                      weights=if (lIwgt) lwgto[ligi] else NULL,
                      par=lparband, iterations=liter)
      if (length(lsmh)) lysmb[lii] <- lsmh^2
      lii <- lsmr<=0
      lsmrl <- - lsmr[lii]
      ligi <- lig[lii]
      lsml <- lsmfunc(lxo[ligi], sqrt(lsmrl),
                      weights=if (lIwgt) lwgto[ligi] else NULL,
                      par=lparband, iterations=liter)
      if (length(lsml)) lysmb[lii] <- - lsml^2
      lysmband[lig] <- lysmb + lsm
      lsmrpos[lig] <- !lii
    }
  }
  lysmin <- matrix(NA, lnx, ncol(lyo), dimnames=list(names(x),colnames(lyo)))
  lysmin[lio,] <- lysm
  lres <- if (resid==2) ly/lysmin else ly-lysmin
  rr <- list(x = lxo, y = lysm, group = if(!lInogrp) factor(lgrpo) else NULL,
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
function(x, y=NULL, data=NULL, panel=panelDefault, 
         nrow=0, ncol=nrow, reduce=TRUE, keeppar=FALSE,
         xaxmar=NULL, yaxmar=NULL, xlabmar=NULL, ylabmar=NULL,
         oma=NULL, mar=rep(0.2,4), cex.diaglabels = 1.5, ...) 
{
## Purpose:    pairs  with different plotting characters, marks and/or colors
##             showing submatrices of the full scatterplot matrix
##             possibly on several pages
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: 23 Jul 93; minor bug-fix+comments:
  lf.axis <- function(k, axm, labm, txt, ...) {
    if (k %in% axm) axis(k)
    if (k %in% labm)
      mtext(txt, side=k, line=(0.5+1.2*(k %in% axm)), ...)
  }
  ##
  oldpar <- par(c("mfrow","mar","cex","mgp")) ##, "ask"
  lmfg <- par("mfg")
  on.exit(if (!keeppar) par(oldpar))
##---------------------- preparations --------------------------
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
    lcall$mar <- mar
    lcall$.environment. <- parent.frame()
    ##
    plargs <- eval(lcall, parent.frame())
##    plargs$main <- lmain
  } else    plargs <- eval(lcall[["plargs"]], sys.frame(1))
  ploptions <- plargs$ploptions
  ## -----------------------------------------------
  ## data and bookkeeping
  ldt <- plargs$pldata
  if (is.null(plargs$main)) plargs$main <- plargs$dataLabel
  xvar <- i.def(attr(ldt,"xvar"), names(ldt))
  nv1 <- length(xvar)
  lv1 <- lv2 <- 0
  if (is.null(y)) {
    xvar <- c(xvar, attr(ldt,"yvar"))
    nv1 <- length(xvar)
    if (reduce) { nv1 <- nv1-1; lv2 <- 1 }
    nv2 <- nv1
    ldata <- ldt[,xvar]
  } else { # cbind y to data for easier preparations
    reduce <- FALSE
    if (!is.null(dim(y))) {
      yvar <- colnames(y)
      ldata <- cbind(ldt[,xvar], as.data.frame(y))
    } else {
      yvar <- attr(ldt,"yvar")
      ldata <- ldt[,c(xvar,yvar)]
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
  ## number of panels
  lsub <- plargs$sub
  lmain <- plargs$main
  lImain <- !( (length(lmain)==0||lmain=="") & is.logical(lsub)&&lsub )
  lpin <- par("pin")
  lnmf <- if (nv1==6 && nv2==6) c(6,6) else {
    if (lpin[1]>lpin[2]) c(6,8) else c(8,6) }
  lnr <- nrow
  if (u.nuna(lnr)||lnr<1)
    lnr <- min(nv2,ceiling(nv2/((nv2-1)%/%lnmf[2]+1)))
  lnc <- ncol
  if (u.nuna(lnc)||lnc<1)
    lnc <- min(nv1,ceiling(nv1/((nv1-1)%/%lnmf[1]+1)))
  npgr <- ceiling(nv2/lnr)
  npgc <- ceiling(nv1/lnc)
  ## --- position of tick marks and axis labels, oma
  if (u.nuna(xaxmar)) xaxmar <- 1+(nv1*nv2>1)
  xaxmar <- ifelse(xaxmar>1,3,1)
  if (u.nuna(yaxmar)) yaxmar <- 2+(nv1*nv2>1)
  yaxmar <- ifelse(yaxmar>2,4,2)
  if (u.nuna(xlabmar)) xlabmar <- if (nv1*nv2==1) xaxmar else 4-xaxmar
  if (u.nuna(ylabmar)) ylabmar <- if (nv1*nv2==1) yaxmar else 6-yaxmar
  lcexmain <- i.getploption("title.cex")
  if (length(oma)!=4)
    oma <- c(2+(xaxmar==1)+(xlabmar==1), 2+(yaxmar==2)+(ylabmar==2),
             1.5+(xaxmar==3)+(xlabmar==3)+lcexmain[1]*lImain,
             2+(yaxmar==4)+(ylabmar==4))
  ## set par
  if (!keeppar) {
    par(mfrow=c(lnr,lnc))
    par(oma=oma, mar=mar, mgp=c(1,0.5,0))
  } else par(mfg=lmfg, new=FALSE)
  ## cex
  lcex <- ploptions$cex
  if (is.function(lcex)) lcex <- lcex(lnobs)
  plargs$ploptions$cex <- lcex
##
##-   ## log
##-   if (length(grep("x",log))>0) ldata[ldata[,1:nv1]<=0,1:nv1] <- NA
##-   if (length(grep("y",log))>0) ldata[ldata[,lv2+1:nv2]<=0,lv2+1:nv2] <- NA
  ##----------------- plots ----------------------------
  for (ipgr in 1:npgr) {
    lr <- (ipgr-1)*lnr ## start at row index  lr
  for (ipgc in 1:npgc) {
    lc <- (ipgc-1)*lnc
    if (reduce&&((lr+lnr)<=lc)) break
  for (jr in 1:lnr) { #-- plot row [j]
    jd2 <- lr+jr  ##  index for  y  axis
    j2 <- lv2 + jd2
    if (jd2<=nv2)  v2 <- ldata[,j2]
    lylab <- i.def(attr(v2,"label"), paste("V",j2,sep="."), valuefalse="")
    for (jc in 1:lnc) { #-- plot column  [j2-lv2] = 1:nv2
      jd1 <- lc+jc
      j1 <- lv1 + jd1
    if (jd2<=nv2 & jd1<=nv1) {
      v1 <- ldata[,j1]
      lxlab <- i.def(attr(v1,"label"), paste("V",j1,sep="."), valuefalse="")
      if (is.character(all.equal(v1,v2))) { # not diagonal
        plframe(v1, v2, plargs, axes=NULL, xlab="", ylab="")
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
        if (jr==lnr||jd2==nv2) lf.axis(1, xaxmar, xlabmar, lxlab)
        if (jc==1) lf.axis(2, yaxmar, ylabmar, lylab)
        if (jr==1) lf.axis(3, xaxmar, xlabmar, lxlab)
        if (jc==lnc||jd1==nv1) lf.axis(4, yaxmar, ylabmar, lylab)
      ## }
      ## box()
      if (lImain)   pltitle(plargs=plargs, sure=FALSE)
    }}
  }}
    ## stamp(sure=FALSE, outer.margin=TRUE) 
  }
  "plmatrix: done"
}

## ====================================================================
panelDefault <-
  function(xx, yy, indx, indy, plargs = NULL, ...)
{
  plargs <- plargs ## evaluate!
  ploptions <- plargs$ploptions
  plargs$ploptions$stamp <- FALSE
  mbox <- i.getploption("factor.show")=="mbox"
  lIsm <- i.def(ploptions$smooth, FALSE)
  if (is.character(xx)) xx <- factor(xx) ## !!! attributes!
  if (is.character(yy)) yy <- factor(yy)
  if (is.factor(xx)) { lIsm <- FALSE
    if (!is.factor(yy) & mbox) {
      plargs$pldata <- data.frame(xx=xx,yy=yy)
      plmboxes.default(xx, yy, plargs=plargs, add=TRUE, ...)
      return()
    }
  } else { 
    if (is.factor(yy) & mbox) {
      plargs$pldata <- data.frame(xx=yy, yy=xx)
      plmboxes.default(yy, xx, plargs=plargs, add=TRUE, horizontal=TRUE, ...)
      return()
    }
  }
  ## ---
  if (lIsm)  plsmooth(xx,yy, plargs=plargs)
  ## reflines
  if (length(lrfl <- i.getploption("reflines")))
    plreflines(lrfl, x=xx, y=yy, plargs=plargs)
  plpoints(xx, yy, plargs=plargs)
}
## ====================================================================
panelSmooth <-
  function(xx, yy, indx, indy, plargs, ...)
    graphics::panel.smooth(xx, yy, pch=plargs$pch, col=plargs$col,
                           cex=plargs$cex.pch, ...)
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
    return(structure(NULL, q=lq, width=lwid))
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
  ldt <- getvariables(x, data)
  l1asymbox <- length(x[[3]])>1 && as.character(x[[3]][[2]])=="1"
  if (l1asymbox)  ldt <- data.frame(ldt[,1],0,ldt[,2])
  plmboxes(ldt[,attr(ldt,"xvar")],ldt[,attr(ldt,"yvar")], data,
           .subdefault = format(x), ...)
}

plmboxes.default <-
  function(x, y, data, width=1, at=NULL, horizontal = FALSE,
    probs=NULL, outliers=TRUE, na=FALSE, # l1asymbox=FALSE,
    refline=NULL, add=FALSE, 
    xlim=NULL, ylim=NULL, axes=TRUE, xlab=NULL, ylab=NULL, 
    labelsvert=FALSE, mar=NULL,
    widthfac=NULL, minheight=NULL, colors=NULL, lwd=NULL,
    plargs=NULL,
    ...)
{
  ## Purpose:    multibox plot
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 14 Dec 2013, 23:38
  f.ylim <- function(ylm, ext)
    c((1+ext)*ylm[1]-ext*ylm[2], (1+ext)*ylm[2]-ext*ylm[1])

  lcall <- match.call()
##-   lIplargs <- is.null(lcall[["plargs"])
##-   plargs <- if (lIplargs)  { pl.control(x, y, data, ...)
  ##-             } else  eval(lcall[["plargs"]]) ##, sys.frame(1))
  if (is.null(plargs)) plargs <- pl.control(x, y, data, ...)
  ploptions <- plargs$ploptions
  ## ----------------------------------------------------------------------
  pdt <- plargs$pldata
  x <- pdt[,i.def(attr(pdt,"xvar"),1), drop=FALSE] ## may have two columns
  y <- pdt[,i.def(attr(pdt,"yvar")[1],2)]
  lhoriz <- as.logical(i.def(horizontal, FALSE, valuetrue=TRUE))
  ## widths
  lwfac <- modarg(widthfac, c(max=2, med=1.3, medmin=0.3, outl=NA, sep=0.003))
  ## colors, line widths
  lcol <- modarg(colors,
                 c(box="lightblue",med="blue",na="gray90",refline="magenta") )
  llwd <- modarg(lwd, c(med=3, range=2))
  ## data
  ## preliminary 
  lx <- factor(x[,1]) # unused levels are dropped
  llr <- ncol(x)>=2 ## asymmetrix mboxes required for binary factor
  if (llr && length(unique(x[,2]))!=2) {
    warning(":plmboxes: second x-variable must be binary. I ignore it.")
    x <- x[,1, drop=FALSE]
    llr <- FALSE
  }
  llist <- split(y,x)
  l1asymbox <- length(levels(lx))==1
  if (l1asymbox) {
    llev <- ""
    llev2 <- c(levels(x[,2]),"","")[1:2]
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
  xlab <- i.def(xlab, attr(x[,1],"label"), valuefalse="")
  if (length(xlab)>1) xlab <- xlab[2]
##  if (xlab=="1") xlab <- ""
  ylab <- i.def(ylab, attr(y, "label"), valuefalse="")
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
  na.pos <- i.def(na, c(min(ly, na.rm=TRUE)*(1-0.3)-0.3*max(ly, na.rm=TRUE)) )
  if (length(na.pos)==1)
    na.pos <- na.pos+ 0.03*diff(range(ly, na.rm=TRUE))*c(-1,1)
  lusr <- par("usr")
  ## plot range
  ##  lrg <- if (add) lusr[3:4] else attr(ly, "plrange")
  lirg <- attr(y, "innerrange")
  ljlim <- any(c(attr(y, "nmod"),0)>0)
  if (add) {
    xlim <- lusr[1:2]
    ylim <- lusr[3:4]
  } else {
    if(u.nuna(xlim)) xlim <- ## better: NAs -> default value
      range(at, na.rm=TRUE)+ max(width[c(1,length(width))])*c(-1,1)*0.5
    if (lIna <- !is.null(na.pos)) {
      lyat <- attr(y,"axisat")
      lyat <- lyat[lyat>max(na.pos)]
      if (is.null(lyat)||anyNA(lyat))
        lyat <- pretty(ly, n=i.def(ploptions$tickintervals, 7))
      attr(y, "plrange") <- range(c(attr(y, "plrange"), na.pos))
    }
    ylim <- i.def(ylim, attr(y, "plrange"))
    ## margins
    lmar1def <- i.getploption("mar")[1]
    ## the next statement defines the maximal label length as 10
    lmaxnchar <- ifelse(is.numeric(labelsvert), min(max(0,labelsvert),20), 10)
    lmar1 <- ifelse(labelsvert,
                    2 + 0.6*min(max(nchar(llev)), lmaxnchar), lmar1def)
    mar <- i.def(mar, c(lmar1, par("mar")[-1]), valuefalse=NULL)
    if (is.na(mar[1])) mar[1] <- lmar1
    if (!is.null(mar)) {
      oldpar <- par(mar=mar)
      on.exit(par(oldpar))
    }
    ## ---------------------------------
    if (lhoriz)
      plframe(y, xlim, plargs, axes=NULL, xlab=xlab, ylab=ylab) else
      plframe(xlim, y, plargs, axes=NULL, xlab=xlab, ylab=ylab)
    if (axes) {
      axis(1+lhoriz, at=at, labels=llev, las=1+2*as.logical(labelsvert))
      if(lIna && anyNA(y)) 
        mtext("NA", 2-lhoriz,1, at=mean(na.pos), las=1) ## bug prevention
      if (l1asymbox) {
        mtext(llev2[1], 1+lhoriz,1, at=0.75)
        mtext(llev2[2], 1+lhoriz,1, at=1.25)
      }
      axis(2-lhoriz, at=lyat)
    }
  } # if (!add)
  ## ---
  if (!is.null(refline)) 
    if(lhoriz)
      abline(v=refline, col=ploptions$reflines.col,
             lty=ploptions$reflines.lty, lwd=ploptions$reflines.lwd)
    else
      abline(h=refline, col=ploptions$reflines.col,
             lty=ploptions$reflines.lty, lwd=ploptions$reflines.lwd)
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
    if (length(lli )) {
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
  pltitle(plargs=plargs)
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
  lpldt <- plargs$pldata
  lxvar <- plargs$xvar
  if (length(lz)==0) lz <- plargs$pldata[,plargs$yvar]
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
  llwd <- plargs$ploptions$default$lwd[1]
  lcol <- i.def(plargs$pldata$"(pcol)", plargs$ploptions$basic.col[1])
  plargs$pldata$"(pcol)" <- colorpale(lcol, pale=pale)
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
               lwd = llwd, col=lcol, xpd=TRUE)
      ##--- mark restricted observations: ---
      if(any(restr)) {
        points((xx - lsxz)[restr], (yy - lsyz)[restr], pch= 8, mkh = 1/40,
               xpd=TRUE)
        points((xx + lsxz)[restr], (yy + lsyz)[restr], pch= 8, mkh = 1/40,
               xpd=TRUE)
      }
    } 
  ljx <- lxvar[1]
  lx <- lpldt[,ljx]
##-   llrg <- attr(lx, "plrange")
##-   attr(lx, "plrange") <- llrg + c(-1,1)*diff(llrg)* 2*lratio
  ljy <- lxvar[2]
  ly <- lpldt[,ljy]
##-   llrg <- attr(ly, "plrange")
##-   attr(ly, "plrange") <- llrg + c(-1,1)*diff(llrg)* 2*lratio
  ##--
  ##---------------
  lplxx <- i.getploption("plextext")
  plframe(lx, ly, plargs=plargs, xlab=xlab[1], ylab=ylab[1], plextext=lplxx)
  ##--- draw symbols: ---
  lpanel(lx, ly, zz=lzj, lwd=llwd, col=lcol)
  lmain <- plargs$main
  if (length(lmain)==0)
    lmain <- paste(attr(lzj, "label"), "~", plargs$formula[2])
  pltitle(lmain, sub=plargs$sub, outer.margin=FALSE, cex=plargs$cex.main)
  stamp(sure=FALSE)
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
  "plfitpairs done"
}
## ============================================================================
plmfg <-
function(mfrow=NULL, mfcol=NULL, mft=NULL, row=TRUE, oma=NULL,
                 mar=NULL, mgp=NULL, ...)
{
## Purpose:    par(mfrow...)
  ## Author: Werner Stahel, 1994 / 2001
  if (length(mfrow)==2) {
    mfcol <- mfrow[2]
    mfrow <- mfrow[1]
  }
  if (is.null(mft)) {
    if (is.null(mfrow)) mfrow <- 1
    if (is.null(mfcol)) mfcol <- 1
  } else {
    t.din <- par("din")
    if (is.null(mfrow))
      mfrow <- max(1,ceiling(sqrt(mft*t.din[2]/t.din[1])))
    mfcol <- ceiling(mft/mfrow)
    mfrow <- ceiling(mft/mfcol)
  }
  mar <- rep(i.getplopt(mar), length=4)
  if (anyNA(mar)) mar <- ifelse(is.na(mar), par("mar"), mar)
  if (length(mgp)!=3) mgp <- c(2,0.8,0)
  mfrow <- max(1,mfrow)
  mfcol <- max(1,mfcol)
  oma <- if (mfrow*mfcol>1) oma else rep(0,4)
  if (is.null(oma)) {
    oma <- i.getploption("oma")
    if (mfrow>1) oma[c(1,3)] <- oma[c(1,3)] + c(max(mgp[1]-mar[1],0), 0)
    if (mfcol>1) oma[c(2,4)] <- oma[c(2,4)] + pmax(mgp[1]-mar[c(2,4)], 0)
  }
  invisible(if(row)
            par(mfrow=c(mfrow,mfcol), oma=oma, mar=mar, mgp=mgp, ...) else
            par(mfcol=c(mfrow,mfcol), oma=oma, mar=mar, mgp=mgp, ...) )
}
## ==========================================================================
stamp <- function(sure=TRUE, outer.margin = NULL,
                  project=getOption("project"), step=getOption("step"),
                  stamp=NULL, line=NULL, ...)
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
      i.def(line, ( if(outer.margin) par("oma") else par("mar") )[4] - 0.8 )
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
  lop <- c(default, names(largs))
  ##
  if (!is.null(x)) {  ## asking for options
    if(!is.character(x)) {
      warning(":ploptions: First argument 'x' must be of mode character")
      return(NULL)
    }
    return(if(length(x)==1) loldo[[x]] else loldo[x])
  }
  ## ---
  if (!is.null(default)) { ## get default values
    if (is.logical(default)) default <- ifelse(default, "all", "")
    if (!is.character(default))
      stop("!ploptions! Argument 'default' must be of mode character or logical")
    if (default=="all") return(ploptions(list=ploptionsDefault), assign=assign)
    ## resets all available components
    if (default=="unset")
      return(ploptions(list=ploptionsDefault[names(ploptionsDefault)%nin%names(loldo)],
              assign=assign) )
    llopt <- ploptionsDefault[default[default%in%names(ploptionsDefault)]]
    return( ploptions(list=llopt) )
  }
  ## set options
  if (length(lop))  lnewo[lop] <- largs[lop]
  if (assign) assign(".ploptions", lnewo, pos=1)
  ## assignInMyNamespace does not work
  invisible( structure(lnewo, old=loldo[intersect(lop,names(loldo))] ) )
}
## ====================================================================
i.getploption <- function(opt, plo=NULL) {
  if (is.null(plo))
    plo <- get("ploptions", envir=parent.frame())
  if (is.function(plo)) plo <- NULL
  lopt <- plo[[opt]]
  if (is.null(lopt)||(!is.function(lopt)&&all(is.na(lopt)))) 
    lopt <- ploptions(opt)
  if (is.null(lopt)) lopt <- ploptionsDefault[[opt]]
##  names(lopt) <- opt
  ## check
  lopt
}
i.getplopt <- function(opt, plo=ploptions) {
  lnam <- as.character(substitute(opt))
  lopt <- opt 
  if (is.function(plo)) plo <- NULL
  if (is.null(lopt)||(!is.function(lopt)&&all(is.na(lopt))))
    lopt <- plo[[opt]]
  if (is.null(lopt)||(!is.function(lopt)&&all(is.na(lopt)))) {
    lopt <- ploptions(lnam)
    if (is.null(lopt)) lopt <- ploptionsDefault[[lnam]]
  }
##  names(lopt) <- opt
  ## check
  lopt
}
##- i.getplopt <- function(opt, name=NULL) {
##-   lanm <- i.def(name, as.character(substitute(opt)))
##-   lopt <-
##-     if (is.null(opt)||(!is.function(opt)&&all(is.na(opt)))) {
##-       la <- ploptions(lanm)
##-       if (is.null(la)) ploptionsDefault[[lanm]] else la
##-     } else opt
##-   ## check
##-   lopt
##- }
## ----------------------------------------------------------------------
c.colors <- c("black","red","blue","darkgreen","brown","orange","purple",
              "olivedrab", "burlywood", "violet")
## -----------------------------------------------------------------------
cexSize <- function(n)  min(1.5/log10(n),2)
markextremes <- function(n) ceiling(sqrt(n)/2)/n
smoothpar <- function(n) 5*n^log10(1/2)
## -----------------------------------------------------------------------
check.colors <- function(x) {x}
vn.range <- function(x, range) {
  if (class(x)!="numeric")
    return(i.warn(":argCheck: argument is not numeric"))
  else if (any(li <- x<range[1]|x>range[2]))
    return(i.warn(":argCheck: argument not in range",
            if (length(x)>1) c("for elements ",paste(li, collapse=", "))) )
  x
}
vn.values <- function(x, values) {x}

i.warn <- function(text) {
  warning("\n",text)
  NULL
}
## ----------------------------------------------------------------------
.ploptions <- ploptionsDefault <-
  list(
    colors = c.colors, linewidth = c(1,1.3,1.7,1.3,1.2),
    cex = cexSize,
    ## basic
    basic.pch = 1, basic.cex.pch=1, basic.cex.plab=1,
    basic.lty=1, basic.lwd=1, basic.col=c.colors, basic.linecol=c.colors,
    ## group
    group.pch=2:18, group.col=c.colors[-1], group.lty=2:6,
    group.linecol=c.colors[-1],
    ## variables
    variables.pch=1:18, variables.col=c.colors[-1], variables.lty=1:6,
    variables.linecol=c.colors,
    ## censored
    censored.pch = c(1,3), censored.cex=1,
    ## frame
    mar=c(3.1,3.1,2.1,1.1), oma=c(2,2,3,2), mgp=c(2,0.8,0),
    tickintervals = 7, stamp=1, 
    innerrange = TRUE, innerrange.factor=4, innerrange.ext=0.1,
    innerrange.function = robrange, 
    plext=0.05, plextext=0.03,
    markextremes = markextremes, ## is a function...
    ## title (mtext)
    title.cex=c(1.2,1,1), title.cexmin=0.7,
    ## gridlines
    gridlines = TRUE, gridlines.lty = 1, gridlines.lwd = 1,
    gridlines.col = "gray90",
    zeroline = TRUE, zeroline.lty = 1, zeroline.lwd = 1,
    zeroline.col = "gray50",
    ## reflines
    reflines.lty = c(4,6), reflines.lwd = c(1,0.7),
    reflines.col = "darkgreen",
    ## smoothlines
    smoothlines.lty = 2, smoothlines.lwd = c(2, 0.7),
    smoothlines.col = "blue", smoothlines.pale = 0.3,
    ## smooth
    smooth = TRUE, 
    smooth.function = "smoothRegr", smooth.par = NA, smooth.iter = 50,
    smooth.minobs = 8,
    ## bars
    bars.lty = 1, bars.lwd = c(2,0.5), bars.col = "burlywood4",
    stamp = TRUE,
    ## factors
    factor.show = "mbox", jitter = TRUE, jitter.factor = 2,
    ## censoring
    condprob.range = c(0,1),
    ## plotregr
    functionxvalues = 51, leveragelim = c(0.99,0.5)
  )
.plargs <- list(ploptions=.ploptions)
  ## makes sure that  .plargs  extists when starting
## -----------------------------------------------------------------------
##- ploptionsConditions <-
##-   list(
##-     colors=check.colors,
##-     linewidth = vn.range(c(0.1,5)), cex = vn.range(c(0.1,5)),
##-     cex.pch=vn.range(c(0.1,5)), cex.plab=vn.range(c(0.1,5)),
##-     ## basic
##-     basic.pch = vn.values(c.pchvalues),
##-     basic.lty=vn.values(c.ltyvalues), basic.lwd=vn.range(c(0.1,5)),
##-     basic.col=check.colors, basic.linecol=check.colors,
##-     ## group
##-     group.pch=vn.values(c.pchvalues),
##-     group.col=c.colors, group.lty=1:6,
##-     group.linecol=c.colors,
##-     ## variables
##-     variables.pch=1:18, variables.col=c.colors, variables.lty=1:6,
##-     variables.linecol=c.colors,
##-     ## censored
##-     censored.pch = c(1,3), censored.cex=1,
##-     ## frame
##-     mar=c(3.1,3.1,2.1,1.1), oma=c(2,2,3,2), mgp=c(2,0.8,0),
##-     tickintervals = 7, stamp=1, 
##-     innerrange = TRUE, innerrange.factor=4,
##-     innerrange.ext=vn.range(c(0,0.5)),
##-     plext=0.05, plextext=0.03,
##-     ## title (mtext)
##-     title.cex=c(1.2,1,1), title.cexmin=0.7,
##-     ## gridlines
##-     gridlines = TRUE, gridlines.lty = 1, gridlines.lwd = 1,
##-     gridlines.col = "gray90",
##-     ## reflines
##-     reflines.lty = c(4,6), reflines.lwd = c(1,0.7),
##-     reflines.col = "darkgreen",
##-     ## smoothlines
##-     smoothlines.lty = 2, smoothlines.lwd = c(2, 0.7),
##-     smoothlines.col = "blue", smoothlines.pale = vn.range(c(0,1)),
##-     smooth.function = "smoothRegr", smooth.minobs = 8,
##-     ## bars
##-     bars.lty = 1, bars.lwd = c(2,0.5), bars.col = "burlywood4",
##-     bars.midpointwidth = 1,
##-     ## factors
##-     factor.show = "mbox", jitter.factor = 2
##-   )
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
## ==========================================================================
vn.range <- function(range) list(class="numeric", range=range)
vn.values <- function(values) list(class="numeric", values=values)
check.colors <- function(colors) NULL
## ====================================================================
u.debug <- function() u.true(getOption("debug"))
u.true <- function(x) length(x)>0 && (!is.na(lx <- as.logical(x[1]))) && lx


## check: warning, return  ploptionsDefault
## checkarg <- function(arg, conditions, default=NULL)
