##- Functions in  general.R
##- 
##- showd
##- sumNA
##- logst
##- asinp
##- nainf.exclude
##- dropdata
##- last
##- dropNA
##- shorten
##- 
##- robrange
##- quinterpol
##- quantilew
##- 
##- modarg
##- 
##- doc
##- "doc<-"
##- tit
##- "tit<-"
##-
##- RNAMES
##- 
##- ##- Other objects in  general.R
##- c.weekdays
##- c.wkd
##- c.months
##- c.mon
##- c.quarters
## ======================================================================
showd <-
  function (data, first=3, nrow.=4, ncol.=NULL, digits=getOption("digits"))
{
  ## print some rows (and columns) of a matrix or data.frame
  if (length(data)==0) {
    cat("Null object: ", str(data), "\n")
    return()
  }
  ldoc <- getOption("doc")
  if (length(ldoc)>0 && ldoc && length(tit(data))>0) {
    cat("tit: ",tit(data),"\n")
  }
  lldim <- length(dim(data))
  if (lldim>2) stop("!showd not yet programmed for arrays")
  if (lldim>0) cat("dim: ",dim(data),"\n") else
    if (is.factor(data)) data <- as.character(data)
  ldata <- as.data.frame(data)
  l.nr <- nrow(ldata)
  l.nc <- ncol(ldata)
  if (is.null(colnames(ldata))) colnames(ldata) <- paste("c",1:l.nc,sep=".")
  ## select columns
  l.ic <- if (length(ncol.)==0) 1:l.nc  else {
    if (length(ncol.)==1) {
      if (l.nc>ncol.)
        c(seq(1,by=l.nc%/%ncol.,length=ncol.-1),l.nc) else 1:l.nc
    } else  {
      lic <- ncol.[ncol.>0&ncol<=l.nc]
      if (length(lic)>0) lic else 1:l.nc
    }
  }
  ## select rows
  if (l.nr<=nrow.+first)
    l.dc <- format(ldata[,l.ic, drop=FALSE], digits=digits)  else {
    l.ir <- c(1:first,round(seq(first,l.nr,length=nrow.+1))[-1])
    l.ir <- unique(c(last(l.ir,-1),l.nr))
    l.dc <- data.frame(u.merge(format(ldata[l.ir,l.ic]),"",after=first),
                       stringsAsFactors=FALSE)
    names(l.dc) <- colnames(ldata)[l.ic]
    lrn <- row.names(ldata)
    if (is.null(lrn)) lrn <- paste("r",1:l.nr,sep=".")
    row.names(l.dc) <- c(lrn[1:first],"...", lrn[l.ir[-(1:first)]])
  }
  ## was vector or array with only 1 column
  if (l.nc==1) {
    if (lldim>0) cat("     transposed column\n")
    row.names(l.dc) <-
      format(rbind(row.names(l.dc),l.dc[,1]),justify="right", digits=digits)[1,]
    l.dc <- t(l.dc)
  }
  print(l.dc,quote=FALSE, digits=digits)
  if (length(ldoc)&&ldoc&&length(doc(data)))
    cat("\ndoc:  ",paste(doc(data),collapse="\n  "),"\n")
  invisible(l.dc)
}
## ===================================================
sumNA <- function (object,inf=TRUE)
{
  ## Purpose:   count NAs along columns
  ## ----------------------------------------------------------------------
  ## Arguments:
  ##   object   data.frame, matrix or vector
  ##   inf      treat Inf as NA
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 10 Oct 2007, 08:18
  ff <- if(inf) {
    function(x)
    if(is.numeric(x)) sum(!is.finite(x)) else sum(is.na(x)) }
      else function(x) sum(is.na(x))
  if (is.matrix(object)) apply(object,2,ff)  else {
    if (is.list(object)) sapply(object,ff)
    else if(is.atomic(object)) ff(object)
  }
}
## ==========================================================================
logst <- function (data, calib=data, threshold=NULL, mult=1)
{
  ## Purpose:   logs of data, zeros and small values treated well
  ## -------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  3 Nov 2001, 08:22
  data <- cbind(data)
  calib <- cbind(calib)
  lncol <- ncol(calib)
  ljthr <- length(threshold)>0
  if (ljthr) {
    if (is.logical(threshold)&&threshold) 
      threshold <- attr(data, "threshold")
    if (!length(threshold)%in%c(1, lncol))
      stop("!logst! argument 'threshold' is inadequate")
    lthr <- rep(threshold, length=lncol)
    ljdt <- !is.na(lthr)
  } else {
    ljdt <- rep(TRUE, lncol)
    lthr <- rep(NA, lncol)
    for (lj in 1:lncol) {
      lcal <- calib[,lj]
      ldp <- lcal[lcal>0&!is.na(lcal)]
      if(length(ldp)==0) ljdt[lj] <- FALSE else {
        lq <- quantile(ldp,probs=c(0.25,0.75),na.rm=TRUE)
        if(lq[1]==lq[2]) lq[1] <- lq[2]/2
        lthr[lj] <- lc <- lq[1]^(1+mult)/lq[2]^mult
      }
    }
  }
  ## transform data
  for (lj in 1:lncol) {
    if (ljdt[lj]) {
      ldt <- data[,lj]
      lc <- lthr[lj]
      li <- which(ldt<lc)
      if (length(li))
        ldt[li] <- lc * 10^((ldt[li]-lc)/(lc*log(10)))
      data[,lj] <- log10(ldt)
  } }
  if (length(colnames(data)))
    lnmpd <- names(ljdt) <- names(lthr) <- colnames(data)  else
    lnmpd <- as.character(1:lncol)
  if (ncol(data)==1) data <- data[,1]
  attr(data,"threshold") <- unname(lthr)
  if (any(!ljdt)) {
    warning(":logst: no positive data",
            if(lncol>1) paste(" for variables ",lnmpd[!ljdt],
            ". These are not transformed") else ". No transformation")
    attr(data,"transformed") <- unname(ljdt)
  }
  data
}
## ===========================================================================
asinp <-
  structure( function (x) asin(sqrt(x/100))/asin(1),
            inverse = function(x) 100* sin(x*asin(1))^2
            )
## asinperc <- asinp  ## compatibility
## ==============================================================
nainf.exclude <- function (object, ...)
  ## na.omit, modified to omit also Inf and NaN values
{
  if (is.atomic(object)) {
    i <- is.finite(object)
    if (length(dim(i))) ## matrix 
      return( object[apply(i,1,all),,drop=FALSE] )
    else return( object[i] )
  }
  ## list
  n <- length(object)
    omit <- FALSE
    vars <- seq_len(n)
    for (j in vars) {
        x <- object[[j]]
        if (!is.atomic(x))
            next
##-         x <- is.na(x)
        x <- if (is.numeric(x)) !is.finite(x) else is.na(x)
        d <- dim(x)
        if (is.null(d) || length(d) != 2)
            omit <- omit | x
        else for (ii in 1:d[2]) omit <- omit | x[, ii]
    }
  xx <- object[!omit, , drop = FALSE]
  if (any(omit > 0L)) {
        temp <- seq(omit)[omit]
        names(temp) <- attr(object, "row.names")[omit]
        attr(temp, "class") <- "exclude"
        attr(xx, "na.action") <- temp
    }
    xx
  }
## ===========================================================================
dropdata <- function (data, rowid=NULL, incol="row.names", colid=NULL)
{
  ## Purpose:   drop observations from a data frame
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel
  li <- lj <- NULL
  lattr <- attributes(data)
  lattr <- lattr[is.na(match(names(lattr),
                             c("dim","dimnames","row.names","names")))]
  ln <- NROW(data)
  if (!is.null(rowid)) {
    lrn <- RNAMES(data)
    if (is.null(lrn)) lrn <- as.character(1:NROW(data))
    if (incol=="row.names")
      li <- match(as.character(rowid),lrn,nomatch=0)
    else {
      incol <- if (is.numeric(incol)) (1:ncol(data))[incol] else
      match(incol, colnames(data))
      if (is.na(incol)) stop("misspecified argument 'incol'")
      li <- match(rowid,data[,incol],nomatch=0)
    }
    if (any(li==0)) warning(":dropdata: observations",
              paste(rowid[li==0],collapse=", "),"not found")
    li <- li[li>0]
    if (!is.null(li)) {
      data <- cbind(data)[-li,]
      names(li) <- lrn[li]
    }
  }
  ## drop variables
  if (!is.null(colid)) {
    lj <- match(as.character(colid),names(data),nomatch=0)
    if (any(lj==0)) warning(":dropdata: variables  ",
              paste(colid[lj==0],collapse=", "),"  not found")
    lj <- lj[lj>0]
    if (!is.null(lj)) data <- data[,-lj,drop=FALSE]
  }
  if (length(li)==0&length(lj)==0) {
      warning(":dropdata: no data to be dropped")
      return(data)
    }
  if (length(li)) {
    if (length(li)==NROW(data)) warning(":dropobs: no observations left")
    if (length(lattr$na.action))  {
      lin <- which(naresid(lattr$na.action, 1:ln%in%li))
      names(lin) <- lrn[li]
      li <- c(lattr$na.action, lin)
    }
    class(li) <- "exclude"
    lattr$na.action <- li
  }
  attributes(data) <- c(attributes(data),lattr)
  data
}
## ===================================================
last <- function (data,n = NULL, ncol=NULL, drop=is.matrix(data))
{
  ldim <- dim(data)
  if (is.null(ldim)) {
    if (is.null(n)) n <- 1
    ldt <- length(data)
    return(data[sign(n)*((ldt-abs(n)+1):ldt)])
  }
  if (length(ldim)!=2)
    stop ("!last! not programmed for arrays of dimension >2")
  if (is.null(n)&is.null(ncol)) n <- 1
  if (is.null(n)) n <- ldim[1]
  if (is.null(ncol)) ncol <- ldim[2]
  data[sign(n)*((ldim[1]-abs(n)+1):ldim[1]),
       sign(ncol)*((ldim[2]-abs(ncol)+1):ldim[2]), drop=drop]
}
## -----------------------------------------------------------------
dropNA <- function (x, inf=TRUE) {
  if (length(dim(x))) {
    if (is.numeric(x)&inf) x[apply(is.finite(x),1,all)]
    else x[!apply(as.matrix(is.na(x)), 1, any),] ## as.matrix needed for Surv obj
  } else if (is.numeric(x)&inf) x[is.finite(x)] else x[!is.na(x)]
}
shorten <- function (x, n=50, endstring="...")
  if (nchar(x)>n)
    paste(substring(x, 1, n-nchar(endstring)),endstring, sep="") else x
## =================================================================
## ===========================================================================
robrange <-
  function (data, trim=0.2, fac=5.0, na.rm=TRUE)
{
  fac <- i.def(fac, 5, valuefalse=5)
  c.trimcoef0 <- 0.74
  c.trimcoef1 <- -0.87
  lna <- any(!is.finite(data))
  if (lna) {
    if(!na.rm) stop("!robrange! 'data' contains NAs")
    data <- data[is.finite(data)]
  }
  ln <- length(data)
  lsdexpected <- (c.trimcoef0 + c.trimcoef1*trim)
  if (is.character(data)|length(data)==0) {
    warning(":robrange: invalid data. Returning NA")
    return(c(NA,NA))
  }
  trim <- i.def(trim, 0.2)
  lmn <- mean(data,trim=trim)
  lds <- sort(abs(data-lmn))
  lnt <- ceiling((1-trim)*ln)
  if (lnt<3 | lnt==ln) {
    warning(":robrange: not enough valid data. returning ordinary range")
    lsd <- Inf } else {
    lsd <- fac*sum(lds[1:lnt]/(lnt-1)) / lsdexpected
    if (lsd==0) {
      warning(":robrange: robust range has width 0. returning ordinary range")
      lsd <- Inf }
               }
  structure(c(max(lmn-lsd,min(data)), min(lmn+lsd,max(data))),
            location=lmn, scale=lsd)
}
## =======================================================================
quinterpol <- function (x, probs = c(0.25,0.5,0.75), extend=TRUE)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 15 Nov 2014, 16:04
  lx <- x[!is.na(x)]
  ltb <- table(lx)
  ln <- length(lx)
  lnn <- length(ltb)
  ln1 <- lnn+1
  lxx <- sort(unique(lx))
  lxm <- (lxx[-1]+lxx[-lnn])/2
  lx0 <- lxx[1]
  lxn <- lxx[lnn]
  if(extend) {
    lx0 <- 2*lxx[1]-lxm[1]
    lx1 <- 2*lxx[lnn]-lxm[lnn-1]
  }
  lxe <- c(rbind(c(lx0,lxm),lxx),lx1)
  lp <- c(0,cumsum(ltb)/ln)
  lpp <- (lp[-1]+lp[-ln1])/2
  lpe <- c(rbind(lp,c(lpp,1)))  ## last element (1) is ineffective
  ld <- outer(probs,lpe,"-")
  li <- apply(ld>0,1,sum)
  lii <- 1:length(probs)
  ldd <- cbind(ld[cbind(lii,li)],ld[cbind(lii,li+1)])
  lh <- ldd[,1]/(ldd[,1]-ldd[,2])
  lxe[li]*(1-lh) + lxe[li+1]*lh
}
## =======================================================================
quantilew <- function (x, probs=c(0.25,0.5,0.75), weights=1, na.rm=FALSE)
{
  ## Purpose:   quantile with weights, crude version
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## -------------------------------------------------------------------------
  ## Author: KSK Projekt, Date: 14 Dec 1999, 12:02
  probs <- probs[!is.na(probs)]
  if (length(weights)==1) return(quantile(x, probs))
  if (length(weights)!=length(x))
    stop("!quantilew! lengths of 'x' and 'weights' must be equal")
  if (any(t.ina <- is.na(x))) {
    if (!na.rm) stop("!quantilew! NAs not allowed while 'na.rm' is FALSE")
    x <- x[!t.ina]
    weights <- weights[!t.ina]
  }
  t.i <- order(x)
  t.d <- x[t.i]
  t.wg <- cumsum(weights[t.i])/sum(weights)
  t.i1 <- apply(outer(t.wg,probs,"<"),2,sum)+1
  t.i2 <- pmin(apply(outer(t.wg,probs,"<="),2,sum)+1,length(t.d))
  (t.d[t.i1]+t.d[t.i2])/2
}
## ==========================================================================
modarg <- function (arg=NULL, default) {
  if (is.null(arg)) return(default)
  if (is.null(names(arg))) { ## unnamed
    if (length(arg)>length(default)) {
      warning(":modarg: argument too long. I use default")
      return(default)
    }
    names(arg) <- names(default)[1:length(arg)]
  }
  if (any(i <- names(arg)%nin%names(default))) {
    warning(":modarg: argument has unsuitable names: ", names(arg)[i])
    arg <- arg[!i]
  }
  if (length(arg)==0) return(default)
  if (is.list(default)) arg <- as.list(arg)
  default[names(arg)] <- arg
  default
}
## ===========================================================================
doc <- function (x) attr(x,"doc")
## ---
"doc<-" <- function (x, value)
{
  ##-- Create doc attribute or  PREpend  new doc to existing one.
  value <- as.character(value)
  attr(x, "doc") <- if (length(value)==0) NULL else
  if(value[1]=="^") value[-1] else c(value, attr(x, "doc"))
  x
}
## ---
tit <- function (x) attr(x,"tit")
## ---
"tit<-" <- function (x, value) ## ! argument must be 'value'. demanded by attr
{
  attr(x, "tit") <- value
  x
}
## ---
## ===========================================================================
getvarnames <-
  function(formula, data = NULL, transformed=FALSE)
{ ## get  varnames 
  if (is.character(formula))
    return(list(varnames=formula, xvar=formula, yvar=NULL))
  if (is.null(formula)) return(list(varnames=NULL, xvar=NULL, yvar=NULL))
  ##    formula <- as.formula(paste("~",paste(formula,collapse="+")))
  ## if (is.list(formula))
  formula <- formula(formula)
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
    if ("." %in% lxv) {
      if (length(data)==0)
        stop("!getvarnames! '.' in formula and no 'data'")
      lform <- formula(terms(formula, data=data))
      lxv <- if (transformed) rownames(attr(terms(lform[-2]), "factors"))
             else all.vars(lform[-2])
    }
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
    lvnm <- getvarnames(formula[[1]], data=data, transformed=transformed)
    lvnmy <- if (length(lfoy <- formula[[2]]))
               getvarnames(lfoy, data=data, transformed=transformed)  else NULL
    lvarnames <- unique(c(lvnm$varnames, lvnmy$varnames))
##    lenv <- environment(formula[[1]])
  } else {
    lvnm <- if (length(formula))
              getvarnames(formula, data=data, transformed=transformed)
    lvnmy <- NULL
    lvarnames <- lvnm$varnames
##    lenv <- environment(formula)
  }
  ## data
  if (length(data)==0) 
    data <- environment(formula)
  else {
    if (!is.data.frame(data) && !is.environment(data) && 
           !is.null(attr(data, "class"))) 
      data <- as.data.frame(data)
    else {
      if (is.array(data)) 
        stop("!getvariables! 'data' must be a data.frame, not a matrix or an array")
    }
  }##
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
    variables <- try(eval(inp, data, lenv), silent=TRUE)
    if (class(variables)=="try-error") {
      lvnmiss <- setdiff(lvn, names(data))
      stop(sub("object", "!getvariables! variable (or object)",
               attr(variables, "condition")$message),
           if (length(lvnmiss)>1)
             paste(". \n    All of ",
                   paste(lvnmiss, collapse=", "), "may be unavailable.")
           )
    }
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
           "  not found (or have wrong mode)")  ## !!! convert into 'delayed error' as below
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
i.getvarattribute <- function (attr, value=NULL, data, ploptionscomp, drop=0)
{ ## get plot property  attr  either from a direct argument  value ,
  ## or an attribute of a variable in data, or a stored property  ploptionscomp
  ## The function avoids a duplicated use of the property
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
    warning(":i.setvarattribute: name  ",paste(lnm[ljvar==0], collapse=", "),
            "  of 'value' not in names(data). Attribute '",attr,
            "' not set.")
  }
  lvar <- var[ljvar]
  if (length(lvar))
    for (lv in lvar) {
      if (attr=="plrange") {
        lattr <- i.adjustticks(value[[lv]], data[,lv])
        attributes(data[,lv])[names(lattr)] <- lattr
      } else  attr(data[,lv], attr) <- value[[lv]]
    }
  data
}
## -------------------------------------------------------------------
i.adjustticks <- function (range, data){
  if (is.null(lta <- attr(data, "ticksat")))
    return( list(plrange=range) )
  lrgtl <- diff(range(lta))
  lrgl <- diff(range)
  if (abs(lrgtl-lrgl)>0.4*lrgl) ## do not use given ticks
    return(list( plrange=range,
                ticksat=NULL, ticklabelsat=NULL, ticklabels=NULL ))
  if (any((ltd <- diff(lta))!= (ltad <- ltd[1]))) ## unequal distances ->
    ## do not know how to extend
    return( list(plrange=range) )
  lnint <- lrgl/ltad
  ltan <- clipat(lta[1]+ltad*seq(-lnint, lnint), range)
  if (is.null(ltal <- attr(data, "ticklabelsat")))
    return(list( plrange=range, ticksat=ltan) )
  if (any((ltd <- diff(ltal))!= (ltad <- ltd[1]))) ## unequal distances 
    return( list(plrange=range, ticksat=ltan) )
  lnint <- lrgl/ltad
  ltaln <- clipat(ltal[1]+ltad*seq(-lnint, lnint), range)
  ## ticklabels are too difficult to extend
  list(plrange=range, ticksat=ltan, ticklabelsat=ltaln, ticklabels=NULL)
}
## ============================================================
## additional useful functions
## ===========================================================================
## auxiliary functions
## ============================================================
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
RNAMES <- function (x) if (!is.null(dim(x))) row.names(x) else names(x)

is.formula <- function (object)
  length(class(object))>0 && inherits(object, "formula")
u.true <- function (x) length(x)>0 && is.logical(x) && (!is.na(x)) && all(x)

nafalse <- function (x) if (is.null(x)) FALSE else ifelse(is.na(x), FALSE, x)
u.nuna <- function (x)  length(x)==0 || (is.atomic(x)&&any(is.na(x)))
"%nin%" <- function (x,y) !x%in%y
u.notfalse <-
  function (x) !(length(x)==1 && is.logical(x) && (!is.na(x)) && !x)
u.debug <- function () u.true(ploptions("debug"))
u.asformula <- function (x) {
  if (is.formula(x)) return(x)
  if (is.character(x)) as.formula(paste("~",paste(x,collapse="+")))
  else stop("!u.asformula! inadequate argument")
}
## ----------------------------------------------------------------
warn <- function ()
  table(paste(names(lw <- warnings()),"@",substr(unlist(lw),1,10)))

getmeth <- function (fn,mt)  getS3method(as.character(substitute(fn)),
                                        as.character(substitute(mt)))

DB <- function (on=TRUE) options(error=if(on) recover else NULL, warn=on)

## ---------------------------------------------------------------------
transferAttributes <- function (x, xbefore, except=NULL)
{
  lattr <- attributes(xbefore)
  latnm <- c("class", "names", "dim", "dimnames", "row.names", except)
  attributes(x) <-
    c(attributes(x)[latnm],lattr[setdiff(names(lattr), latnm)])
  if (is.list(x) & is.list(xbefore) && all(names(x)==names(xbefore)))
    for (lnm in names(x))
      x[[lnm]] <- transferAttributes(x[[lnm]], xbefore[[lnm]])
  x
}
## ==========================================================================
u.merge <- function (dd1, dd2 = NA, which=NULL, after=NULL,
                    length=NULL, names=NULL)
{
## Purpose:   merge two vectors or expand a vector by NA s
## -------------------------------------------------------------------------
## Arguments:
##   dd1      first vector or matrix or data.frame (?),
##   dd2      second vector, ...
##   which    is T for indices for which first vector is used
##   after    elements of  dd2  will be inserted after "after" in dd1
##   length   length of the result (will be expanded if necessary)
##   names    names of the result (if length is adequate)
## -------------------------------------------------------------------------
## Author: Werner Stahel, Date: 11 Mar 93, 13:50, and later
  llen <- length
  n1 <- length(dd1)
  nc1 <- ncol(dd1)
  nc2 <- ncol(dd2)
  if (length(nc1)>0) {
    n1 <- nrow(dd1)
    if (!( length(dd2)==1 || is.null(nc2) || nc2==nc1 ))
      stop("unsuitable second argument")
    }
## --- generate  which  vector for all cases
  if (length(which)==0) {
## - after  specified
      if (length(after)==0) stop("specify either  which  or  after")
      if (is.logical(after))  after <- which(after)
      wh <- rep(TRUE,n1+length(after))
      wh[after+1:length(after)] <- FALSE }
  else {
## - which  specified
    if(is.logical(which)) wh <- which
    else {
      if (length(llen)==0)  llen <- n1+length(which)
        wh <- rep(TRUE, llen)
        wh[which] <- FALSE }
  }
## --- merge
  nn <- length(wh)
  n2 <- nn-n1
  if (!(is.null(names)|length(names)==nn))
    warning("argument  names  not used (unsuitable length)")
  if (length(nc1)>0) {
    if (!(length(dd2)==1 || NROW(dd2)==n2))
      stop("unsuitable number of rows")
    rr <- matrix(NA,nn,nc1)
    rr[wh,] <- as.matrix(dd1)
    rr[!wh,] <- if (is.data.frame(dd2)) as.matrix(dd2) else dd2
##-     if (length(names)>0) row.names(rr) <- names else {
##-       if (length(lrn1 <- row.names(dd1))>0)
  }
  else {
    rr <- rep(NA,nn)
    rr[wh] <- dd1
    rr[!wh] <- dd2
    if (length(names)>0) names(rr) <- names
  }
  rr
}
## ==========================================================================
i.naresid.exclude <- function (omit, x, ...) 
{ ## can be replaced by  naresid  when naresid allows for data.frames
  if(length(omit) == 0) return(x) ##!!!
  if(!is.numeric(omit)) 
    stop("invalid argument 'omit'")
  if (is.null(x)) 
    return(x)
  if (is.matrix(x)|is.data.frame(x)) { ##!!!
    n <- nrow(x)
    temp <- rownames(x)
    keep <- rep.int(NA, n + length(omit))
    keep[-omit] <- 1L:n
    x <- x[keep, , drop = FALSE]
    if (length(temp)) {
      keep[omit] <- n+1:length(omit)
      rown <- c(temp, names(omit))[keep]
      if (!u.notfalse(any(duplicated(rown)))) rownames(x) <- rown
    }
  }
  else {
    n <- length(x)
    keep <- rep.int(NA, n + length(omit))
    keep[-omit] <- 1L:n
    x <- x[keep]
    temp <- names(x)
    if (length(temp)) {
      temp[omit] <- names(omit)
      names(x) <- temp
    }
  }
  x
}
## ==========================================================================
c.weekdays <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
            "Saturday", "Sunday")
c.wkd <- substring(c.weekdays,1,3)
c.months <- c("January", "February", "March", "April", "May", "June",
        "July", "August", "September", "October", "November",
        "December")
c.mon <- substring(c.months,1,3)
c.quarters <- c("I","II","III","IV")
## ==========================================================================
