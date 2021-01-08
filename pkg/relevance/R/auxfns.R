
u.true <- function (x) length(x)>0 && is.logical(x) && (!is.na(x)) && all(x)
u.notfalse <-
  function (x) !(length(x)==1 && is.logical(x) && (!is.na(x)) && !x)
u.isnull <- function(x)  length(x)==0||all(is.na(x))
"%nin%" <- function (x,y) !x%in%y
i.last <- function(data, n=1) data[sign(n)*(((ldt <- length(data))-abs(n)+1):ldt)]
## u.debug <- function () u.true(rlvoptions("debug"))
DB <- function (on=TRUE) options(error=if(on) recover else NULL, warn=on)
## -----------------------------------------------------------
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
## ====================================================================
checkOption <- function(optname, value, list = NULL) {
  if (is.null(list)) list <- setNames(list(value), optname)
  lnl <- length(list)
  loptnames <- names(list)
  for (lil in seq_len(lnl)) {
    lnm <- loptnames[lil]
    lvalue <- list[[lnm]]
    lcheck <- rlvoptionsCheck[[lnm]]
    if (length(lcheck)) {
##    if (is.list(lcheck)) {
      if (!is.list(lcheck[[1]])) lcheck <- list(lcheck)
      lnopt <- length(lcheck)
      lmsg <- rep("", lnopt)
      for (lj in seq_len(lnopt)) {
        lch <- lcheck[[lj]]
        lfn <- get(lch[[1]])
        lmsg[lj] <- lmsgj <-
          switch(paste("v",length(lch)-1,sep=""), v0="", v1=lfn(lvalue),
                 v2=lfn(lvalue, lch[[2]]), v3=lfn(lvalue, lch[[2]], lch[[3]]),
                 "something wrong, please report")
        if (lmsgj=="") break
      }
##    }
      if (all(lmsg!="")) {
        warning(":check.option: argument '", lnm,
                "' not suitable. It should\n    ",
                paste(lmsg, collapse=" -- or \n  "),
                "\n  instead of (str())")
        str(lvalue)
        list[lnm] <- NULL
      }
    }
  }
  list
}
## ----------------------------------------------------------------------
cnr <- function(range=NA, na.ok=TRUE, length=NA)
  list("check.numrange", range=range, na.ok=na.ok, length=length)
cnv <- function(values=NA) list("check.numvalues", values=values)
cch <- function(values=NA) list("check.char", values=values)
ccl <- function() list("check.color", NULL)
clg <- function(na.ok=TRUE) list("check.logical", na.ok=na.ok)
cfn <- function() list("check.function", NULL)
cls <- function(names=NULL) list("check.list", names=names)
cln <- function(values=NA) list("check.listnum", values=values)

## -----------------------------------------------------------
check.numrange <- #f
  function(x, range, na.ok=TRUE, length=NA) {
  if (is.function(x)) return("be numeric (not a function)")
  if (!is.na(length)) {
    if (length > (lnx <- length(x)))
      return(paste("have length at least ",length))
  ##  if (length < lnx) 
  }
  if (na.ok && (u.isnull(x) || all(is.na(x))) ) return("")
  if (!(is.numeric(x)|is.logical(x))) return("be numeric")
  if ((!na.ok) && any(is.na(x))) return("not contain NAs")
  if (all(is.na(range))) return("")
  range <- ifelse(is.na(range), c(-Inf,Inf), range)
  if (!any(li <- x<range[1]|x>range[2], na.rm=TRUE)) return("")
  paste("be within [",paste(range, collapse=", "),"]",
        if (length(x)>1) paste("\n  violated for element(s) ",
                           paste(which(li), collapse=", ")))
}
check.numvalues <- #f
  function(x, values=NA, na.ok=TRUE) {
  if (na.ok && (u.isnull(x) || all(is.na(x))) ) return("")
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
check.char <- #f
  function(x, values=NA, na.ok=TRUE) {
  if (na.ok && (u.isnull(x) || all(is.na(x))) ) return("")
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
check.logical <- #f
  function(x, values, na.ok=TRUE) {
  if (na.ok && (u.isnull(x) || (is.atomic(x)&&all(is.na(x)))) ) return("")
  if ((is.logical(x) | (is.numeric(x))) && !all(is.na(x)) )  return("")
  "be of mode logical (or interpretable as such)"
}
check.list <- #f
  function(x, names = NULL) {
    if (is.list(x)) {
      if (length(names)) {
        if (all(names%in%names(x))) return("")
        else return(paste("be a list with names ",paste(names, collapse=", ")))
      }
    }
    "be a list"
  }
check.listnum <- #f
  function(x, values=NA, na.ok=TRUE) {
  if (is.list(x)) {
    lchk <- lapply(x, function(xx) check.numvalues(xx, values, na.ok) )
    if (all(lchk=="")) return("")
    return(paste("if a list, all components must be numeric"))
  }
  "be a list"
  }
check.function <- #f
  function(x, values, na.ok=TRUE) {
  if (is.function(x)) return("")
  if (na.ok && (u.isnull(x) || all(is.na(x))) ) return("")
  if (is.character(x))  {
    lfn <- try(get(x), silent=TRUE)
    if (inherits(lfn, "try-error"))
      return(paste("be a function or the name of an existing function.\n   '",
                   lfn, "' is not available.") )
    else return("")
  }
  "be a function or the name of an existing function."
}
