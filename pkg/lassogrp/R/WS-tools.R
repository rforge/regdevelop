#### WS-tools : Some of Werner Stahel's R Tools, not related much to the package:
#### ========

descr <- function(x) attr(x,"descr")
## ---
"descr<-" <- function(x, value)
{
  ##-- Create descr attribute or  PREpend  new descr to existing one.
  value <- as.character(value)
  attr(x, "descr") <- if (length(value)==0) NULL else
  if(value[1]=="^") value[-1] else c(value, attr(x, "descr"))
  x
}
## ---
tit <- function(x) attr(x,"tit")
## ---
"tit<-" <- function(x, value) ## ! argument must be 'value'. demanded by attr
{
  attr(x, "tit") <- value
  x
}
## ---
stamp <- function(sure=TRUE, outer.margin = NULL,
                  project=getOption("project"), step=getOption("step"),
                  stamp=getOption("stamp"), ...)
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
  if (length(stamp)==0) {
    message(":stamp: setting options(stamp=1)")
    options(stamp=1)
    stamp <- 1
  }
  if (length(outer.margin)==0) outer.margin <- par('oma')[4]>0
  t.txt <- date()
  t.txt <- paste(substring(t.txt,5,10),",",substring(t.txt,22,23),"/",
                 substring(t.txt,13,16),sep="")
  if (length(project)>0) t.txt <- paste(t.txt,project,sep=" | ")
  if (length(step)>0) t.txt <- paste(t.txt,step,sep=" | ")
  if (sure || stamp==2 || ( stamp==1 && (
                              ##     last figure on page
                              { t.tfg <- par("mfg") ; all(t.tfg[1:2]==t.tfg[3:4]) }
                              || isTRUE(outer.margin) ))  )
    mtext(t.txt, 4, cex = 0.6, adj = 0, outer = outer.margin, ...)
  invisible(t.txt)
}
