
u.true <- function (x) length(x)>0 && is.logical(x) && (!is.na(x)) && all(x)
u.notfalse <-
  function (x) !(length(x)==1 && is.logical(x) && (!is.na(x)) && !x)
u.isnull <- function(x)  length(x)==0||all(is.na(x))
"%nin%" <- function (x,y) !x%in%y
i.last <- function(data, n=1) data[sign(n)*(((ldt <- length(data))-abs(n)+1):ldt)]
u.debug <- function () u.true(rloptions("debug"))
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
## ---------------------------------------
DB <- function (on=TRUE) options(error=if(on) recover else NULL, warn=on)
