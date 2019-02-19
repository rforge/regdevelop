extractNames <- function (x, orig=NULL) {
  if (length(grep("cbind", x))) 
    eval( parse(text=
            paste('c("',
            gsub(" *, *",'","', sub("cbind *\\((.*)\\)", "\\1", x)),
                  '")', sep="") ) )  else  {
                    if(is.null(orig)) NULL else
                    paste("V",seq_along(orig), sep="")
                  }
}

RNAMES <- function (x) if (!is.null(dim(x))) row.names(x) else names(x)
