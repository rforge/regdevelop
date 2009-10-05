naresid.exclude <- 
function (omit, x, ...) 
{
    if (length(omit) == 0 || !is.numeric(omit)) 
        stop("invalid argument 'omit'")
    if (length(x) == 0) 
        return(x)
    if (is.matrix(x)) {
        n <- nrow(x)
        keep <- rep.int(NA, n + length(omit))
        keep[-omit] <- 1:n
        x <- x[keep, , drop = FALSE]
        temp <- rownames(x)
        if (length(temp)) {
            temp[omit] <- names(omit)
            rownames(x) <- temp
        }
    }
    else {
        n <- length(x)
        keep <- rep.int(NA, n + length(omit))
        keep[-omit] <- 1:n
        x <- x[keep]
        temp <- names(x)
        if (length(temp)) {
            temp[omit] <- names(omit)
            names(x) <- temp
        }
    }
    BR()
    x
}
