
cexpcov = function(W, S, cw, cov = FALSE) {
    if (cov)
        RET = .Call("ec", W, S, cw)
    else 
        RET = .Call("ev", W, S, cw)
    names(RET) = c("ET", "VT")
    RET
}


ct <- function(x, contrasts=FALSE) {
  a <- diag(length(x))
  colnames(a) <- x
  a
}  
   
   
logical2dual = function(x) {
    if (x[1]) 
        as.numeric(paste(as.numeric(x), collapse = ""))
    else 
        -as.numeric(paste(as.numeric(!x), collapse = ""))
}


dual2logical = function(x) {
    bin = c()
    i = 10   
    xsign = sign(x)
    x = x*xsign
    while(TRUE) {
        bin = c(bin, (x / i - floor(x / i)) > 0)
        x = floor(x / i)
        if (x == 0) break()
    } 
    if (xsign < 0)
        return(!rev(bin))
    else 
        return(rev(bin))
}
 
