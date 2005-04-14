
# $Id$

### Wrapper for functions defined in ./src/Utilsc

qsvd <- function(x) {
    if (!is.matrix(x) || ncol(x) != nrow(x))
        stop(sQuote("x"), " is not a quadratic matrix")

    svdmem <- new("svd_mem", ncol(x)) 
    dummy <- .Call("CR_svd", x, svdmem, PACKAGE = "partylab2")
    return(svdmem@svd)
}

MPinv <- function(x, tol = sqrt(.Machine$double.eps)) {
    if (!is.matrix(x) || ncol(x) != nrow(x))
        stop(sQuote("x"), " is not a quadratic matrix")

    svdmem <- new("svd_mem", ncol(x))
    RET <- .Call("R_MPinv", x, tol, svdmem, PACKAGE = "partylab2")
    return(RET@MPinv)
}

### Median Survival Time, see survival:::print.survfit

mst <- function(x) {
    minmin <- function(y, xx) {

        ### don't complain about Inf 
        ww <- getOption("warn")
        on.exit(options(warn = ww))
        options(warn = -1)

        if (any(!is.na(y) & y==.5)) {
            if (any(!is.na(y) & y <.5))
                .5*(min(xx[!is.na(y) & y==.5]) + min(xx[!is.na(y) & y<.5]))
            else
                .5*(min(xx[!is.na(y) & y==.5]) + max(xx[!is.na(y) & y==.5]))
        } else   min(xx[!is.na(y) & y<=.5])
    }
    med <- minmin(x$surv, x$time)
    return(med)
}

### Logrank scores, see exactRankTests:::cscores.Surv

logrank_scores <- function(x) {

  if (class(x) != "Surv")
      stop(sQuote("x"), " is not of class ", sQuote("Surv"))

  time <- x[,1]
  event <- x[,2]
  N <- length(time)
  ot <- order(time)
  rt <- rank(time, ties.method = "max")
  fact <- event/(N - rt + 1)
  RET <- event - cumsum(fact[ot])[rt]
  attr(RET, "scores") <- "LogRank"
  return(RET)
}
