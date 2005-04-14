
# $Id$

Split <- function(x, y, weights, splitctrl) {

    if (is.factor(y))
        ym <- sapply(levels(y), function(l) as.numeric(y == l))
     else 
        ym <- matrix(y, ncol = 1)
    storage.mode(ym) <- "double"

    if (is.factor(x)) {
        xm <- sapply(levels(x), function(l) as.numeric(x == l))
        storage.mode(xm) <- "double"
        xc <- as.numeric(x)
        storage.mode(xc) <- "integer"
        lecxy <- new("LinStatExpectCovar", ncol(xm), ncol(ym))
        lec <- new("LinStatExpectCovar", as.integer(1), ncol(ym))
        eci <- ExpectCovarInfluence(ym, weights)
        split <- .Call("R_splitcategorical", xm, xc, ym, weights, lec, lecxy,
                       eci, splitctrl, PACKAGE = "partylab2")
    } else {
        ox <- order(x)
        storage.mode(ox) <- "integer"
        xm <- matrix(x, ncol = 1)
        storage.mode(xm) <- "double"
        lec <- new("LinStatExpectCovar", as.integer(1), ncol(ym))
        eci <- ExpectCovarInfluence(ym, weights)
        split <- .Call("R_split", xm, ym, weights, ox, lec,
                       eci, splitctrl, PACKAGE = "partylab2")
    }
    split
}
