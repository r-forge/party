
# $Id$

### Wrapper for functions defined in ./src/LinearStatistic.c

LinearStatistic <- function(x, y, weights) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_LinearStatistic", x, y, weights, PACKAGE = "partylab2")
}

ExpectCovarInfluence <- function(y, weights) {
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_ExpectCovarInfluence", y, weights, PACKAGE = "partylab2")
}

ExpectCovarLinearStatistic <- function(x, y, weights) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    expcovinf <- ExpectCovarInfluence(y, weights)
    .Call("R_ExpectCovarLinearStatistic", x, y, weights, expcovinf,
          PACKAGE = "partylab2")
}

PermutedLinearStatistic <- function(x, y, indx, perm) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"

    if (any(indx < 1 || indx > nrow(y)))
        stop("wrong indices")

    if (any(perm < 1 || perm > nrow(y)))
        stop("wrong indices")

    # C indexing 
    indx <- indx - 1
    perm <- perm - 1
    storage.mode(indx) <- "integer"
    storage.mode(perm) <- "integer"
    .Call("R_PermutedLinearStatistic", x, y, indx, perm, 
          PACKAGE = "partylab2")
}
