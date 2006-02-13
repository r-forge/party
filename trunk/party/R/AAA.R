
# $Id$

.onLoad <- function(lib, pkg) {
    if (!require("methods")) stop("cannot load methods")
    if (!require("grid")) stop("cannot load grid")
    if (!require("survival")) stop("cannot load survival")
    if (!require("modeltools")) stop("cannot load modeltools")
    if (!require("strucchange")) stop("cannot load sandwich")
    # if (!require("coin")) stop("cannot load coin")
    GCtorture <<- FALSE
    .Call("party_init", PACKAGE = "party")
    return(TRUE)
}
