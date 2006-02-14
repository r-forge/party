
# $Id$

.onLoad <- function(lib, pkg) {
    if (!require("methods")) stop("cannot load methods")
    GCtorture <<- FALSE
    .Call("party_init", PACKAGE = "party")
    return(TRUE)
}
