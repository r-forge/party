
# $Id$

if (FALSE) {
.First.lib <- function(libname, pkgname) {
    if (!require("methods")) stop("cannot load methods")
    if (!require("grid")) stop("cannot load grid")
    if (!require("survival")) stop("cannot load survival")
    library.dynam("party", pkgname, libname)
    GCtorture <<- FALSE
    .Call("party_init", PACKAGE = "party")
}
}

if (TRUE) {
.onLoad <- function(lib, pkg) {
    if (!require("methods")) stop("cannot load methods")
    if (!require("grid")) stop("cannot load grid")
    if (!require("survival")) stop("cannot load survival")
    if (!require("modeltools")) stop("cannot load modeltools")
    GCtorture <<- FALSE
    .Call("party_init", PACKAGE = "party")
    return(TRUE)
}
}
