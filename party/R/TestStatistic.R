
# $Id$

### Wrapper for functions defined in ./src/TestStatistic.c

maxabsTestStatistic <- function(t, mu, Sigma, tol = sqrt(.Machine$double.eps)) {
    storage.mode(t) <- "double"
    storage.mode(mu) <- "double"
    storage.mode(Sigma) <- "double"
    storage.mode(tol) <- "double"
    
    if (length(t) != length(mu) || length(t) != nrow(Sigma)) 
        error("dimensions don't match")
        
    .Call("R_maxabsTestStatistic", t, mu, Sigma, tol, PACKAGE = "partylab2")
}

quadformTestStatistic <- function(t, mu, Sigma, tol = sqrt(.Machine$double.eps)) {
    storage.mode(t) <- "double"
    storage.mode(mu) <- "double"
    storage.mode(Sigma) <- "double"
    storage.mode(tol) <- "double"
    
    if (length(t) != length(mu) || length(t) != nrow(Sigma)) 
        error("dimensions don't match")
        
    SigmaPlus <- MPinv(Sigma, tol = tol)
    .Call("R_quadformTestStatistic", t, mu, SigmaPlus, PACKAGE = "partylab2")
}

