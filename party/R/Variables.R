
# $Id$

initVariableFrame.df <- function(obj, fun = NULL) {

    RET <- new("VariableFrame", nrow(obj), ncol(obj))
    
    is_ordinal <- sapply(obj, is.ordered)
    is_nominal <- sapply(obj, is.factor) & !is_ordinal

    jt <- c()

    ### speedup: extract the slots and re-assign them afterwards
    scores <- RET@scores
    levels <- RET@levels
    variables <- obj
    transformations <- RET@transformations
    ordering <- RET@ordering
    has_missings <- RET@has_missings
    whichNA <- RET@whichNA
 
    for (j in 1:ncol(obj)) {
        x <- variables[[j]]
        if (is_ordinal[j]) {
            sc <- attr(x, "scores")
            if (is.null(sc))
                sc <- 1:nlevels(x)
            storage.mode(sc) <- "double"
            scores[[j]] <- sc
        }

        if (is.factor(x)) {
            xt <- sapply(levels(x), function(l) as.numeric(x == l))
            storage.mode(xt) <- "double"
            levels[[j]] <- attr(x, "levels")

            # <FIXME> storage mode of nominal and ordered factors are different!!!       
            if (!is_ordinal[j]) {
                storage.mode(RET@variables[[j]]) <- "integer"
                tordering <- NULL
            } else {
                storage.mode(variables[[j]]) <- "double"
                tordering <- order(as.numeric(x))
            }
            # </FIXME>
        } else {
            if (!is.null(fun)) 
                xt <- matrix(fun(x), ncol = 1)
            else 
                xt <- matrix(x, ncol = 1)
            storage.mode(xt) <- "double"
            tordering <- order(x)
            storage.mode(tordering) <- "integer"
            storage.mode(variables[[j]]) <- "double"
        }
        nas <- is.na(x)
        xt[nas,] <- 0
        transformations[[j]] <- xt
        ordering[[j]] <- tordering
        has_missings[j] <- any(nas)
        whichNA[[j]] <- which(nas)
        ### this is suboptimal
        jt <- cbind(jt, xt)
    }            
    RET@jointtransf <- jt
    RET@is_nominal <- is_nominal
    RET@is_ordinal <- is_ordinal
    RET@is_censored <- rep(FALSE, ncol(obj))
    RET@variables <- variables
    RET@scores <- scores
    RET@levels <- levels
    RET@transformations <- transformations
    RET@ordering <- ordering
    RET@has_missings <- has_missings
    RET@whichNA <- whichNA
    RET
}

setGeneric(name = "initVariableFrame",
           def = function(obj, fun = NULL)
               standardGeneric("initVariableFrame")
)

setMethod("initVariableFrame", signature = c("data.frame", "ANY"), 
    definition = initVariableFrame.df
)

initVariableFrame.Surv <- function(obj, fun = logrank_scores) {

    RET <- new("VariableFrame", nrow(obj), as.integer(1))
    RET@variables <- list(obj)

    RET@is_nominal <- FALSE
    RET@is_ordinal <- FALSE
    
    RET@transformations <- list(matrix(fun(obj), ncol = 1))

    RET@ordering <- list(order(RET@transformations[[1]]))
    RET@has_missings <- any(is.na(obj))
    RET@whichNA <- list(which(is.na(obj)))
    RET@jointtransf <- RET@transformations[[1]]
    RET@is_censored <- TRUE
    RET
}

setMethod("initVariableFrame", signature = c("Surv", "ANY"),
    definition = initVariableFrame.Surv
)

initVariableFrame.matrix <- function(obj, fun = NULL) {

    RET <- new("VariableFrame", nrow(obj), ncol(obj))
    objDF <- as.data.frame(obj)
    class(objDF) <- "list"
    RET@variables <- objDF

    RET@is_nominal <- rep(FALSE, ncol(obj))
    RET@is_ordinal <- rep(FALSE, ncol(obj))
   
    if (!is.null(fun)) 
        RET@transformations <- lapply(objDF, function(x) 
                                      matrix(fun(x), ncol = 1))
    else
        RET@transformations <- lapply(objDF, function(x) matrix(x, ncol = 1))

    RET@ordering <- lapply(RET@transformations, order)
    RET@has_missings <- unlist(lapply(RET@transformations, function(x) 
                               any(is.na(x))))
    if (any(is.na(obj)))
        RET@whichNA <- lapply(objDF, function(x) which(is.na(x)))
    RET@jointtransf <- matrix(unlist(RET@transformations), 
                              nrow = nrow(obj), byrow = TRUE)
    RET@is_censored <- rep(FALSE, ncol(obj))
    RET
}

setMethod("initVariableFrame", signature = c("matrix", "ANY"),
    definition = initVariableFrame.matrix
)
