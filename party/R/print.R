
# $Id$

prettysplit <- function(x, inames = NULL, ilevels = NULL) {
    names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics")
    if (x$ordered) {
        class(x) <- "orderedSplit"
    } else {
        class(x) <- "nominalSplit"
    }
    if (!is.null(ilevels)) {
        if (!is.null(ilevels[x[["variableID"]]]))
            attr(x$splitpoint, "levels") <- ilevels[[x[["variableID"]]]]
    }
    if (!is.null(inames)) x$variableName <- inames[x[["variableID"]]]
    return(x)
}

prettytree <- function(x, inames = NULL, ilevels = NULL) {
    names(x) <- c("nodeID", "weights", "criterion", "terminal",
                  "psplit", "ssplits", "prediction", "left", "right")
    if (is.null(inames) && extends(class(x), "BinaryTree"))
        inames <- x@inputnames
    if (x$terminal) {
        class(x) <- "TerminalNode"
        return(x)
    }

    x$psplit <- prettysplit(x$psplit, inames = inames, ilevels = ilevels)

    class(x) <- "SplittingNode"
    x$left <- prettytree(x$left, inames = inames, ilevels = ilevels)   
    x$right <- prettytree(x$right, inames = inames, ilevels = ilevels)    
    return(x)
}
 
print.TerminalNode <- function(x, n = 1, ...) {
    cat(paste(paste(rep(" ", n), collapse = ""), x$nodeID, ")* ", sep=""),
        "weights = ", sum(x$weights), ", prediction = ", round(x$prediction, 3), "\n")
}
 
print.SplittingNode <- function(x, n = 1, ...) {
    cat(paste(paste(rep(" ", n), collapse = ""), x$nodeID, ") ", sep=""))
    print(x$psplit, left = TRUE)
    cat("; max(criterion) = ", round(x$criterion[[3]], 3), ", max(test statistic) = ", 
        round(max(x$criterion[[1]]), 3), "\n")
    print(x$left, n + 2)
    cat(paste(paste(rep(" ", n), collapse = ""), x$nodeID, ") ", sep=""))
    print(x$psplit, left = FALSE)
    cat("\n")
    print(x$right, n + 2)
}

print.orderedSplit <- function(x, left = TRUE, ...) {
    if (!is.null(attr(x$splitpoint, "levels"))) {
        sp <- attr(x$splitpoint, "levels")[x$splitpoint]
    } else {
        sp <- x$splitpoint
    }
    if (left) {
        cat(x$variableName, " <= ", sp)
    } else {
        cat(x$variableName, " > ", sp)
    }
}

print.nominalSplit <- function(x, left = TRUE, ...) {
    if (left) {
        txt <- paste("\{", paste(attr(x$splitpoint, "levels")[as.logical(x$splitpoint)], 
                                 collapse = ", "), "\}")
    } else {
        txt <- paste("\{", paste(attr(x$splitpoint, "levels")[!as.logical(x$splitpoint)], 
                                 collapse = ", "), "\}")
    }
    cat(x$variableName, " == ", txt)
}


print.BinaryTreePartition <- function(x, ...)
    print(x@tree)

print.BinaryTree <- function(x, ...)
    print(x@tree)
