
# $Id$

predict.BinaryTree <- function(object, ...) {
    conditionalTree@predict(object, ...)
}

setGeneric("weights", function(object, ...) standardGeneric("weights"))

setMethod("weights", signature = signature(object = "BinaryTree"),
    definition = function(object, newdata = NULL, ...)
        object@prediction_weights(newdata = newdata, ...)
)


setGeneric("where", function(object, ...) standardGeneric("where"))

setMethod("where", signature = signature(object = "BinaryTree"),
    definition = function(object, newdata = NULL, ...)
        object@get_where(newdata = newdata, ...)
)
