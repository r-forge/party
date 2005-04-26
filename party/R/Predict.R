
# $Id$

oldpredict <- function(object, mincriterion = 0, 
    newdata = NULL, type = c("class", "prob"), ...) {

    type <- match.arg(type)

    ### mimic old-style predict functions
    response <- object@responses
    if (all(response@is_nominal || response@is_ordinal)) {
        if (type == "class") {
            return(object@predict_response(newdata = newdata, 
                   mincriterion = mincriterion, ...))
        } else {
            return(object@cond_distr_response(newdata = newdata,
                   mincriterion = mincriterion, ...))
        }
    }

    if (any(response@is_censored))
        return(object@cond_distr_response(newdata = newdata,                
               mincriterion = mincriterion, ...))

    return(object@predict_response(newdata = newdata,                
           mincriterion = mincriterion, ...)) 
}


predict.BinaryTree <- function(object, newdata = NULL, mincriterion = 0, 
    type = c("class", "prob"), ...) {

    oldpredict(object = object, mincriterion = mincriterion,
               newdata = newdata, type = type, ...)
}
