treedesign = function(MFF) {
    responseterms = terms(MFF@formulae@response)
    attr(responseterms, "intercept") = 0
    MFF@responseMatrix = model.matrix(responseterms, MFF@response)
    designterms = terms(MFF@formulae@input)
    attr(designterms, "intercept") = 0
    fact = unlist(lapply(MFF@input, function(x) 
                     is.factor(x) && !is.ordered(x)))
    if (any(fact)) {
        eval(parse(text=paste("thisctrs <- list(",
                         paste(colnames(MFF@input)[fact], 
                        "=\"ct\"", collapse=", "), ")")))
        MFF@designMatrix = model.matrix(designterms, MFF@input, 
                                    contrasts.arg = thisctrs)
    } else {
        MFF@designMatrix = model.matrix(designterms, MFF@input)
    }

    varindex = attr(MFF@designMatrix, "assign")

    inputs = vector(length = ncol(MFF@input), mode = "list")
    for (i in 1:ncol(MFF@input)) {
        if (fact[i]) {
            var = new("CategoricalVariable")
            var@name = colnames(MFF@input)[i]
            var@columns = which(varindex == i)
            var@whichNA = which(is.na(MFF@input[[i]]))
            var@levels = levels(MFF@input[[i]])
        } else {
            var = new("OrderedVariable")
            var@name = colnames(MFF@input)[i]
            var@columns = which(varindex == i)
            var@whichNA = which(is.na(MFF@input[[i]]))
            var@order = order(MFF@input[[i]])
        }
        inputs[[i]] = var
    }

    if (is.factor(MFF@response[[1]])) {
         response = new("CategoricalVariable")
         response@name = colnames(MFF@response)[1]
         response@columns = 1:nlevels(MFF@response[[1]])
         response@whichNA = NULL
         response@levels = levels(MFF@response[[1]])
     } else {
         response = new("OrderedVariable")
         response@name = colnames(MFF@response)[1]
         response@columns = as.integer(1)
         response@whichNA = NULL
         response@order = order(MFF@response[[1]])
     }
                      
     S = MFF@responseMatrix 
     W = t(MFF@designMatrix) 
     W[is.na(W)] = 0

     VarList = new("TreeGrow", p = length(inputs), inputs = inputs, 
                    response = response, Scores = S, Weights = W, control =
                    new("GrowControl", minsplit = 10, minstat = 1.96))
     return(VarList)
}

