treedesign = function(MFF) {
    fact = unlist(lapply(MFF@input, function(x) 
                     is.factor(x) && !is.ordered(x)))

    inputs = vector(length = ncol(MFF@input), mode = "list")
    for (i in 1:ncol(MFF@input)) {
        if (fact[i]) {
            var = new("CategoricalVariable")
            var@coding = as.integer(MFF@input[[i]])
            var@name = colnames(MFF@input)[i]
            var@values = .Call("dummyvar", MFF@input[[i]], 
                               as.logical(TRUE)) 
            var@whichNA = which(is.na(MFF@input[[i]]))
            var@levels = levels(MFF@input[[i]])
        } else {
            var = new("OrderedVariable")
            var@name = colnames(MFF@input)[i]
            var@values = matrix(as.numeric(MFF@input[[i]]), nrow = 1)
            var@whichNA = which(is.na(MFF@input[[i]]))
            var@order = order(MFF@input[[i]])
        }
        var@values[,var@whichNA] = 0
        inputs[[i]] = var
    }

    if (is.factor(MFF@response[[1]])) {
         response = new("CategoricalVariable")
         response@name = colnames(MFF@response)[1]
         response@coding = as.integer(MFF@response[[1]])
         response@values = .Call("dummyvar", MFF@response[[1]],
                               as.logical(FALSE))
         response@whichNA = NULL
         response@levels = levels(MFF@response[[1]])
     } else {
         response = new("OrderedVariable")
         response@name = colnames(MFF@response)[1]
         response@values = matrix(MFF@response[[1]], ncol = 1)
         response@whichNA = NULL
         response@order = order(MFF@response[[1]])
     }
                      
     VarList = new("TreeGrow", p = length(inputs), inputs = inputs, 
                    response = response, control =
                    new("GrowControl", minsplit = 10, minstat = 1.96))
     return(VarList)
}

