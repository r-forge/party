
treedesign = function(MFF, rtrans = NULL) {

    inputs = vector(length = ncol(MFF@input), mode = "list")
    for (i in 1:ncol(MFF@input)) {
        if (is.factor(MFF@input[[i]])) {
            if (is.ordered(MFF@input[[i]])) {
                var = new("OrderedCategoricalVariable")
                if (is.null(scores)) 
                    var@scores = 1:nlevels(MFF@input[[i]])
                else
                    var@scores = scores
                var@order = order(MFF@input[[i]])
            } else {
                var = new("CategoricalVariable")
            }
            var@values = .Call("dummyvar", MFF@input[[i]], 
                                as.logical(TRUE)) 
            var@coding = as.integer(MFF@input[[i]])
            var@name = colnames(MFF@input)[i]
            var@whichNA = which(is.na(MFF@input[[i]]))
            var@levels = levels(MFF@input[[i]])
            scores = attr(MFF@input[[i]], "scores")
        } else {
            var = new("ContinuousVariable")
            var@name = colnames(MFF@input)[i]
            var@values = matrix(as.numeric(MFF@input[[i]]), nrow = 1)
            var@whichNA = which(is.na(MFF@input[[i]]))
            var@order = order(MFF@input[[i]])
        }
        var@values[,var@whichNA] = 0
        inputs[[i]] = var
    }

    for (i in 1:ncol(MFF@response)) {
        if (i == 2) stop("Multivariate Responses not yet implemented!")
        if (is.factor(MFF@response[[i]])) {
            if (is.ordered(MFF@response[[i]])) {
                var = new("OrderedCategoricalVariable")
                if (is.null(scores)) 
                    var@scores = 1:nlevels(MFF@response[[i]])
                else
                    var@scores = scores
                var@order = order(MFF@response[[i]])
            } else {
                var = new("CategoricalVariable")
            }
            var@values = .Call("dummyvar", MFF@response[[i]], 
                                as.logical(FALSE)) 
            var@coding = as.integer(MFF@response[[i]])
            var@name = colnames(MFF@response)[i]
            var@whichNA = NULL
            var@levels = levels(MFF@response[[i]])
        } else {
            var = new("ContinuousVariable")
            var@name = colnames(MFF@response)[i]
            var@values = matrix(as.numeric(MFF@response[[i]]), ncol = 1)
            var@whichNA = NULL
            var@order = order(MFF@response[[i]])
        }
        response = var
    }
              
    workingresponse = var
    if (!is.null(rtrans)) {
        workingresponse = rtrans(var)
    }
        
    varlist = new("PartyStarters", p = length(inputs), nobs = nrow(response@values), 
                  inputs = inputs, 
                  workingresponse = workingresponse, 
                  response = response,
                  control = new("PartyControl"))
    return(varlist)
}

