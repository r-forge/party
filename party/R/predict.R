
getNodeNumber = function(x, VarList, obs) {
    node = 1
    while(TRUE) {
        if (class(x@nodes[[node]]) == "TerminalNode") {
            return(x@nodes[[node]]@number)
        }
        pselect = x@nodes[[node]]@primarysplit@variable
        if (class(VarList@inputs[[pselect]]) == "OrderedVariable") {
            left = obs[pselect] <= x@nodes[[node]]@primarysplit@cutpoint
        } else {
            lev = as.integer(obs[pselect])
            left = dual2logical(x@nodes[[node]]@primarysplit@cutpoint)[lev] 
        }
        node = ifelse(left, x@nodeindex[node,2], x@nodeindex[node,3])
    }
}

getWeights = function(x, number) {
    x@nodes[[number]]@weights
}

treepredict = function(x, VarList, newobs) {
    y = VarList@Scores
    f = sapply(newobs, is.factor)
    for (i in 1:ncol(newobs))
      if (f[i]) newobs[,i] = as.numeric(newobs[,i])
    newobs = as.matrix(newobs)

    pred = matrix(0, nr = nrow(newobs), nc = ncol(y))
    
    for (i in 1:nrow(newobs)) {
        node = getNodeNumber(x, VarList, newobs[i,])
        pred[i,] = (x@nodes[[node]]@weights %*% y) / 
                           sum(x@nodes[[node]]@weights) 
    }
    if (class(VarList@response) == "CategoricalVariable")
        colnames(pred) = VarList@response@levels
    rownames(pred) = rownames(newobs)
    pred
}  
