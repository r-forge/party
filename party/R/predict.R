
getNodeNumber = function(x, obs) {
    node = getNode(x, 1)
    while(TRUE) {
        if (class(node) == "TerminalNode") {
            return(node@number)
        }
        pselect = node@primarysplit@variable
        split = node@primarysplit
        if (class(split) == "ContinuousSplit") {
            left = obs[pselect] <= split@cutpoint
        } else {
            lev = as.integer(obs[pselect])
            left = any(split@levelset == lev)
        }
        if (left) node = getLeftNode(x, node@number) else
                  node = getRightNode(x, node@number)
    }
}

getNode = function(x, number)
    x@nodes[[which(x@nodeindex[,1] == number)]]

getWeights = function(x, number)
    getNode(x, number)@weights

getLeftNode = function(x, number) {
    indx = which(x@nodeindex[,1] == number)
    indx = which(x@nodeindex[,1] == x@nodeindex[indx,2])
    x@nodes[[indx]]
}

getRightNode = function(x, number) {
    indx = which(x@nodeindex[,1] == number)
    indx = which(x@nodeindex[,1] == x@nodeindex[indx,3])
    x@nodes[[indx]]
}

treepredict = function(x, VarList, newobs) {
    y = VarList@response@values
    f = sapply(newobs, is.factor)
    for (i in 1:ncol(newobs))
      if (f[i]) newobs[,i] = as.numeric(newobs[,i])
    newobs = as.matrix(newobs)

    pred = matrix(0, nr = nrow(newobs), nc = ncol(y))
    
    for (i in 1:nrow(newobs)) {
        nodenr = getNodeNumber(x, newobs[i,])
        pred[i,] = (getWeights(x, nodenr) %*% y) / 
                           sum(getWeights(x, nodenr))
    }
    if (class(VarList@response) == "CategoricalVariable")
        colnames(pred) = VarList@response@levels
    rownames(pred) = rownames(newobs)
    pred
}  
