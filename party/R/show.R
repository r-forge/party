
pre = function(ni) {
    pres = vector(length=nrow(ni), mode = "list")
    pres[[1]] = NULL
    for (i in 2:nrow(ni)) {
        dp = which(ni[,2] == i | ni[,3] == i)
        pres[[i]] = c(pres[[dp]], dp)
    }
    pres
}

show.BinaryTree = function(object) {
    
    depth = pre(object@nodeindex)
    depth = unlist(lapply(depth, length))
    
    for (n in 1:nrow(object@nodeindex)) {
        if (is.na(object@nodeindex[n, 2])) next
        nd = object@nodes[[n]]
        if (class(nd@primarysplit) == "ContinuousSplit") 
        cat(rep(" ", depth[n]), n, " Var: ",
            object@treegrow@inputs[[nd@primarysplit@variable]]@name, 
            " cutpoint: ", nd@primarysplit@cutpoint, "\n")
        else 
        cat(rep(" ", depth[n]), n, " Var: ",
            object@treegrow@inputs[[nd@primarysplit@variable]]@name, 
            " levels: ", nd@primarysplit@levelset, "\n")
    }
}

setMethod("show", "BinaryTree", function(object) show.BinaryTree(object))
