
foo = function(object, h = 0) {
    if (length(object@nodes) == 1) {
        z = list(object@nodes[[1]]@number)
        attr(z, "label") = paste("Node", object@nodes[[1]]@number)
        attr(z, "members") = 1
        attr(z, "height") = h
        attr(z, "leaf") = TRUE
    } else {
        two = 2:2
        nrleft = object@nodeindex[1,2]
        nrright = object@nodeindex[1,3]
        ileft = (object@nodeindex[,1] >= nrleft & 
                 object@nodeindex[,1] < nrright)
        iright = (object@nodeindex[,1] >= nrright)
        z = vector(length = 2, mode = "list")
        split = object@nodes[[1]]@primarysplit
        nv = object@treegrow@inputs[[split@variable]]@name
        if (class(split) == "ContinuousSplit") {
            attr(z, "edgetext") = as.expression(bquote(.(nv) <= .(split@cutpoint)))
        } else {
            lev = object@treegrow@inputs[[split@variable]]@levels
            attr(z, "edgetext") = paste(nv, "=",
                paste("\{", paste(lev[split@levelset], collapse = ","), "\}", 
                      collapse = ""))
       }
        attr(z, "label") = paste("Node", object@nodes[[1]]@number)

        leftt = object
        leftt@nodeindex = leftt@nodeindex[ileft,,drop=FALSE]
        leftt@nodes = leftt@nodes[ileft]
        rightt = object
        rightt@nodeindex = rightt@nodeindex[iright,,drop=FALSE]
        rightt@nodes = rightt@nodes[iright]

        z[[1]] = foo(leftt, h - 1)
        z[[2]] = foo(rightt, h - 1)
        attr(z, "members") = attr(z[[1]], "members") + attr(z[[2]], "members")
        attr(z, "height") = h
        if (attr(z[[1]], "members") == 1 || attr(z[[2]], "members") == 1) 
            attr(z, "midpoint") = 0.5 
        else 
            attr(z, "midpoint") = (attr(z[[1]], "members") + attr(z[[1]],
                                   "midpoint") + attr(z[[2]], "midpoint"))/2
    }
    class(z) = "dendrogram"
    z
}

as.dendrogram.BinaryTree = function(object) {
    maxheight = max(unlist(sapply(pre(object@nodeindex), length)))
    foo(object, maxheight)
}

setMethod("as.dendrogram", "BinaryTree", function(object, ...) {
    as.dendrogram.BinaryTree(object)
})

    
