    

standstat = function(W, S, cw) {
    L = as.vector(W %*% (S * cw))
    ap = cexpcov(W, S, cw, cov = FALSE)
    zeros = ap$VT < 1e-10
    L[zeros] = 0
    T = L
    T[!zeros] = ((L[!zeros] - ap$ET[!zeros])/sqrt(ap$VT[!zeros]))
    T
}

splitordered = function(v, i, cw, S = NULL) {
    svar = v@inputs[[i]]
    x = v@Weights[svar@columns,]

    if (is.null(S)) S = v@Scores
    if (length(svar@whichNA) > 0) cw[svar@whichNA] = 0

    ux = sort(unique(x[cw > 0]))
    Wsplit = .Call("contcat", x, ux) 
    if (ncol(S) == 1 || all(rowSums(S) == 1))
        psplit = matrix(standstat(Wsplit, S[,1,drop=FALSE], cw), nc = 1)
    else 
        psplit = matrix(standstat(Wsplit, S, cw), 
                        ncol = ncol(S))
    cutpoint = ux[which.max(apply(abs(psplit), 1, max))]
    sp = new("Split", variable = i, cutpoint = cutpoint, totheleft = TRUE)
    sp
}

splitcategorical = function(v, i, cw, S = NULL) {
    svar = v@inputs[[i]]
    x = v@Weights[svar@columns,]

    if (is.null(S)) S = v@Scores
    if (length(svar@whichNA) > 0) cw[svar@whichNA] = 0

    splits = vector(length = ncol(S), mode = "list")
    stats = rep(0, ncol(S))

    for (j in 1:ncol(S)) {
        T = standstat(x, S[,j,drop = FALSE], cw)
        ox = order(T)
        W = apply(x[ox,], 2, cumsum)
        s = matrix(standstat(W, S, cw), nc = ncol(S))
        s = apply(abs(s), 1, max)
        splits[[j]] = ox[-(1:which.max(abs(s)))]
        stats[j] = max(abs(s))
    }
    
    cutpoint = logical2dual((1:nrow(x)) %in% splits[[which.max(stats)]])
    sp = new("Split", variable = i, cutpoint = cutpoint, totheleft = TRUE)
    sp

}

surrogates = function(v, pselect, leftw, cw, n = 1) {
    criterion = rep.int(0, v@p)
    myS = matrix(leftw, ncol = 1)
    for (p in 1:v@p) {
        if (p == pselect) next
        tcw = cw
        if (length(v@inputs[[p]]@whichNA) > 0)
            tcw[v@inputs[[p]]@whichNA] = 0
        ri = v@inputs[[p]]@columns
        criterion[p] = max(abs(standstat(v@Weights[ri,,drop = FALSE], myS, tcw)))
    }
    surr = rev(order(criterion))[1:n]
    
    ss = c()
    for (p in surr) {
        varselect = v@inputs[[p]]
        if (class(varselect) == "OrderedVariable") {
            ri = varselect@columns
            cutpoint = splitordered(v, pselect, cw, S = myS)
            leftcw = cw * (v@Weights[ri,] <= cutpoint)
        } else {
            ri = varselect@columns
            cutpoint = splitcategorical(v, pselect, cw, S = myS)
            leftlevels = dual2logical(cutpoint)
            Wtmp = v@Weights[ri,]
            leftcw = cw * (colSums(Wtmp[leftlevels,,drop=FALSE]))
        }
        if (sum(leftcw * leftw * cw) < sum((1-leftcw) * leftw * cw))
            direction = 1  ### right
        else
            direction = 0  ### left
        ss = rbind(ss, c(p, cutpoint, direction))
    }
    ss
}
    
best = function(v, cw) {
    criterion = rep.int(0, v@p)
    for (p in 1:v@p) {
        tcw = cw
        if (length(v@inputs[[p]]@whichNA) > 0)
            tcw[v@inputs[[p]]@whichNA] = 0
        ri = v@inputs[[p]]@columns
        criterion[p] = max(abs(standstat(v@Weights[ri,,drop = FALSE], v@Scores, tcw)))
    }
    if (max(criterion) < v@control@minstat) return(NULL)
    criterion
}

node = function(v, cw) {

    if (sum(cw) < v@control@minsplit) {
        tnode = new("TerminalNode", number = 0, weights = cw)
        return(tnode) 
    }

    criterion = best(v, cw)
    if (is.null(criterion)) {
        tnode = new("TerminalNode", number = 0, weights = cw)
        return(tnode) 
    }
    pselect = which.max(criterion)

    varselect = v@inputs[[pselect]]

    if (class(varselect) == "OrderedVariable") {
       ri = varselect@columns
       split = splitordered(v, pselect, cw)
    } else {
       ri = varselect@columns
       split = splitcategorical(v, pselect, cw)
    }

    node = new("Node", number = 0,
               primarysplit = split, criterion = criterion,
               surrogatesplits = list(), weights = cw)
    return(node)
}

treegrow = function(v, sn, sni, gi, gt, cw = NULL, nr = 1) {

    if (is.null(cw)) cw = rep(1, ncol(v@Weights))

    nd = node(v, cw)
    nd@number = nr
    sn(nd)
    if (class(nd) == "TerminalNode") {
        sni(c(nr, NA, NA))
        return(NULL)
    }

    nd@number = nr
    pselect = nd@primarysplit@variable
    cutpoint = nd@primarysplit@cutpoint
    varselect = v@inputs[[pselect]]
    ri = varselect@columns
#    cat("Split in ", varselect@name, " <= ", 
#        cutpoint, " with ", sum(cw), "obs \n")

    if (class(varselect) == "OrderedVariable") {
        leftcw = cw * (v@Weights[ri,] <= cutpoint)
        rightcw = cw * (v@Weights[ri,] >  cutpoint)
    } else {
       leftlevels = dual2logical(cutpoint)
       Wtmp = v@Weights[ri,]
       leftcw = cw * (colSums(Wtmp[leftlevels,,drop=FALSE]))
       rightcw = cw * (colSums(Wtmp[!leftlevels,,drop=FALSE]))
    }

    if (any(is.na(leftcw))) {
        prob = mean(leftcw * cw, na.rm = TRUE)
        leftcw[is.na(leftcw)] = rbinom(sum(is.na(leftcw)), 
            size = 1, prob = prob)
        rightcw[is.na(rightcw)] = 1 - leftcw[is.na(rightcw)]
    }

    treegrow(v, sn, sni, gi, gt, leftcw, nr = nr + 1)
    nextnr = gi() + 1    
    treegrow(v, sn, sni, gi, gt, rightcw, nr = nextnr)
    sni(c(nr, nr+1, nextnr))
}

stree = function(v) {
    nobs = ncol(v@Weights)
    tree = new("BinaryTree", 
               nodes = vector(length = floor(nobs / 2), mode = "list"), 
               nodeindex = matrix(0, nrow = floor(nobs)/2, ncol = 3),
               treegrow = v)
    sn = function(node) tree@nodes[[node@number]] <<- node
    sni = function(indx) tree@nodeindex[indx[1],] <<- indx
    gi = function() max(tree@nodeindex[,1])
    gt = function() tree
    treegrow(v, sn, sni, gi, gt)
    indx = unlist(lapply(tree@nodes, is.null))
    tree@nodes = tree@nodes[!indx]
    tree@nodeindex = tree@nodeindex[!indx,]
    tree
}

treepredict = function(x, VarList, newobs) {
    y = VarList@Scores
    node = 1
    f = sapply(newobs, is.factor)
    for (i in 1:ncol(newobs))
      if (f[i]) newobs[,i] = as.numeric(newobs[,i])

    pred = matrix(0, nr = nrow(newobs), nc = ncol(y))
    for (i in 1:nrow(newobs)) {
        node = 1
        while(TRUE) {
            if (class(x@nodes[[node]]) == "TerminalNode") {
                pred[i,] = (x@nodes[[node]]@weights %*% y) / sum(x@nodes[[node]]@weights) 
                break
            }
            pselect = x@nodes[[node]]@primarysplit@variable
            if (class(VarList@inputs[[pselect]]) == "OrderedVariable") {
                left = newobs[i,pselect] <= x@nodes[[node]]@primarysplit@cutpoint
            } else {
                lev = as.integer(newobs[i,pselect])
                left = dual2logical(x@nodes[[node]]@primarysplit@cutpoint)[lev] 
            }
            node = ifelse(left, x@nodeindex[node,2], x@nodeindex[node,3])
        }
    }
    if (class(VarList@response) == "CategoricalVariable")
    colnames(pred) = VarList@response@levels
    rownames(pred) = rownames(newobs)
    pred
}  

