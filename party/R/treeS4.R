    

standstat = function(W, S, cw) {
    L = as.vector(W %*% (S * cw))
    ap = cexpcov(W, S, cw, cov = FALSE)
    zeros = ap$VT < 1e-10
    L[zeros] = 0
    T = L
    T[!zeros] = ((L[!zeros] - ap$ET[!zeros])/sqrt(ap$VT[!zeros]))
    T
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
    
    levelset = splits[[which.max(stats)]]
    #    cutpoint = logical2dual((1:nrow(x)) %in% splits[[which.max(stats)]])
    sp = new("CategoricalSplit", variable = i, levelset = levelset, totheleft = TRUE)
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
            split = splitordered(v, pselect, cw, S = myS)
            leftcw = cw * (v@Weights[ri,] <= split@cutpoint)
        } else {
            ri = varselect@columns
            csplit = splitcategorical(v, pselect, cw, S = myS)
            leftlevels = csplit@levelset
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
    sevS = .Call("evS", v@Scores, cw)
    S = v@Scores
    Scw = S * cw
    for (p in 1:v@p) {
        tcw = cw
        if (length(v@inputs[[p]]@whichNA) > 0) {
            tcw[v@inputs[[p]]@whichNA] = 0
            es = .Call("evS", v@Scores, cw)
            Sw = S * tcw
        } else {
            es = sevS
            Sw = Scw
        }
        ri = v@inputs[[p]]@columns
        W = v@Weights[ri,,drop = FALSE]
        L = as.vector(W %*% Sw)
        ap = .Call("evL", W, S, cw, es)
        zeros = ap[[2]] < 1e-10
        L[zeros] = 0
        T = L
        T[!zeros] = ((L[!zeros] - ap[[1]][!zeros])/sqrt(ap[[2]][!zeros]))
        criterion[p] = max(abs(T))
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
    split = nd@primarysplit
    varselect = v@inputs[[pselect]]
    ri = varselect@columns

    if (class(split) == "OrderedSplit") {
        leftcw = cw * (v@Weights[ri,] <= split@cutpoint)
        rightcw = cw * (v@Weights[ri,] >  split@cutpoint)
    } else {
       leftlevels = (1:length(ri)) %in% split@levelset
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

