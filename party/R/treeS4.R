    

standstat = function(W, S, cw) {
    L = as.vector(W %*% (S * cw))
    ap = cexpcov(W, S, cw, cov = FALSE)
    zeros = ap$VT < 1e-10
    L[zeros] = 0
    T = L
    T[!zeros] = ((L[!zeros] - ap$ET[!zeros])/sqrt(ap$VT[!zeros]))
    T
}

surrogates = function(v, pselect, leftw, cw, n = 1) {

    criterion = rep.int(0, v@p)
    myS = matrix(leftw, ncol = 1)
    for (p in 1:v@p) {
        if (p == pselect) next
        tcw = cw
        if (length(v@inputs[[p]]@whichNA) > 0)
            tcw[v@inputs[[p]]@whichNA] = 0
        criterion[p] = max(abs(standstat(v@inputs[[p]]@values, myS, tcw)))
    }

    surr = rev(order(criterion))[1:n]
    ss = c()
    for (p in surr) {
        svar = v@inputs[[p]]
        if (class(svar) == "ContinuousVariable" || 
            class(svar) == "OrderedCategoricalVariable") {
            split = splitordered(v, pselect, cw, S = myS)
            leftcw = cw * (as.vector(svar@values) <= split@cutpoint)
        } else {
            csplit = splitcategorical(v, pselect, cw, S = myS)
            leftlevels = csplit@levelset
            Wtmp = svar@values
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

    svar = v@inputs[[pselect]]

    if (class(svar) == "ContinuousVariable")
       split = splitcontinuous(v, pselect, cw)
    if (class(svar) == "OrderedCategoricalVariable")
       split = splitordered(v, pselect, cw) 
    if (class(svar) == "CategoricalVariable")
       split = splitcategorical(v, pselect, cw)

    if (is.null(split)) {
        tnode = new("TerminalNode", number = 0, weights = cw)
        return(tnode)
    }

    node = new("Node", number = 0,
               primarysplit = split, criterion = criterion,
               surrogatesplits = list(), weights = cw)
    return(node)
}

treegrow = function(v, sn, sni, gi, gt, cw = NULL, nr = 1) {

    if (is.null(cw)) cw = rep(1, nrow(v@workingresponse@values))

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
    svar = v@inputs[[pselect]]

    if (class(split) == "ContinuousSplit") {
        leftcw = cw * (as.vector(svar@values) <= split@cutpoint)
        rightcw = cw * (as.vector(svar@values) >  split@cutpoint)
    } else { 
       Wtmp = svar@values 
       leftlevels = (1:nrow(Wtmp)) %in% split@levelset
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
    nobs = v@nobs
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

