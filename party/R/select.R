
best = function(v, cw) {
    criterion = rep.int(0, v@p)
    S = v@response@values
    sevS = .Call("evS", S, cw)
    Scw = S * cw
    for (p in 1:v@p) {

        svar = v@inputs[[p]]
        tcw = cw
        if (length(svar@whichNA) > 0) {
            tcw[svar@whichNA] = 0
            es = .Call("evS", S, tcw)
            Sw = S * tcw
        } else {
            es = sevS
            Sw = Scw
        }
        L = as.vector(svar@values %*% Sw)
        if (class(svar) == "OrderedCategoricalVariable") {
            ap = .Call("ec", svar@values, S, cw)
            scores = matrix(rep(svar@scores, 
                     length(L)/length(svar@scores)), nr = 1)
            L = L %*% scores
            EL = ap[[1]] %*% scores
            VL = scores %*% ap[[2]] %*% t(scores)
        } else { 
            ap = .Call("evL", svar@values, S, cw, es)
            EL = ap[[1]]
            VL = ap[[2]]
        } 
        zeros = VL < v@control@varnull
        L[zeros] = 0
        T = L
        T[!zeros] = ((L[!zeros] - EL[!zeros])/sqrt(VL[!zeros]))
        criterion[p] = max(abs(T))
    }
    if (max(criterion) < v@control@minstat) return(NULL)
    criterion
}
