
splitordered = function(v, pselect, cw, S = NULL) {
    svar = v@inputs[[pselect]]
    x = v@Weights[svar@columns,]
    ox = svar@order

    if (is.null(S)) S = v@Scores
    if (length(svar@whichNA) > 0) cw[svar@whichNA] = 0

    evS = .Call("evS", S, cw) 
    ES = evS[[1]]
    VS = evS[[2]]
    m = evS[[3]]

    mmin = ceiling(m * v@control@minprob) 
    mmax = floor(m * (1 - v@control@minprob)) 

    f1 = m / (m - 1)
    f2 = 1 / (m - 1)

    L = 0
    w = 0
    tmax = 0
    imax = 1

    criterion = rep(0, nrow(S))

    for (i in 1:(nrow(S)-1)) {
        j = ox[i]
        if (i > 1 && x[j] == tx) criterion[ox[i-1]] = 0
        tx = x[j]
        if (cw[j] == 0) next
        w = w + 1
        if (w > mmax) break
        L = L + S[j,] * cw[j]
        E = w * evS[[1]]
        V = f1 * VS * w - f2 * VS * w^2
        if (w >= mmin) 
            criterion[j] = max(abs((L - E)/sqrt(V)))
    }

    cutpoint = x[which.max(criterion)]
    sp = new("OrderedSplit", variable = pselect, 
             cutpoint = cutpoint, totheleft = TRUE)
    sp
}
