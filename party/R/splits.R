
splitcontinuous = function(v, pselect, cw) {
    svar = v@inputs[[pselect]]
    x = as.vector(svar@values)
    ox = svar@order
    S = v@response@values
    if (length(svar@whichNA) > 0) cw[svar@whichNA] = 0
    split = spoWH(x, ox, S, cw, v)
    if (max(split$crit) == 0) return(NULL)
    sp = new("ContinuousSplit", 
             variable = pselect,
             cutpoint = split$cut,
             criterium = split$crit,
             totheleft = TRUE)
    sp
}

splitordered = function(v, pselect, cw) {
    svar = v@inputs[[pselect]]
    x = as.vector(svar@coding)
    ox = svar@order
    S = v@response@values
    if (length(svar@whichNA) > 0) cw[svar@whichNA] = 0
    split = spoWH(x, ox, S, cw, v)
    if (max(split$crit) == 0) return(NULL)
    sp = new("OrderedSplit", 
             variable = pselect,
             levelset = 1:split$cut,
             criterium = split$crit,
             totheleft = TRUE)
    sp
}


spoWH = function(x, ox, S, cw, v) {
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
        w = w + cw[j]
        if (w > mmax) break
        L = L + S[j,] * cw[j]
        E = w * ES
        V = f1 * VS * w - f2 * VS * w^2
        if (V < v@control@varnull) next
        if (w >= mmin) 
            ### ordered response factors???
            criterion[j] = max(abs((L - E)/sqrt(V)))
    }

    cutpoint = x[which.max(criterion)]
    list(cut = cutpoint, crit = criterion)
}

splitcategorical = function(v, pselect, cw) {
    svar = v@inputs[[pselect]]
    x = svar@values 
    xc = svar@coding

    S = v@workingresponse@values
    if (length(svar@whichNA) > 0) cw[svar@whichNA] = 0

    splits = vector(length = ncol(S), mode = "list")
    stats = rep(0, ncol(S))

    evS = .Call("evS", S, cw) 
    ES = evS[[1]]
    VS = evS[[2]]
    m = evS[[3]]

    mmin = ceiling(m * v@control@minprob) 
    mmax = floor(m * (1 - v@control@minprob)) 

    f1 = m / (m - 1)
    f2 = 1 / (m - 1)

    T = matrix(standstat(x, S, cw), ncol = ncol(S))

    for (i in 1:ncol(T)) {
        ox = order(T[,i])
        L = 0
        E = 0
        V = 0
        w = 0
        criterion = rep(0, nrow(T))
        for (k in 1:(nrow(T)-1)) {
            tlev = ox[k]
            for (j in 1:v@nobs) {
                if (cw[j] == 0) next
                w = w + (xc[j] == tlev)*cw[j]
                L = L + (xc[j] == tlev)*S[j,] * cw[j]
            }
            if (w > mmax) break
            E = w * ES
            V = f1 * VS * w - f2 * VS * w^2
            if (V < v@control@varnull) next
            if (w >= mmin) 
                ### ordered response factors???
                criterion[k] = max(abs((L - E)/sqrt(V)))
        }
        splits[[i]] = ox[-(1:which.max(abs(criterion)))]
        stats[i] = max(abs(criterion))
    }
    
    levelset = splits[[which.max(stats)]]
    sp = new("CategoricalSplit", variable = pselect, 
             levelset = levelset, totheleft = TRUE, 
             criterium = max(stats))
    sp
}
