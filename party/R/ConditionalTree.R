
# $Id$

ctreefit <- function(object, controls, weights = NULL, fitmem = NULL, ...) {

    if (!extends(class(object), "LearningSample"))
        stop(sQuote("object"), " is not of class ", sQuote("LearningSample"))
    if (!extends(class(controls), "TreeControl"))
        stop(sQuote("controls"), " is not of class ", sQuote("TreeControl"))

    if (is.null(fitmem)) 
        fitmem <- TreeFitMemory(object, TRUE)
    if (!extends(class(fitmem), "TreeFitMemory"))
        stop(sQuote("fitmem"), " is not of class ", sQuote("TreeFitMemory"))

    if (is.null(weights))
        weights <- object@weights
    if (length(weights) != object@nobs || storage.mode(weights) != "double")
        stop(sQuote("weights"), " are not a double vector of ", 
             object@nobs, " elements")
    if (max(abs(floor(weights) -  weights)) > sqrt(.Machine$double.eps))
        stop(sQuote("weights"), " contains real valued elements; currently
             only integer values are allowed") 


    where <- rep(0, object@nobs)
    storage.mode(where) <- "integer"

    ### grow the tree
    tree <- .Call("R_TreeGrow", object, weights, fitmem, controls, where,
                  PACKAGE = "party")

    ### create S3 classes and put names on lists
    tree <- prettytree(tree, names(object@inputs@variables), object@inputs@levels)

    ### prepare the returned object
    RET <- new("BinaryTree")
    RET@tree <- tree
    RET@where <- where
    RET@responses <- object@responses
    RET@inputnames <- names(object@inputs@variables)
    RET@levels <- object@inputs@levels

    ### get terminal node numbers
    RET@get_where <- function(newdata = NULL, mincriterion = 0, ...) {

        if (is.null(newdata) && mincriterion == 0) return(where)

        if (is.null(newdata)) {
            newinp <- object@inputs
        } else {
            penv <- new.env()
            object@menv@set("input", data = newdata, env = penv)
            newinp <- initVariableFrame(get("input", envir = penv))
        }

        .Call("R_get_nodeID", tree, newinp, mincriterion, PACKAGE = "party")
    }

    ### (estimated) conditional distribution of the response given the
    ### covariates
    RET@cond_distr_response <- function(newdata = NULL, mincriterion = 0, ...) { 
        
        wh <- RET@get_where(newdata = newdata, mincriterion = mincriterion)

        response <- object@responses

        ### survival: estimated Kaplan-Meier
        if (any(response@is_censored)) {
            swh <- sort(unique(wh))
            w <- .Call("R_getweights", tree, swh,
                       PACKAGE = "party")
            RET <- vector(mode = "list", length = length(wh))
            resp <- response@variables[[1]]
            for (i in 1:length(swh))
                RET[wh == swh[i]] <- list(survival:::survfit(resp,
                    weights = w[[i]], subset = w[[i]] > 0))
            return(RET)
        }

        ### classification: estimated class probabilities
        ### regression: the means, not really a distribution
        RET <- .Call("R_getpredictions", tree, wh, PACKAGE = "party")
        return(RET)
    }

    ### predict in the response space, always!
    RET@predict_response <- function(newdata = NULL, mincriterion = 0, ...) { 

        cdresp <- RET@cond_distr_response(newdata = newdata, 
                                          mincriterion = mincriterion, ...)

        response <- object@responses
        ### classification: classes
        if (all(response@is_nominal || response@is_ordinal)) {
            lev <- levels(response@variables[[1]])
            RET <- factor(lev[unlist(lapply(cdresp, which.max))],
                          levels = levels(response@variables[[1]]))
            return(RET)
        }

        ### survival: median survival time
        if (any(response@is_censored)) {
            RET <- sapply(cdresp, mst)
            return(RET)
        }

        ### regression: mean (median would be possible)
        RET <- unlist(cdresp)
        RET <- matrix(unlist(RET),
                      nrow = length(RET), byrow = TRUE)
        ### <FIXME> what about multivariate responses?
        if (response@ninputs == 1)
            colnames(RET) <- names(response@variables)
        ### </FIXME>
        return(RET)
    }

    RET@prediction_weights <- function(newdata = NULL, mincriterion = 0, ...) {

        wh <- RET@get_where(newdata = newdata, mincriterion = mincriterion)

        swh <- sort(unique(wh))
        w <- .Call("R_getweights", tree, swh,
                   PACKAGE = "party")
        RET <- vector(mode = "list", length = length(wh))   
        
        for (i in 1:length(swh))
            RET[wh == swh[i]] <- list(w[[i]])
        return(RET)
    }
    return(RET)
}

ctreedpp <- function(formula, data = list(), subset = NULL, 
                     na.action = NULL, rfun = NULL, ifun = NULL, ...) {

    dat <- ModelEnvFormula(formula = formula, data = data, 
                           subset = subset, ...)
    inp <- initVariableFrame(dat@get("input"), fun = ifun)

    if (has(dat, "censored"))
        censored <- as.logical(dat@get("censored")[[1]])
    response <- dat@get("response")

    if (any(is.na(response)))
        stop("missing values in response variable not allowed")

    if (length(response) != 1 && has(dat, "censored"))
        stop("multiple failure times not yet implemented")

    if (has(dat, "censored")) {
        if (!is.null(rfun))
            resp <- initVariableFrame(Surv(response[[1]], censored), 
                                      fun = rfun)
        else
            resp <- initVariableFrame(Surv(response[[1]], censored))
    } else {
        resp <- initVariableFrame(response, fun = rfun)
    }
    RET <- new("LearningSample", inputs = inp, responses = resp,
               weights = rep(1, inp@nobs), nobs = inp@nobs,
               ninputs = inp@ninputs, menv = dat)
    RET
}


conditionalTree <- new("StatModel",
                   capabilities = new("StatModelCapabilities"),
                   name = "unbiased conditional recursive partitioning",
                   dpp = ctreedpp,
                   fit = ctreefit,
                   predict = function(object, ...) 
                                 object@predict_response(...) )

treecontrols <- function(teststattype = "quadform", 
                         testtype = "Bonferroni", 
                         nresample = 9999, mincriterion = 0.95, 
                         stump = FALSE, minsplit = 20, maxsurrogate = 0) {

    RET <- new("TreeControl")
    if (teststattype %in% levels(RET@varctrl@teststattype)) {
        RET@varctrl@teststattype <- factor(teststattype, 
            levels = levels(RET@varctrl@teststattype))
    } else {
        stop(teststattype, " not defined")
    }

    if (testtype %in% levels(RET@gtctrl@testtype))
        RET@gtctrl@testtype <- factor(testtype, 
            levels = levels(RET@gtctrl@testtype))
    else
        stop(testtype, " not defined")

    if (RET@gtctrl@testtype == "MonteCarlo") RET@varctrl@pvalue <- FALSE
    RET@gtctrl@nresample <- as.integer(nresample)
    RET@gtctrl@mincriterion <- mincriterion
    RET@splitctrl@minsplit <- minsplit
    RET@splitctrl@maxsurrogate <- as.integer(maxsurrogate)
    RET@tgctrl@stump <- stump
    RET
}

ctree <- function(formula, data, subset = NULL, weights = NULL, 
                  teststattype = c("quadform", "maxabs"),
                  testtype = c("Bonferroni", "MonteCarlo", "Raw"),
                  mincriterion = 0.95, minsplit = 20, stump = FALSE,
                  nresample = 9999, maxsurrogate = 0, 
                  ifun = NULL, rfun = NULL) {

    teststattype <- match.arg(teststattype)
    testtype <- match.arg(testtype)
    ctrl <- treecontrols(teststattype, testtype, stump = stump, nresample =
                         nresample, maxsurrogate = maxsurrogate, 
                         mincriterion = mincriterion, minsplit = minsplit)
    ls <- conditionalTree@dpp(formula, data, subset, ifun = ifun, rfun = rfun)
    conditionalTree@fit(ls, ctrl, weights = weights)
}
