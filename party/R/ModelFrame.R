
setClass("ModelMatrix", representation = 
    representation(		
        responseMatrix = "MATRIXorNULL",
        designMatrix = "MATRIXorNULL"), 
    validity = function(object) {
             NULLcheck = !is.null(object@responseMatrix) &&
                             !is.null(object@designMatrix)
             if (NULLcheck) 
                 return(nrow(responseMatrix) == nrow(designMatrix))
             else
                 return(TRUE)
    }
)

setClass("ModelFrame", representation = 
    representation(
        response = "DForNULL",
        input    = "DForNULL",
        blocks   = "FACTORorNULL",	# interactions of 2 or more factors???
        cweights = "INTEGERorNULL",
        censored = "logical", 
        blockobs = "logical"),
    prototype = list(response = NULL, input = NULL, cweights = NULL, 
                     censored = FALSE, blockobs = FALSE),  
    contains = "ModelMatrix", 
    validity = function(object) {
        if (is.null(object@response) && is.null(object@input)) {
            warning("object contains not data")
            return(TRUE)
        }
        if (is.null(object@response)) {
            warning("at least response must be specified")
            return(FALSE)
        }
        n = nrow(object@response) 

        if (!is.null(object@input)) {
            if (nrow(object@input) != n) {
                warning("sample size of", sQuote("input"), "does not match")
                return(FALSE)
            }
        }

        if (!is.null(object@blocks)) {
            if (length(object@blocks) != n && !object@blockobs) {
                warning("sample size of", sQuote("blocks"), "does not match")
                return(FALSE)
            }
            if (any(is.na(object@blocks))) {
                warning(sQuote("blocks"), " contains NA's")
                return(FALSE)
            }

        } 
        if (!is.null(object@cweights)) {
            if (length(object@cweights) != n)  {
                warning("sample size of", sQuote("cweights"), "does not match")
                return(FALSE)
            }
            if (any(object@cweights < 0)) {
                warning("values of ", sQuote("cweights"), "must be positive integers")
                return(FALSE)
            }
            if (any(is.na(object@cweights))) {
                warning(sQuote("cweights"), " contains NA's")
                return(FALSE)
            }

        } 
        if (!is.null(object@censored)) {
            if (any(object@censored) && length(object@censored) != n) {
                warning("sample size of", sQuote("censored"), "does not match")
                return(FALSE)
            }
            if (any(is.na(object@censored))) {
                warning(sQuote("censored"), " contains NA's")
                return(FALSE)
            }
        }
        return(TRUE)
    }
)

setClass("ModelFormulae", representation = 
    representation(
        response = "FORMULAorNULL",
        input = "FORMULAorNULL",
        blocks = "FORMULAorNULL",
        weights = "FORMULAorNULL",
        censored = "FORMULAorNULL")
)


setClass("ModelFrameFormula", representation(
         formulae  = "ModelFormulae"), contains = "ModelFrame")

setGeneric("ModelFrame", function(object, ...) standardGeneric("ModelFrame"))

ModelFrameFormula = function(formula, data = list(), subset = NULL, 
                             na.action = NULL, frame = NULL, ...) {
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    if (is.null(subset)) mf$subset <- NULL

    ### NA-handling will for the ModelFrame objects later on...
    mf$na.action = na.pass

    MF = new("ModelFrameFormula")

    if (!missing(formula)) {
        formula = terms(formula, data = data)
        attributes(formula) = NULL
        if (length(formula) == 3) {
            fresponse = formula[c(1,2)]
            frhs = formula[c(1,3)]
            if (frhs[[2]] == "1")
                frhs = NULL
        }
        if (length(formula) == 2) {
            fresponse = NULL   
            frhs = formula
        }
        finput = frhs
        fblocks = frhs


        # <FIXME>
        # will fail for `y ~ . | blocks' constructs
        # </FIXME>

        if (!is.null(frhs) && length(frhs[[2]]) > 1) {
            if (deparse(frhs[[2]][[1]]) == "|") {
                finput[[2]] = frhs[[2]][[2]]
                fblocks[[2]] = frhs[[2]][[3]]
            } else {
                fblocks = NULL
            }
        } else {
            fblocks = NULL
        }

        fcensored = NULL

        if (!is.null(fresponse) && length(fresponse[[2]]) == 3) {
            if (fresponse[[2]][[1]] == "Surv") {
                fcensored = formula(paste("~", fresponse[[2]][[3]]))
                fresponse = formula(paste("~", fresponse[[2]][[2]])) 
            }
        }

        MF@formulae = new("ModelFormulae", response = formula(fresponse), 
                                           input    = formula(finput),
                                           blocks   = formula(fblocks),
                                           censored = formula(fcensored),
                                           weights  = NULL)

        if (is.null(frame)) frame = parent.frame()


        if (!is.null(fresponse)) {
            mf$formula = fresponse
            MF@response = eval(mf, frame)
        }

        if (!is.null(fcensored)) {
            mf$formula = fcensored
            censored = eval(mf, frame)
            if (ncol(censored) != 1) 
                stop("more than one censoring variable")
            else 
                censored = censored[,1]
            if (!is.logical(censored)) censored = as.logical(censored)
            MF@censored = censored
        }

        if (!is.null(finput)) {
            mf$formula = finput
            MF@input <- eval(mf, frame)   
        }

        if (!is.null(fblocks)) {
            mf$formula = fblocks
            blocks <- eval(mf, frame)
            if (ncol(blocks) != 1) 
                stop("more than one block factor currently not allowed")
            else 
                blocks = blocks[[1]]
            if (all(class(blocks) %in% c("AsIs", "factor")))
                class(blocks) = "factor"
            MF@blocks = blocks
            MF@blockobs = (nlevels(MF@blocks) == nrow(MF@response))
        }
    }
    
    ### handle NA's
    if (!is.null(na.action))
        MF = na.action(MF)
    MF
}

setMethod("ModelFrame", "formula", function(object, ...)
    ModelFrameFormula(formula = object, ...))

ModelFrameTable = function(object, paired = FALSE) {
    tdim = dim(object)
    if (!any(length(tdim) == c(2,3))) 
        stop("cannot handle tables with dimensions unequal two or three")
    dftab = as.data.frame(object)

    if (paired) {
        cn = colnames(dftab)
        pxtab = c()
        for (i in 1:nrow(dftab)) {
            pxtab = rbind(pxtab, c(1, dftab[i,1], dftab[i,3]),
                                 c(2, dftab[i,2], dftab[i,3]))
        }
        pxtab = as.data.frame(pxtab)
        pxtab[,1] = factor(pxtab[,1], labels = cn[1:2])
        pxtab[,2] = factor(pxtab[,2], labels = levels(dftab[,2]))
        dftab = pxtab
        colnames(dftab) = c(paste("V", 1:2, sep=""), "Freq")
        blocks = gl(nrow(dftab)/2, 2)
        MF = new("ModelFrame", response = dftab[,2,drop = FALSE],
                 input = dftab[,1,drop = FALSE], blocks = blocks, 
                 blockobs = TRUE, cweights = as.integer(dftab[["Freq"]]))

    } else {
        if (length(tdim) == 2) 
            MF = new("ModelFrame", response = dftab[,1:2], 
                     cweights = as.integer(dftab[["Freq"]]))
        else
            MF = new("ModelFrame", response = dftab[,1:2], blocks = dftab[,3],
                     cweights = as.integer(dftab[["Freq"]]))
    }
    return(MF)
}

setMethod("ModelFrame", "table", function(object, paired = FALSE) 
    ModelFrameTable(object = object, paired = paired))
        

setMethod("subset", "ModelFrame",
    function(x, i, drop=FALSE)
    {
    if (!is.null(x@responseMatrix))
        x@responseMatrix = x@responseMatrix[i,,drop = drop]
    if (!is.null(x@designMatrix))
        x@designMatrix = x@designMatrix[i,,drop = drop]
    if (!is.null(x@response))
        x@response = x@response[i,,drop = drop]
    if (!is.null(x@input))
        x@input = x@input[i,,drop = drop]
    if (!is.null(x@blocks))
        x@blocks = x@blocks[i]
    if (!is.null(x@cweights))
        x@cweights = x@cweights[i]
    if (!is.null(x@censored) && length(x@censored) > 1)
        x@censored = x@censored[i]
    x
}
)


setMethod("na.fail", "ModelFrame", function(object, ...) {
    ccresponse = complete.cases(object@response)
    ccinput = complete.cases(object@input)
    if(all(ccresponse && ccinput)) object else stop("missing values in object");
})

setMethod("na.fail", "ModelFrame", function(object, ...) {
    ccresponse = complete.cases(object@response)
    ccinput = complete.cases(object@input)
    if(all(ccresponse && ccinput)) object else stop("missing values in object");
})

setMethod("na.omit", "ModelFrame", function(object, ...) {
    ### <FIXME> this is just a prototype, we need na.action attributes
    ### the way na.omit.default etc. implement
    ccresponse = complete.cases(object@response)
    ccinput = complete.cases(object@input)
    cc = (ccresponse & ccinput)
    subset(object, cc)
    ### </FIXME>
})

setMethod("show", "ModelFrame", function(object) {
    cat("\n")
    cat("\t A ModelFrame with \n\n")
    n = NULL
    if (!is.null(object@response)) {
        cat("\t Response Variable(s): ", colnames(object@response), "\n")
        n = nrow(object@response)
    }
    if (!is.null(object@input)) {
        cat("\t Input Variable(s): ", colnames(object@input), "\n")
        n = nrow(object@input)
    }
    if (is.null(n)) 
        cat("\t no observations\n")
    else
        cat("\t with n = ", n, " observations\n")
})
