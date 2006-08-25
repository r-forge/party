
/**
    Functions for variable selection in each node of a tree
    *\file IndependenceTest.c
    *\author $Author$
    *\date $Date$
*/
                
#include "party.h"


/**
    Computes the test statistic and, if requested, the corresponding 
    P-value for a linear statistic \n
    *\param linexpcov an object of class `LinStatExpectCovar'
    *\param varctrl an object of class `VariableControl'
    *\param ans_teststat; return value, the test statistic
    *\param ans_pvalue; return value, the p-value
*/
        
void C_TeststatPvalue(const SEXP linexpcov, const SEXP varctrl, 
                      double *ans_teststat, double *ans_pvalue) {
    
    double releps, abseps, tol;
    int maxpts;
    
    maxpts = get_maxpts(varctrl);
    tol = get_tol(varctrl);
    abseps = get_abseps(varctrl);
    releps = get_releps(varctrl);
    
    /* compute the test statistic */
    ans_teststat[0] = C_TestStatistic(linexpcov, get_teststat(varctrl), 
                                  get_tol(varctrl));

    /* compute the p-value if requested */                                  
    if (get_pvalue(varctrl))
        ans_pvalue[0] =  C_ConditionalPvalue(ans_teststat[0], linexpcov, 
                                         get_teststat(varctrl),
                                         tol, &maxpts, &releps, &abseps);
    else
        ans_pvalue[0] = 1.0;
}

/**
    Computes the test statistic and the node criterion \n
    *\param linexpcov an object of class `LinStatExpectCovar'
    *\param varctrl an object of class `VariableControl'
    *\param ans_teststat; return value, the test statistic
    *\param ans_criterion; return value, thep-value
*/
        
void C_TeststatCriterion(const SEXP linexpcov, const SEXP varctrl, 
                         double *ans_teststat, double *ans_criterion) {
    
    C_TeststatPvalue(linexpcov, varctrl, ans_teststat, ans_criterion);
    
    /* the node criterion is to be MAXIMISED, 
       i.e. 1-pvalue or test statistic \in \[0, \infty\] */
    if (get_pvalue(varctrl))
        ans_criterion[0] = 1 - ans_criterion[0];
    else
        ans_criterion[0] = ans_teststat[0];
    
}


/**
    Test of independence between x and y \n
    *\param x values of the transformation
    *\param y values of the influence function
    *\param weights case weights
    *\param linexpcov an object of class `VariableControl' for T
    *\param varctrl an object of class `VariableControl'
    *\param ans; return value, a double vector (teststat, pvalue)
*/

void C_IndependenceTest(const SEXP x, const SEXP y, const SEXP weights, 
                        SEXP linexpcov, SEXP varctrl, 
                        SEXP ans) {
    
    /* compute linear statistic and its conditional expectation and
       covariance
    */
    C_LinStatExpCov(REAL(x), ncol(x), REAL(y), ncol(y), 
                    REAL(weights), nrow(x), 1, 
                    GET_SLOT(linexpcov, PL2_expcovinfSym), linexpcov);

    /* compute test statistic */
    if (get_teststat(varctrl) == 2) 
        C_LinStatExpCovMPinv(linexpcov, get_tol(varctrl));
    C_TeststatPvalue(linexpcov, varctrl, &REAL(ans)[0], &REAL(ans)[1]);
}


/**
    R-interface to C_IndependenceTest \n
    *\param x values of the transformation
    *\param y values of the influence function
    *\param weights case weights
    *\param linexpcov an object of class `VariableControl' for T
    *\param varctrl an object of class `VariableControl'
*/

SEXP R_IndependenceTest(SEXP x, SEXP y, SEXP weights, SEXP linexpcov, SEXP varctrl) {
                        
    SEXP ans;
    
    PROTECT(ans = allocVector(REALSXP, 2));
    C_IndependenceTest(x, y, weights, linexpcov, varctrl, ans);
    UNPROTECT(1);
    return(ans);
}


/**
    Perform a global test on independence of a response and multiple inputs \n
    *\param learnsample an object of class `LearningSample'
    *\param weights case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param varctrl an object of class `VariableControl'
    *\param gtctrl an object of class `GlobalTestControl'
    *\param minsplit minimum sum of weights to proceed
    *\param ans_teststat return value; vector of test statistics
    *\param ans_criterion return value; vector of node criteria 
            (adjusted) pvalues or raw test statistics
*/

void C_GlobalTest(const SEXP learnsample, const SEXP weights, 
                  SEXP fitmem, const SEXP varctrl, 
                  const SEXP gtctrl, const double minsplit, 
                  double *ans_teststat, double *ans_criterion) {

    int ninputs, nobs, j, i, k, RECALC = 1, type;
    SEXP responses, inputs, y, x, xmem, expcovinf;
    SEXP thiswhichNA;
    double *thisweights, *dweights, *pvaltmp;
    int *ithiswhichNA, RANDOM, mtry, *randomvar, *index;
    int *dontuse, *dontusetmp;
    
    ninputs = get_ninputs(learnsample);
    nobs = get_nobs(learnsample);
    responses = GET_SLOT(learnsample, PL2_responsesSym);
    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    dweights = REAL(weights);
    
    y = get_transformation(responses, 1);
    
    expcovinf = GET_SLOT(fitmem, PL2_expcovinfSym);
    C_ExpectCovarInfluence(REAL(y), ncol(y), REAL(weights), 
                           nobs, expcovinf);
    
    if (REAL(GET_SLOT(expcovinf, PL2_sumweightsSym))[0] < minsplit) {
        for (j = 0; j < ninputs; j++) {
            ans_teststat[j] = 0.0;
            ans_criterion[j] = 0.0;
        }
    } else {

        dontuse = INTEGER(get_dontuse(fitmem));
        dontusetmp = INTEGER(get_dontusetmp(fitmem));
    
        for (j = 0; j < ninputs; j++) dontusetmp[j] = !dontuse[j];
    
        /* random forest */
        RANDOM = get_randomsplits(gtctrl);
        mtry = get_mtry(gtctrl);
        if (RANDOM & (mtry > ninputs)) {
            warning("mtry is larger than ninputs, using mtry = inputs");
            mtry = ninputs;
            RANDOM = 0;
        }
        if (RANDOM) {
            index = Calloc(ninputs, int);
            randomvar = Calloc(mtry, int);
            GetRNGstate();
            C_SampleNoReplace(index, ninputs, mtry, randomvar);
            PutRNGstate();
            j = 0;
            for (k = 0; k < mtry; k++) {
                j = randomvar[k];
                while(dontuse[j] && j < ninputs) j++;
                if (j == ninputs) 
                    error("not enough variables to sample from");
                dontusetmp[j] = 0;
            }
            Free(index);
            Free(randomvar);
        }

        for (j = 1; j <= ninputs; j++) {

            if ((RANDOM && dontusetmp[j - 1]) || dontuse[j - 1]) {
                ans_teststat[j - 1] = 0.0;
                ans_criterion[j - 1] = 0.0;
                continue; 
            }
        
            x = get_transformation(inputs, j);

            xmem = get_varmemory(fitmem, j);
            if (!has_missings(inputs, j)) {
                C_LinStatExpCov(REAL(x), ncol(x), REAL(y), ncol(y),
                                REAL(weights), nrow(x), !RECALC, expcovinf,
                                xmem);
            } else {
                thisweights = REAL(get_weights(fitmem, j));
                thiswhichNA = get_missings(inputs, j);
                ithiswhichNA = INTEGER(thiswhichNA);
                for (i = 0; i < nobs; i++) thisweights[i] = dweights[i];
                for (k = 0; k < LENGTH(thiswhichNA); k++)
                    thisweights[ithiswhichNA[k] - 1] = 0.0;
                C_LinStatExpCov(REAL(x), ncol(x), REAL(y), ncol(y),
                                thisweights, nrow(x), RECALC, 
                                GET_SLOT(xmem, PL2_expcovinfSym),
                                xmem);
            }

            if (get_teststat(varctrl) == 2)
                C_LinStatExpCovMPinv(xmem, get_tol(varctrl));
            C_TeststatCriterion(xmem, varctrl, &ans_teststat[j - 1], 
                                &ans_criterion[j - 1]);
        }                

        type = get_testtype(gtctrl);
        switch(type) {
            /* Bonferroni: p_adj = 1 - (1 - p)^k */
            case BONFERRONI: 
                    for (j = 0; j < ninputs; j++)
                        ans_criterion[j] = R_pow_di(ans_criterion[j], ninputs);
                    break;
            /* Monte-Carlo */
            case MONTECARLO: 
                    pvaltmp = Calloc(ninputs, double);
                    C_MonteCarlo(ans_criterion, learnsample, weights, fitmem, 
                                 varctrl, gtctrl, pvaltmp);
                    for (j = 0; j < ninputs; j++)
                        ans_criterion[j] = 1 - pvaltmp[j];
                    Free(pvaltmp);
                    break;
            /* aggregated */
            case AGGREGATED: 
                    error("C_GlobalTest: aggregated global test not yet implemented");
                    break;
            /* raw */
            case UNIVARIATE: break;
            case TESTSTATISTIC: break;
            default: error("C_GlobalTest: undefined value for type argument");
                     break;
        }
    }
}


/**
    R-interface to C_GlobalTest \n
    *\param learnsample an object of class `LearningSample'
    *\param weights case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param varctrl an object of class `VariableControl'
    *\param gtctrl an object of class `GlobalTestControl'
*/

SEXP R_GlobalTest(SEXP learnsample, SEXP weights, SEXP fitmem, 
                  SEXP varctrl, SEXP gtctrl) {

    SEXP ans, teststat, criterion;
    
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, 
        teststat = allocVector(REALSXP, get_ninputs(learnsample)));
    SET_VECTOR_ELT(ans, 1, 
        criterion = allocVector(REALSXP, get_ninputs(learnsample)));

    C_GlobalTest(learnsample, weights, fitmem, varctrl, gtctrl, 0, 
                 REAL(teststat), REAL(criterion));
    
    UNPROTECT(1);
    return(ans);
}
