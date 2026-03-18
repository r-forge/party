
/**
    Functions for variable selection in each node of a tree
    *\file IndependenceTest.c
    *\author $Author$
    *\date $Date$
*/
                
#include "party.h"
#include <libcoinAPI.h>

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
    /* GET_SLOT and ncol are assumed NOT to return a fresh object so
       we don't PROTECT here */
    C_LinStatExpCov(REAL(x), ncol(x), REAL(y), ncol(y), 
                    REAL(weights), nrow(x), 1, 
                    PROTECT(GET_SLOT(linexpcov, PL2_expcovinfSym)), linexpcov);
    UNPROTECT(1);

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
    *\param depth an integer giving the depth of the current node
*/

void C_GlobalTest2(const SEXP learnsample, const SEXP weights, 
                  SEXP fitmem, const SEXP varctrl, 
                  const SEXP gtctrl, const double minsplit, 
                  double *ans_teststat, double *ans_criterion, int depth) {

    int ninputs, nobs, j, i, k, RECALC = 1, type;
    SEXP responses, inputs, y, x, xmem, expcovinf;
    SEXP Smtry;
    double *thisweights, *pvaltmp, stweights = 0.0;
    int RANDOM, mtry, *randomvar, *index;
    int *dontuse, *dontusetmp, countvars = 0;
    
    ninputs = get_ninputs(learnsample);
    nobs = get_nobs(learnsample);
    responses = GET_SLOT(learnsample, PL2_responsesSym);
    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    
    /* y = get_transformation(responses, 1); */
    y = get_test_trafo(responses);
    
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
        Smtry = get_mtry(gtctrl);
        if (LENGTH(Smtry) == 1) {
            mtry = INTEGER(Smtry)[0];
        } else {
            /* mtry may vary with tree depth */
            depth = (depth <= LENGTH(Smtry)) ? depth : LENGTH(Smtry);
            mtry = INTEGER(get_mtry(gtctrl))[depth - 1];
        }
        if (RANDOM & (mtry > ninputs)) {
            warning("mtry is larger than ninputs, using mtry = inputs");
            mtry = ninputs;
            RANDOM = 0;
        }
        if (RANDOM) {
            index = R_Calloc(ninputs, int);
            randomvar = R_Calloc(mtry, int);
            C_SampleNoReplace(index, ninputs, mtry, randomvar);
            j = 0;
            for (k = 0; k < mtry; k++) {
                j = randomvar[k];
                while(dontuse[j] && j < ninputs) j++;
                if (j == ninputs) 
                    error("not enough variables to sample from");
                dontusetmp[j] = 0;
            }
            R_Free(index);
            R_Free(randomvar);
        }

        countvars = 0;
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
                thisweights = C_tempweights(j, REAL(weights), fitmem, inputs);

                /* check if minsplit criterion is still met 
                   in the presence of missing values
                   bug spotted by Han Lee <Han.Lee@geodecapital.com>
                       fixed 2006-08-31
                */
                stweights = 0.0;
                for (i = 0; i < nobs; i++) stweights += thisweights[i];
                if (stweights < minsplit) {
                    ans_teststat[j - 1] = 0.0;
                    ans_criterion[j - 1] = 0.0;
                    continue; 
                }

               /* GET_SLOT is assumed NOT to return a fresh object so
                  we don't PROTECT here */
                C_LinStatExpCov(REAL(x), ncol(x), REAL(y), ncol(y),
                                thisweights, nrow(x), RECALC, 
                                PROTECT(GET_SLOT(xmem, PL2_expcovinfSym)),
                                xmem);
                UNPROTECT(1);
            }

            /* count the number of non-constant variables */
            countvars++;

            /* teststat = "quad" 
               ATTENTION: we mess with xmem by putting elements with zero variances last
                          but C_Node() reuses the original xmem for setting up 
                          categorical splits */
            if (get_teststat(varctrl) == 2)
                C_LinStatExpCovMPinv(xmem, get_tol(varctrl));
                
            C_TeststatCriterion(xmem, varctrl, &ans_teststat[j - 1], 
                                &ans_criterion[j - 1]);

            /* teststat = "quad"
               make sure that the test statistics etc match the original order of levels 
               <FIXME> can we avoid to compute these things twice??? */ 
            if (get_teststat(varctrl) == 2) {
                if (!has_missings(inputs, j)) {
                    /* ncol is assumed NOT to return a fresh object so
                       we don't PROTECT here */
                    C_LinStatExpCov(REAL(x), ncol(x), REAL(y), ncol(y),
                                    REAL(weights), nrow(x), !RECALC, expcovinf,
                                    xmem);
                } else {
                    C_LinStatExpCov(REAL(x), ncol(x), REAL(y), ncol(y),
                                    thisweights, nrow(x), RECALC, 
                                    PROTECT(GET_SLOT(xmem, PL2_expcovinfSym)),
                                    xmem);
                    UNPROTECT(1);
                }
            }
            /* </FIXME> */
        }                

        type = get_testtype(gtctrl);
        switch(type) {
            /* Bonferroni: p_adj = 1 - (1 - p)^k */
            case BONFERRONI: 
                    for (j = 0; j < ninputs; j++)
                        ans_criterion[j] = R_pow_di(ans_criterion[j], countvars);
                    break;
            /* Monte-Carlo */
            case MONTECARLO: 
                    pvaltmp = R_Calloc(ninputs, double);
                    C_MonteCarlo(ans_criterion, learnsample, weights, fitmem, 
                                 varctrl, gtctrl, pvaltmp);
                    for (j = 0; j < ninputs; j++)
                        ans_criterion[j] = 1 - pvaltmp[j];
                    R_Free(pvaltmp);
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
    *\param depth an integer giving the depth of the current node
*/

void C_GlobalTest(const SEXP learnsample, const SEXP weights, 
                  SEXP fitmem, const SEXP varctrl, 
                  const SEXP gtctrl, const double minsplit, 
                  double *ans_teststat, double *ans_criterion, int depth) {

    int ninputs, nobs, j, i, k, q, n, l, type;
    SEXP responses, inputs, y, x, xmem, expcovinf, expcovinfj, tmp;
    SEXP Smtry;

    SEXP LECV, PLS, iEMPTY, iFALSE, iTRUE, ans;

    double *thisweights, stweights = 0.0;
    int RANDOM, mtry, *randomvar, *index;
    int *dontuse, *dontusetmp, countvars = 0;
    
    ninputs = get_ninputs(learnsample);
    nobs = get_nobs(learnsample);

    PROTECT(responses = GET_SLOT(learnsample, PL2_responsesSym));
    PROTECT(inputs = GET_SLOT(learnsample, PL2_inputsSym));
    
    PROTECT(iEMPTY = allocVector(INTSXP, 0));
    PROTECT(iTRUE = allocVector(INTSXP, 1));
    INTEGER(iTRUE)[0] = 1;
    PROTECT(iFALSE = allocVector(INTSXP, 1));
    INTEGER(iFALSE)[0] = 0;
    
    /* y = get_transformation(responses, 1); */
    PROTECT(y = get_test_trafo(responses));

    thisweights = REAL(weights);
    for (i = 0; i < nobs; i++) stweights += thisweights[i];

    expcovinf = GET_SLOT(fitmem, PL2_expcovinfSym);
    REAL(GET_SLOT(expcovinf, PL2_sumweightsSym))[0] = stweights;
    
    if (stweights < minsplit) {
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
        Smtry = get_mtry(gtctrl);
        if (LENGTH(Smtry) == 1) {
            mtry = INTEGER(Smtry)[0];
        } else {
            /* mtry may vary with tree depth */
            depth = (depth <= LENGTH(Smtry)) ? depth : LENGTH(Smtry);
            mtry = INTEGER(get_mtry(gtctrl))[depth - 1];
        }
        if (RANDOM & (mtry > ninputs)) {
            warning("mtry is larger than ninputs, using mtry = inputs");
            mtry = ninputs;
            RANDOM = 0;
        }
        if (RANDOM) {
            index = R_Calloc(ninputs, int);
            randomvar = R_Calloc(mtry, int);
            C_SampleNoReplace(index, ninputs, mtry, randomvar);
            j = 0;
            for (k = 0; k < mtry; k++) {
                j = randomvar[k];
                while(dontuse[j] && j < ninputs) j++;
                if (j == ninputs) 
                    error("not enough variables to sample from");
                dontusetmp[j] = 0;
            }
            R_Free(index);
            R_Free(randomvar);
        }

        countvars = 0;
        for (j = 1; j <= ninputs; j++) {

            if ((RANDOM && dontusetmp[j - 1]) || dontuse[j - 1]) {
                ans_teststat[j - 1] = 0.0;
                ans_criterion[j - 1] = 0.0;
                continue; 
            }
        
            PROTECT(x = get_transformation(inputs, j));
            PROTECT(xmem = get_varmemory(fitmem, j));
            expcovinfj = GET_SLOT(xmem, PL2_expcovinfSym);

            if (!has_missings(inputs, j)) {
                PROTECT(LECV = libcoin_R_ExpectationCovarianceStatistic(x, y, weights, 
                        iEMPTY, // SEXP subset
                        iEMPTY, // SEXP block 
                        iFALSE, // SEXP varonly 
                        GET_SLOT(gtctrl, PL2_tolSym)));

                if (get_testtype(gtctrl) == MONTECARLO)
                    SET_VECTOR_ELT(LECV, 12, // PermutedLinearStatistic_SLOT
                    libcoin_R_PermutedLinearStatistic(x, y, weights, 
                        iEMPTY, // SEXP subset
                        iEMPTY, // SEXP block
                        GET_SLOT(gtctrl, PL2_nresampleSym)));
            } else {

                PROTECT(LECV = libcoin_R_ExpectationCovarianceStatistic(x, y, weights, 
                    get_subset(inputs, j), // SEXP subset
                    iEMPTY, // SEXP block
                    iFALSE, // SEXP varonly
                    GET_SLOT(gtctrl, PL2_tolSym)));

                if (get_testtype(gtctrl) == MONTECARLO)
                    SET_VECTOR_ELT(LECV, 12, // PermutedLinearStatistic_SLOT
                    libcoin_R_PermutedLinearStatistic(x, y, weights, 
                        get_subset(inputs, j), // SEXP subset
                        iEMPTY, // SEXP block
                        GET_SLOT(gtctrl, PL2_nresampleSym)));


                stweights = REAL(VECTOR_ELT(LECV, 15))[0]; // Sumweights_SLOT
                if (stweights < minsplit) {
                    ans_teststat[j - 1] = 0.0;
                    ans_criterion[j - 1] = 0.0;
                    continue; 
                }
            }

            tmp = VECTOR_ELT(LECV, 7); // ExpectationInfluence_SLOT
            for (q = 0; q < LENGTH(tmp); q++) {
                REAL(GET_SLOT(expcovinf, PL2_expectationSym))[q] = REAL(tmp)[q];
                REAL(GET_SLOT(expcovinfj, PL2_expectationSym))[q] = REAL(tmp)[q];
            }

            tmp = VECTOR_ELT(LECV, 8); // CovarianceInfluence_SLOT, packed 
            n = (int) (sqrt(0.25 + 2 * LENGTH(tmp)) - 0.5);
            k = 0;
            for (i = 0; i < n; i++) {
                REAL(GET_SLOT(expcovinf, PL2_covarianceSym))[i * n + i] = REAL(tmp)[k];     /* diagonal */
                REAL(GET_SLOT(expcovinfj, PL2_covarianceSym))[i * n + i] = REAL(tmp)[k];     /* diagonal */
                k++;
                for (l = i + 1; l < n; l++) {
                    REAL(GET_SLOT(expcovinf, PL2_covarianceSym))[i * n + l] = REAL(tmp)[k]; /* lower triangular */
                    REAL(GET_SLOT(expcovinf, PL2_covarianceSym))[l * n + i] = REAL(tmp)[k]; /* upper triangular */
                    REAL(GET_SLOT(expcovinfj, PL2_covarianceSym))[i * n + l] = REAL(tmp)[k]; /* lower triangular */
                    REAL(GET_SLOT(expcovinfj, PL2_covarianceSym))[l * n + i] = REAL(tmp)[k]; /* upper triangular */
                    k++;
                }
            }

            tmp = VECTOR_ELT(LECV, 15); // Sumweights_SLOT
            REAL(GET_SLOT(expcovinf, PL2_sumweightsSym))[0] = REAL(tmp)[0];

            tmp = VECTOR_ELT(LECV, 0); // LinearStatistic_SLOT
            for (q = 0; q < LENGTH(tmp); q++) {
                REAL(GET_SLOT(xmem, PL2_linearstatisticSym))[q] = REAL(tmp)[q];
            }

            tmp = VECTOR_ELT(LECV, 1); // Expectation_SLOT
            for (q = 0; q < LENGTH(tmp); q++) {
                REAL(GET_SLOT(xmem, PL2_expectationSym))[q] = REAL(tmp)[q];
            }

            tmp = VECTOR_ELT(LECV, 2); // Covariance_SLOT, packed 
            n = (int) (sqrt(0.25 + 2 * LENGTH(tmp)) - 0.5);
            k = 0;
            for (i = 0; i < n; i++) {
                REAL(GET_SLOT(xmem, PL2_covarianceSym))[i * n + i] = REAL(tmp)[k];     /* diagonal */
                k++;
                for (l = i + 1; l < n; l++) {
                    REAL(GET_SLOT(xmem, PL2_covarianceSym))[i * n + l] = REAL(tmp)[k]; /* lower triangular */
                    REAL(GET_SLOT(xmem, PL2_covarianceSym))[l * n + i] = REAL(tmp)[k]; /* upper triangular */
                    k++;
                }
            }

            INTEGER(GET_SLOT(xmem, PL2_dimensionSym))[0] = n;

            tmp = VECTOR_ELT(LECV, 15); // Sumweights_SLOT
            REAL(GET_SLOT(expcovinf, PL2_sumweightsSym))[0] = REAL(tmp)[0];

            /* count the number of non-constant variables */
            countvars++;

            if (get_teststat(varctrl) == 2)
                PROTECT(ans = libcoin_R_QuadraticTest(LECV, 
                    GET_SLOT(varctrl, PL2_pvalueSym), 
                    iTRUE,  // SEXP lower
                    iFALSE, // SEXP give_log
                    iFALSE  // SEXP PermutedStatistics 
                    ));
            else
                PROTECT(ans = libcoin_R_MaximumTest(LECV,
                    iTRUE, // SEXP alternative (1 = "two.sided") 
                    GET_SLOT(varctrl, PL2_pvalueSym), 
                    iTRUE, // SEXP lower
                    iFALSE, // SEXP give_log
                    iFALSE, // SEXP PermutedStatistics
                    GET_SLOT(varctrl, PL2_maxptsSym), 
                    GET_SLOT(varctrl, PL2_relepsSym),
                    GET_SLOT(varctrl, PL2_absepsSym)));

            ans_teststat[j - 1] = REAL(VECTOR_ELT(ans, 0))[0]; 
            if (get_pvalue(varctrl))
                ans_criterion[j - 1] = REAL(VECTOR_ELT(ans, 1))[0]; // 1 - pvalue
            else
                ans_criterion[j - 1] = ans_teststat[j - 1];

            UNPROTECT(4);
        }                

        type = get_testtype(gtctrl);
        switch(type) {
            /* Bonferroni: p_adj = 1 - (1 - p)^k */
            case BONFERRONI: 
                    for (j = 0; j < ninputs; j++)
                        ans_criterion[j] = R_pow_di(ans_criterion[j], countvars);
                    break;
            /* Monte-Carlo */
            case MONTECARLO: 
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
    UNPROTECT(6);
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

    GetRNGstate();

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, 
        teststat = allocVector(REALSXP, get_ninputs(learnsample)));
    SET_VECTOR_ELT(ans, 1, 
        criterion = allocVector(REALSXP, get_ninputs(learnsample)));

    C_GlobalTest(learnsample, weights, fitmem, varctrl, gtctrl, 0, 
                 REAL(teststat), REAL(criterion), 1);
                 
    PutRNGstate();
    
    UNPROTECT(1);
    return(ans);
}
