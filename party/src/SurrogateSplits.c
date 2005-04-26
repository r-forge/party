
/**
    Suggorgate splits
    *\file $RCSfile$
    *\author $Author$
    *\date $Date$
*/
                
#include "PL2_common.h"

SEXP R_surrogates(SEXP node, SEXP learnsample, SEXP weights, SEXP controls, SEXP fitmem) {

    SEXP x, y, ym, expcovinf; 
    SEXP splitctrl, inputs; 
    SEXP ans, cutpoint, maxstat, splitstat, order; 
    SEXP variableID, tcutpoint, tmp, split;
    int nobs, ninputs, j, jselect, maxsurr;
    double ms;
    
    nobs = get_nobs(learnsample);
    ninputs = get_ninputs(learnsample);
    splitctrl = get_splitctrl(controls);
    maxsurr = get_maxsurrogate(splitctrl);
    if (maxsurr == 0) return(R_NilValue);

    ans = S3get_surrogatesplits(node);
    if (maxsurr != LENGTH(ans))
        error("nodes does not have %s surrogate splits", maxsurr);

    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    jselect = S3get_variableID(S3get_primarysplit(node));
    y = S3get_nodeweights(VECTOR_ELT(node, 7));
    PROTECT(ym = allocMatrix(REALSXP, nobs, 1));
    for (j = 0; j < nobs; j++) REAL(ym)[j] = REAL(y)[j];
    PROTECT(expcovinf = R_ExpectCovarInfluence(ym, weights));
    PROTECT(tmp = allocVector(VECSXP, ninputs));
    PROTECT(maxstat = allocVector(REALSXP, ninputs));
    PROTECT(order = allocVector(INTSXP, ninputs));
    PROTECT(splitstat = allocVector(REALSXP, nobs));    
    
    for (j = 0; j < ninputs; j++) {
    
         INTEGER(order)[j] = j + 1;
         REAL(maxstat)[j] = 0.0;
         SET_VECTOR_ELT(tmp, j, cutpoint = allocVector(REALSXP, 1));

         if ((j + 1) == jselect) continue;

         x = get_variable(inputs, j + 1);
         C_split(REAL(x), 1, REAL(ym), 1, REAL(weights), nobs,
                 INTEGER(get_ordering(inputs, j + 1)), 0, 0, splitctrl,
                 GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
                 expcovinf, REAL(cutpoint), &ms, REAL(splitstat));

         REAL(maxstat)[j] = -ms;
    }

    rsort_with_index(REAL(maxstat), INTEGER(order), ninputs);
    
    for (j = 0; j < maxsurr; j++) {
        SET_VECTOR_ELT(ans, j, split = allocVector(VECSXP, 2));
        SET_VECTOR_ELT(split, 0, variableID = allocVector(INTSXP, 1));
        SET_VECTOR_ELT(split, 1, tcutpoint = allocVector(REALSXP, 1));
        INTEGER(variableID)[0] = INTEGER(order)[j];
        REAL(tcutpoint)[0] = REAL(VECTOR_ELT(tmp, INTEGER(order)[j] - 1))[0];
    }

    UNPROTECT(6);
    return(ans);
}
 