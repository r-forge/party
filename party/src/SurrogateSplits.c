
/**
    Suggorgate splits
    *\file $RCSfile$
    *\author $Author$
    *\date $Date$
*/
                
#include "PL2_common.h"

SEXP R_surrogates(SEXP node, SEXP learnsample, SEXP weights, SEXP controls, SEXP fitmem) {

    SEXP x, y, expcovinf; 
    SEXP splitctrl, inputs; 
    SEXP ans;
    SEXP split, thiswhichNA;
    int nobs, ninputs, i, j, k, jselect, maxsurr, *order;
    double ms, cp, *thisweights, *cutpoint, *maxstat, *splitstat;
    
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

    expcovinf = GET_SLOT(fitmem, PL2_expcovinfssSym);
    C_ExpectCovarInfluence(REAL(y), 1, REAL(weights), nobs, expcovinf);
    
    maxstat = Calloc(ninputs, double);
    cutpoint = Calloc(ninputs, double);
    splitstat = Calloc(nobs, double);
    order = Calloc(ninputs, int);
    
    for (j = 0; j < ninputs; j++) {
    
         order[j] = j + 1;
         maxstat[j] = 0.0;
         cutpoint[j] = 0.0;

         /* ordered input variables only (for the moment) */
         if ((j + 1) == jselect || is_nominal(inputs, j + 1))
             continue;

         x = get_variable(inputs, j + 1);

         if (has_missings(inputs, j + 1)) {

             thisweights = REAL(get_weights(fitmem, j + 1));
             for (i = 0; i < nobs; i++) thisweights[i] = REAL(weights)[i];
             thiswhichNA = get_missings(inputs, j + 1);
             for (k = 0; k < LENGTH(thiswhichNA); k++)
                 thisweights[INTEGER(thiswhichNA)[k] - 1] = 0.0;
                 
             C_ExpectCovarInfluence(REAL(y), 1, thisweights, nobs, expcovinf);
             
             C_split(REAL(x), 1, REAL(y), 1, thisweights, nobs,
                     INTEGER(get_ordering(inputs, j + 1)), 0, 0, splitctrl,
                     GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
                     expcovinf, &cp, &ms, splitstat);
         } else {
         
             C_split(REAL(x), 1, REAL(y), 1, REAL(weights), nobs,
             INTEGER(get_ordering(inputs, j + 1)), 0, 0, splitctrl,
             GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
             expcovinf, &cp, &ms, splitstat);
         }

         maxstat[j] = -ms;
         cutpoint[j] = cp;
    }

    rsort_with_index(maxstat, order, ninputs);
    
    for (j = 0; j < maxsurr; j++) {
        Rprintf("var: %d %f\n", order[j], cutpoint[order[j] - 1]);
        SET_VECTOR_ELT(S3get_surrogatesplits(node), j, split = allocVector(VECSXP, 4));
        C_init_orderedsplit(split, 0);
        S3set_variableID(split, order[j]);
        REAL(S3get_splitpoint(split))[0] = cutpoint[order[j] - 1];
    }
    
    Free(maxstat);
    Free(cutpoint);
    Free(splitstat);
    Free(order);

    return(ans);
}
