
/**
    Suggorgate splits
    *\file $RCSfile$
    *\author $Author$
    *\date $Date$
*/
                
#include "PL2_common.h"

void C_surrogates(SEXP node, SEXP learnsample, SEXP weights, SEXP controls, SEXP fitmem) {

    SEXP x, y, expcovinf; 
    SEXP splitctrl, inputs; 
    SEXP ans;
    SEXP split, thiswhichNA;
    int nobs, ninputs, i, j, k, jselect, maxsurr, *order;
    double ms, cp, *thisweights, *cutpoint, *maxstat, *splitstat, *dweights, *tweights, *dx, *dy;
    double cut, *twotab;
    
    nobs = get_nobs(learnsample);
    ninputs = get_ninputs(learnsample);
    splitctrl = get_splitctrl(controls);
    maxsurr = get_maxsurrogate(splitctrl);

    ans = S3get_surrogatesplits(node);
    if (maxsurr != LENGTH(ans))
        error("nodes does not have %s surrogate splits", maxsurr);

    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    jselect = S3get_variableID(S3get_primarysplit(node));
    y = S3get_nodeweights(VECTOR_ELT(node, 7));

    tweights = Calloc(nobs, double);
    dweights = REAL(weights);
    for (i = 0; i < nobs; i++) tweights[i] = dweights[i];
    if (has_missings(inputs, jselect)) {
        thiswhichNA = get_missings(inputs, jselect);
        for (k = 0; k < LENGTH(thiswhichNA); k++)
            tweights[INTEGER(thiswhichNA)[k] - 1] = 0.0;
    }

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
             for (i = 0; i < nobs; i++) thisweights[i] = tweights[i];
             thiswhichNA = get_missings(inputs, j + 1);
             for (k = 0; k < LENGTH(thiswhichNA); k++)
                 thisweights[INTEGER(thiswhichNA)[k] - 1] = 0.0;
                 
             C_ExpectCovarInfluence(REAL(y), 1, thisweights, nobs, expcovinf);
             
             C_split(REAL(x), 1, REAL(y), 1, thisweights, nobs,
                     INTEGER(get_ordering(inputs, j + 1)), 0, 0, splitctrl,
                     GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
                     expcovinf, &cp, &ms, splitstat);
         } else {
         
             C_split(REAL(x), 1, REAL(y), 1, tweights, nobs,
             INTEGER(get_ordering(inputs, j + 1)), 0, 0, splitctrl,
             GET_SLOT(fitmem, PL2_linexpcov2sampleSym),
             expcovinf, &cp, &ms, splitstat);
         }

         maxstat[j] = -ms;
         cutpoint[j] = cp;
    }

    rsort_with_index(maxstat, order, ninputs);
    
    twotab = Calloc(4, double);
    
    for (j = 0; j < maxsurr; j++) {

        for (i = 0; i < 4; i++) twotab[i] = 0.0;
        cut = cutpoint[order[j] - 1];
        SET_VECTOR_ELT(S3get_surrogatesplits(node), j, split = allocVector(VECSXP, 5));
        C_init_orderedsplit(split, 0);
        S3set_variableID(split, order[j]);
        REAL(S3get_splitpoint(split))[0] = cut;
        dx = REAL(get_variable(inputs, order[j]));
        dy = REAL(y);

        for (i = 0; i < nobs; i++) {
            twotab[0] += ((dy[i] == 1) && (dx[i] <= cut)) * tweights[i];
            twotab[1] += (dy[i] == 1) * tweights[i];
            twotab[2] += (dx[i] <= cut) * tweights[i];
            twotab[3] += tweights[i];
        }
        S3set_toleft(split, (int) (twotab[0] - twotab[1] * twotab[2] / twotab[3]) > 0);
    }
    
    Free(maxstat);
    Free(cutpoint);
    Free(splitstat);
    Free(order);
    Free(tweights);
    Free(twotab);
}

SEXP R_surrogates(SEXP node, SEXP learnsample, SEXP weights, SEXP controls, SEXP fitmem) {

    C_surrogates(node, learnsample, weights, controls, fitmem);
    return(S3get_surrogatesplits(node));
    
}

void C_splitsurrogate(SEXP node, SEXP learnsample) {

    SEXP weights, split, surrsplit;
    SEXP inputs, whichNA;
    double cutpoint, *dx, *dweights, *leftweights, *rightweights;
    int *iwhichNA, k;
    int nobs, i, nna, ns;
                    
    weights = S3get_nodeweights(node);
    dweights = REAL(weights);
    inputs = GET_SLOT(learnsample, PL2_inputsSym);
    nobs = get_nobs(learnsample);
            
    leftweights = REAL(S3get_nodeweights(S3get_leftnode(node)));
    rightweights = REAL(S3get_nodeweights(S3get_rightnode(node)));
    surrsplit = S3get_surrogatesplits(node);

    /* if the primary split has any missings */
    split = S3get_primarysplit(node);
    if (has_missings(inputs, S3get_variableID(split))) {

        /* where are the missings? */
        whichNA = get_missings(inputs, S3get_variableID(split));
        iwhichNA = INTEGER(whichNA);
        nna = LENGTH(whichNA);

        /* for all missing values ... */
        for (k = 0; k < nna; k++) {
            ns = 0;
            i = iwhichNA[k] - 1;
            if (dweights[i] == 0) continue;
            
            /* loop over surrogate splits until an appropriate one is found */
            while(TRUE) {
            
                if (ns >= LENGTH(surrsplit)) break;
            
                split = VECTOR_ELT(surrsplit, ns);
                if (has_missings(inputs, S3get_variableID(split))) {
                    if (INTEGER(get_missings(inputs, S3get_variableID(split)))[i]) {
                        ns++;
                        continue;
                    }
                }

                cutpoint = REAL(S3get_splitpoint(split))[0];
                dx = REAL(get_variable(inputs, S3get_variableID(split)));

                if (S3get_toleft(split)) {
                    if (dx[i] <= cutpoint) {
                        leftweights[i] = dweights[i];
                        rightweights[i] = 0.0;
                    } else {
                        rightweights[i] = dweights[i];
                        leftweights[i] = 0.0;
                    }
                } else {
                    if (dx[i] <= cutpoint) {
                        rightweights[i] = dweights[i];
                        leftweights[i] = 0.0;
                    } else {
                        leftweights[i] = dweights[i];
                        rightweights[i] = 0.0;
                    }
                }
                break;
            }
        }
    }
}
