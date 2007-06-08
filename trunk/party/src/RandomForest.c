
/**
    Random forest with conditional inference trees
    *\file RandomForest.c
    *\author $Author$
    *\date $Date$
*/

#include "party.h"


/**
    An experimental implementation of random forest like algorithms \n
    *\param learnsample an object of class `LearningSample'
    *\param weights a vector of case weights
    *\param fitmem an object of class `TreeFitMemory'
    *\param controls an object of class `TreeControl'
*/


SEXP R_Ensemble(SEXP learnsample, SEXP weights, SEXP fitmem, SEXP controls) {
            
     SEXP nweights, tree, where, ans;
     double *dnweights, *dweights, sw = 0.0, *prob;
     int nobs, i, b, B , nodenum = 1, *iweights, *iweightstmp, 
         *iwhere, replace, fraction, wgrzero = 0;
     
     B = get_ntree(controls);
     nobs = get_nobs(learnsample);
     
     PROTECT(ans = allocVector(VECSXP, B));

     iweights = Calloc(nobs, int);
     iweightstmp = Calloc(nobs, int);
     prob = Calloc(nobs, double);
     dweights = REAL(weights);

     for (i = 0; i < nobs; i++) {
         sw += dweights[i];
         if (dweights[i] > 0) wgrzero++;
     }
     for (i = 0; i < nobs; i++)
         prob[i] = dweights[i]/sw;

     replace = get_replace(controls);
     /* fraction of number of obs with weight > 0 */
     fraction = (int) (get_fraction(controls) * wgrzero);

     if (!replace) {
         if (fraction < 10)
             error("fraction of %f is too small", fraction);
     }

     /* <FIXME> can we call those guys ONCE? what about the deeper
         calls??? </FIXME> */
     GetRNGstate();
  
     for (b  = 0; b < B; b++) {
         SET_VECTOR_ELT(ans, b, tree = allocVector(VECSXP, NODE_LENGTH + 1));
         SET_VECTOR_ELT(tree, NODE_LENGTH, where = allocVector(INTSXP, nobs));
         iwhere = INTEGER(where);
         for (i = 0; i < nobs; i++) iwhere[i] = 0;
     
         C_init_node(tree, nobs, get_ninputs(learnsample), 
                     get_maxsurrogate(get_splitctrl(controls)),
                     ncol(get_predict_trafo(GET_SLOT(learnsample, 
                                                   PL2_responsesSym))));

         /* generate altered weights for perturbation */
         if (replace) {
             /* weights for a bootstrap sample */
             rmultinom((int) sw, prob, nobs, iweights);
         } else {
             /* weights for sample splitting */
             C_SampleSplitting(nobs, prob, iweights, fraction);
         }

         nweights = S3get_nodeweights(tree);
         dnweights = REAL(nweights);
         for (i = 0; i < nobs; i++) dnweights[i] = (double) iweights[i];
     
         C_TreeGrow(tree, learnsample, fitmem, controls, iwhere, &nodenum, 1);
         nodenum = 1;
     }

     PutRNGstate();

     Free(prob); Free(iweights); Free(iweightstmp);
     UNPROTECT(1);
     return(ans);
}
