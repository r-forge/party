
/**
    The tree growing recursion
    *\file $RCSfile$
    *\author $Author$
    *\date $Date$
*/

#include "PL2_common.h"

void C_TreeGrow(SEXP node, SEXP learnsample, SEXP fitmem, 
                SEXP controls, int *where, int *nodenum) {

    SEXP weights;
    int nobs, i;
    double *dweights;
    
    weights = S3get_nodeweights(node);
    
    if ((nodenum[0] == 2 || nodenum[0] == 3) && 
        get_stump(GET_SLOT(controls, PL2_tgctrlSym)))
        C_Node(node, learnsample, weights, fitmem, controls, 1);
    else
        C_Node(node, learnsample, weights, fitmem, controls, 0);
    
    S3set_nodeID(node, nodenum[0]);    
    
    if (!S3get_nodeterminal(node)) {

        C_splitnode(node, learnsample, controls);

        nodenum[0] += 1;
        C_TreeGrow(S3get_leftnode(node), learnsample, fitmem, 
                   controls, where, nodenum);

        nodenum[0] += 1;                                      
        C_TreeGrow(S3get_rightnode(node), learnsample, fitmem, 
                   controls, where, nodenum);
    } else {
        dweights = REAL(weights);
        nobs = get_nobs(learnsample);
        for (i = 0; i < nobs; i++)
            if (dweights[i] > 0) where[i] = nodenum[0];
    } 
}

SEXP R_TreeGrow(SEXP learnsample, SEXP weights, SEXP fitmem, SEXP controls, SEXP where) {
            
     SEXP ans, nweights;
     double *dnweights, *dweights;
     int nobs, i, nodenum = 1;
     
     nobs = get_nobs(learnsample);
     PROTECT(ans = allocVector(VECSXP, 9));
     C_init_node(ans, nobs, get_ninputs(learnsample), get_maxsurrogate(get_splitctrl(controls)),
                 ncol(GET_SLOT(GET_SLOT(learnsample, PL2_responsesSym), 
                      PL2_jointtransfSym)));

     nweights = S3get_nodeweights(ans);
     dnweights = REAL(nweights);
     dweights = REAL(weights);
     for (i = 0; i < nobs; i++) dnweights[i] = dweights[i];
     
     C_TreeGrow(ans, learnsample, fitmem, controls, INTEGER(where), &nodenum);
     UNPROTECT(1);
     return(ans);
}


SEXP R_Ensemble2(SEXP learnsample, SEXP weights, SEXP fitmem, SEXP controls, SEXP ans) {
            
     SEXP nweights, tree, where;
     double *dnweights, *dweights, sw = 0.0, *prob;
     int nobs, i, b, B , nodenum = 1, *iweights, *iwhere;
     
     B = LENGTH(ans);
     nobs = get_nobs(learnsample);

     iweights = Calloc(nobs, int);
     prob = Calloc(nobs, double);
     dweights = REAL(weights);

     for (i = 0; i < nobs; i++)
         sw += dweights[i];
     for (i = 0; i < nobs; i++)
         prob[i] = dweights[i]/sw;

     for (b  = 0; b < B; b++) {
         SET_VECTOR_ELT(ans, b, tree = allocVector(VECSXP, 10));
         SET_VECTOR_ELT(tree, 9, where = allocVector(INTSXP, nobs));
         iwhere = INTEGER(where);
         for (i = 0; i < nobs; i++) iwhere[i] = 0;
     
         C_init_node(tree, nobs, get_ninputs(learnsample), 
                     get_maxsurrogate(get_splitctrl(controls)),
                     ncol(GET_SLOT(GET_SLOT(learnsample, PL2_responsesSym), 
                          PL2_jointtransfSym)));

         /* weights for a bootstrap sample */
         GetRNGstate();
         rmultinom((int) sw, prob, nobs, iweights);
         PutRNGstate();

         nweights = S3get_nodeweights(tree);
         dnweights = REAL(nweights);
         for (i = 0; i < nobs; i++) dnweights[i] = (double) iweights[i];
     
         C_TreeGrow(tree, learnsample, fitmem, controls, iwhere, &nodenum);
         nodenum = 1;
     }
     Free(prob); Free(iweights);
     return(ans);
}

