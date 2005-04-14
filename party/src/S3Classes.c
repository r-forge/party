
/**
    S3 classes for dealing with nodes and splits
    *\file $RCSfile$
    *\author $Author$
    *\date $Date$
*/
                
#include "PL2_common.h"
                
void C_init_node(SEXP node, int nobs, int ninputs, int nsurr, int q) {

    SEXP nodeID, weights, criterion, primarysplit, surrogatesplits, 
         terminal, prediction;

    SET_VECTOR_ELT(node, 0, nodeID = allocVector(INTSXP, 1));
    if (nobs > 0) 
        SET_VECTOR_ELT(node, 1, weights = allocVector(REALSXP, nobs));
    else
        SET_VECTOR_ELT(node, 1, R_NilValue);
    SET_VECTOR_ELT(node, 2, criterion = allocVector(VECSXP, 3));
    /* teststats */
    SET_VECTOR_ELT(criterion, 0, allocVector(REALSXP, ninputs)); 
    /* criterion, aka pvalues */
    SET_VECTOR_ELT(criterion, 1, allocVector(REALSXP, ninputs));
    /* max(criterion) */
    SET_VECTOR_ELT(criterion, 2, allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(node, 3, terminal = allocVector(LGLSXP, 1));
    INTEGER(terminal)[0] = 0;
    SET_VECTOR_ELT(node, 4, primarysplit = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(node, 5, surrogatesplits = allocVector(VECSXP, nsurr));
    SET_VECTOR_ELT(node, 6, prediction = allocVector(REALSXP, q));

}

void S3set_nodeID(SEXP node, int nodeID) {
    INTEGER(VECTOR_ELT(node, 0))[0] = nodeID;
}

int S3get_nodeID(SEXP node) {
    return(INTEGER(VECTOR_ELT(node, 0))[0]);
}

SEXP S3get_nodeweights(SEXP node) {
    SEXP ans;
    
    ans = VECTOR_ELT(node, 1);
    if (ans == R_NilValue)
        error("node has no weights element"); 
    return(VECTOR_ELT(node, 1));
}

SEXP S3get_teststat(SEXP node) {
    return(VECTOR_ELT(VECTOR_ELT(node, 2), 0));
}

SEXP S3get_criterion(SEXP node) {
    return(VECTOR_ELT(VECTOR_ELT(node, 2), 1));
}

SEXP S3get_maxcriterion(SEXP node) {
    return(VECTOR_ELT(VECTOR_ELT(node, 2), 2));
}

void S3set_nodeterminal(SEXP node) {
    INTEGER(VECTOR_ELT(node, 3))[0] = 1;
}

int S3get_nodeterminal(SEXP node) {
    return(INTEGER(VECTOR_ELT(node, 3))[0]);
}

SEXP S3get_primarysplit(SEXP node) {
    return(VECTOR_ELT(node, 4));
}

SEXP S3get_surrogatesplits(SEXP node) {
    return(VECTOR_ELT(node, 5));
}

SEXP S3get_prediction(SEXP node) {
    return(VECTOR_ELT(node, 6));
}

SEXP S3get_leftnode(SEXP node) {
    return(VECTOR_ELT(node, 7));
}

SEXP S3get_rightnode(SEXP node) {
    return(VECTOR_ELT(node, 8));
}

void C_init_orderedsplit(SEXP split, int nobs) {
    
    SEXP variableID, splitpoint, splitstatistics, ordered;
    
    SET_VECTOR_ELT(split, 0, variableID = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(split, 1, ordered = allocVector(LGLSXP, 1));
    INTEGER(ordered)[0] = 1;
    SET_VECTOR_ELT(split, 2, splitpoint = allocVector(REALSXP, 1));
    if (nobs > 0)
        SET_VECTOR_ELT(split, 3, 
                       splitstatistics = allocVector(REALSXP, nobs));
    else
        SET_VECTOR_ELT(split, 3, R_NilValue);
}

void C_init_nominalsplit(SEXP split, int nlevels, int nobs) {
    
    SEXP variableID, splitpoint, splitstatistics, ordered;
    
    SET_VECTOR_ELT(split, 0, variableID = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(split, 1, ordered = allocVector(LGLSXP, 1));
    INTEGER(ordered)[0] = 0;
    SET_VECTOR_ELT(split, 2, splitpoint = allocVector(INTSXP, nlevels));
    if (nobs > 0)
        SET_VECTOR_ELT(split, 3, 
                       splitstatistics = allocVector(REALSXP, nobs));
    else
        SET_VECTOR_ELT(split, 3, R_NilValue);
}

void S3set_variableID(SEXP split, int variableID) {
    INTEGER(VECTOR_ELT(split, 0))[0] = variableID;
}

int S3get_variableID(SEXP split) {
    return(INTEGER(VECTOR_ELT(split, 0))[0]);
}

int S3is_ordered(SEXP split) {
    return(INTEGER(VECTOR_ELT(split, 1))[0]);
}

void S3set_ordered(SEXP split) {
    INTEGER(VECTOR_ELT(split, 1))[0] = 1;
}

void S3set_nominal(SEXP split) {
    INTEGER(VECTOR_ELT(split, 1))[0] = 0;
}

SEXP S3get_splitpoint(SEXP split) {
   return(VECTOR_ELT(split, 2));
}
   
SEXP S3get_splitstatistics(SEXP split) {
   SEXP ans;
   
   ans = VECTOR_ELT(split, 3);
   if (ans == R_NilValue)
       error("split has not splitstatistics element");
   return(ans);
}
