                                                                   
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "party.h"

SEXP kronecker (SEXP A, SEXP B) {

    /* 
    *
    *    Kronecker Product of the 
    *    real (m x n) matrix A and the 
    *    real (r x s) matrix B
    *
    */ 

    /*  The Kronecker product, a real (mr x ns) matrix */
    
    SEXP ans; 

    int i, m, j, n, k, r, l, s, mr, ns, x = 0;
    double y, z;

    if (!isReal(A) || !isReal(B)) 
        error("A and B are not of type REAL");

    if (isMatrix(A)) {
        m = INTEGER(getAttrib(A, R_DimSymbol))[0];
        n = INTEGER(getAttrib(A, R_DimSymbol))[1];
    } else {
        m = LENGTH(A);
        n = 1;
    }
    
    if (isMatrix(B)) {
        r = INTEGER(getAttrib(B, R_DimSymbol))[0];
        s = INTEGER(getAttrib(B, R_DimSymbol))[1];
    } else {
        r = LENGTH(B);
        s = 1;
    }

    mr = m*r;
    ns = n*s;
  
    PROTECT(ans = allocMatrix(REALSXP, mr, ns));

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
      
            y = REAL(A)[aindx(i, j, m)];

            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++) {
                
                    x = aindx(i*r + k, j*s + l, mr);
                    z = REAL(B)[aindx(k, l, r)];
                    REAL(ans)[x] = y * z;
                    
                }
            }
        }
    }
    UNPROTECT(1);
    return(ans);
}  

void setAllZero(SEXP A) {

    int n, m, i, j;
    SEXP tmp;
    
    n = LENGTH(A);
    for (i = 0; i < n; i++) {
        tmp = VECTOR_ELT(A, i);
        m = LENGTH(tmp);
        for (j = 0; j < m; j++) REAL(tmp)[j] = 0.0;
    }
}

SEXP ec(SEXP Weights, SEXP Scores, SEXP cweights) {

    /*
    *
    *   Conditional Expectation and Covariance of 
    *   Linear Statistics of the form
    *
    *       L = vec(W %*% diag(cw) %*% S)
    *
    */
    
    /*  (p x nobs) matrix of weights  */

    SEXP W;		
    
    /*  (nobs x q) matrix of scores   */

    SEXP S;
    
    /*  nobs vector of case weights   */

    SEXP cw;
    
    /*  list of two return values: conditional expectation and covariance  */

    SEXP ans, expL, covL;
    
    /*  dimensions of W and S and corresponding loop variables  */

    int nobs, i;
    int p, k;
    int q, j; 
    int pq;
    
    /*  sum of case weights  */

    double scw;
    
    /*  mothers little helpers  */

    double f1, f2;
    SEXP ES, VS, VT1, VT2;
    SEXP wi, wiT, swi, swiT, wi_k_VS, VTp, wi_k_wiT;
    SEXP helpers;
    
    /* coerce the inputs to REALSXPs */
    
    PROTECT(W  = coerceVector(Weights, REALSXP));
    PROTECT(S  = coerceVector(Scores, REALSXP));
    PROTECT(cw = coerceVector(cweights, REALSXP));

    /* determine the dimensions and some checks */

    nobs = INTEGER(getAttrib(W, R_DimSymbol))[1];
    p    = INTEGER(getAttrib(W, R_DimSymbol))[0];
    q    = INTEGER(getAttrib(S, R_DimSymbol))[1];
    pq   = p * q;
    
    if (INTEGER(getAttrib(S, R_DimSymbol))[0] != nobs)
        error("score matrix does not have %d rows", nobs);
    if (LENGTH(cw) != nobs) 
        error("vector of case weights does not have %d elements", nobs);

    /*  compute the sum of the case weights */
        
    scw = 0;
    for (i = 0; i < nobs; i++) scw = scw + REAL(cw)[i];

    /*  allocate storage: the list of return values */

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, expL = allocVector(REALSXP, pq));
    SET_VECTOR_ELT(ans, 1, covL = allocMatrix(REALSXP, pq, pq));

    /* allocate storage: all helpers */
    
    PROTECT(helpers = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(helpers, 0, ES = allocMatrix(REALSXP, 1, q));
    SET_VECTOR_ELT(helpers, 1, VS = allocMatrix(REALSXP, q, q));
    SET_VECTOR_ELT(helpers, 2, wi = allocMatrix(REALSXP, p, 1));
    SET_VECTOR_ELT(helpers, 3, wiT = allocMatrix(REALSXP, 1, p));
    SET_VECTOR_ELT(helpers, 4, swi = allocMatrix(REALSXP, p, 1));
    SET_VECTOR_ELT(helpers, 5, swiT = allocMatrix(REALSXP, 1, p));
    SET_VECTOR_ELT(helpers, 6, VTp = allocMatrix(REALSXP, p, p));

    /* make sure the helpers are initially zero */

    setAllZero(helpers);

    /* 
    *   ES:  the expectation of the scores 
    *   swi: row sums of the weights
    */
    

    for (i = 0; i < nobs; i++) {

        /*  observations with zero case weights do not contribute */
    
        if (REAL(cw)[i] == 0.0) continue;
    
        for (j = 0; j < q; j++) {
            REAL(ES)[j] = REAL(ES)[j] 
                          + REAL(cw)[i] * REAL(S)[aindx(i, j, nobs)];
        } 
        for (k = 0; k < p; k++) {
            REAL(swi)[k] = REAL(swi)[k] 
                           + REAL(cw)[i] * REAL(W)[aindx(k, i, p)];
        }
    }

    for (j = 0; j < q; j++) {
        REAL(ES)[j] = REAL(ES)[j] / scw;
    }

    /*
    *   expL: expectation of the linear statistic L
    *   swiT: transpose of the row sum of the weights
    */

    for (k = 0; k < p; k++) {
        REAL(swiT)[k] = REAL(swi)[k];
        for (j = 0; j < q; j++) {
            REAL(expL)[aindx(k,j,p)] = REAL(swi)[k] * REAL(ES)[j];
        }
    }

    /*
    *   VS:  covariance of the scores
    */ 

    for (i = 0; i < nobs; i++) {

        if (REAL(cw)[i] == 0.0) continue;
     
        for (j = 0; j < q; j++) {
            for (k = 0; k < q; k++) {
                REAL(VS)[aindx(k, j, q)] = REAL(VS)[aindx(k, j, q)] + 
                    REAL(cw)[i] * (REAL(S)[aindx(i, k, nobs)] - REAL(ES)[k]) * 
                                  (REAL(S)[aindx(i, j, nobs)] - REAL(ES)[j]);
            }
        }
    }

    for (j = 0; j < q*q; j++) {
        REAL(VS)[j] = REAL(VS)[j] / scw;
    }
    
    /* 
    *   covL:  covariance of the linear statistic L
    */
    
    for (i = 0; i < nobs; i++) {

        if (REAL(cw)[i] == 0.0) continue;
        
        for (k = 0; k < p; k++) {
            REAL(wi)[k] = REAL(W)[aindx(k, i, p)];
            REAL(wiT)[k] = REAL(W)[aindx(k, i, p)];
        } 

        wi_k_wiT = kronecker(wi, wiT);

        for (k = 0; k < p*p; k++) {
            REAL(VTp)[k] = REAL(VTp)[k] + REAL(cw)[i] * REAL(wi_k_wiT)[k];
        }
    }
    
    VT1 = kronecker(VS, VTp);

    wi_k_VS = kronecker(VS, swi);
    VT2 = kronecker(wi_k_VS, swiT);

    f1 = scw/(scw - 1);
    f2 = (1/(scw - 1));

    for (k = 0; k < (pq * pq); k++) {
        REAL(covL)[k] = f1 * REAL(VT1)[k] - f2 * REAL(VT2)[k];
    }

    UNPROTECT(5);
    return(ans);
}


SEXP evS(SEXP Scores, SEXP cweights) {

    /*
    *    expectation and variance (!) of the scores only
    */


    SEXP S;
    
    /*  nobs vector of case weights   */

    SEXP cw;
    
    /*  list of two return values: conditional expectation and variance  */

    SEXP ans, ES, VS, scw;
    
    /*  dimensions of W and S and corresponding loop variables  */

    int nobs, i;
    int k;
    int q, j; 
    
    /* coerce the inputs to REALSXPs */
    
    PROTECT(S  = coerceVector(Scores, REALSXP));
    PROTECT(cw = coerceVector(cweights, REALSXP));

    /* determine the dimensions and some checks */

    nobs = INTEGER(getAttrib(S, R_DimSymbol))[0];
    q    = INTEGER(getAttrib(S, R_DimSymbol))[1];
    
    if (LENGTH(cw) != nobs) 
        error("vector of case weights does not have %d elements", nobs);

    /*  allocate storage: the list of return values */

    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, ES = allocVector(REALSXP, q));
    SET_VECTOR_ELT(ans, 1, VS = allocVector(REALSXP, q));
    SET_VECTOR_ELT(ans, 2, scw = allocVector(REALSXP, 1));
    
    setAllZero(ans);

    for (i = 0; i < nobs; i++) {

        if (REAL(cw)[i] == 0.0) continue;

        REAL(scw)[0] = REAL(scw)[0] + REAL(cw)[i];

        for (k = 0; k < q; k++) {
            REAL(ES)[k] = REAL(ES)[k] 
                          + REAL(cw)[i] * REAL(S)[aindx(i, k, nobs)];
        }

    }

    for (k = 0; k < q; k++) {
        REAL(ES)[k] = REAL(ES)[k] / REAL(scw)[0];
    }

    for (i = 0; i < nobs; i++) {
        
        if (REAL(cw)[i] == 0.0) continue;
    
        for (j = 0; j < q; j++) {
                REAL(VS)[j] = REAL(VS)[j] + 
                    REAL(cw)[i] * (REAL(S)[aindx(i, j, nobs)] - REAL(ES)[j]) * 
                                  (REAL(S)[aindx(i, j, nobs)] - REAL(ES)[j]);
        }
    }
    
    for (k = 0; k < q; k++) {
        REAL(VS)[k] = REAL(VS)[k] / REAL(scw)[0];
    }

    UNPROTECT(3);
    return(ans);
}


SEXP evL(SEXP Weights, SEXP Scores, SEXP cweights, SEXP evSans) {

    /*
     *
     *   Conditional Expectation and Variance of
     *   Linear Statistics of the form
     *
     *       L = vec(W %*% diag(cw) %*% S)
     *
     */
                                
    /*  (p x nobs) matrix of weights  */

    SEXP W;	
    
    /*  (nobs x q) matrix of scores   */

    SEXP S;
    
    /*  nobs vector of case weights   */

    SEXP cw;
    
    
    /*  list of two return values: conditional expectation and variance  */

    SEXP ans, expL, varL;
    
    /*  dimensions of W and S and corresponding loop variables  */

    int nobs, i;
    int p, k;
    int q, j; 
    int pq;
    
    double scw;
    
    /*  mothers little helpers  */

    double f1, f2;
    SEXP helpers, ES, VS, wi, wii;
    
    /* coerce the inputs to REALSXPs */
    
    PROTECT(W  = coerceVector(Weights, REALSXP));
    PROTECT(S  = coerceVector(Scores, REALSXP));
    PROTECT(cw = coerceVector(cweights, REALSXP));

    /* determine the dimensions and some checks */

    nobs = INTEGER(getAttrib(W, R_DimSymbol))[1];
    p    = INTEGER(getAttrib(W, R_DimSymbol))[0];
    q    = INTEGER(getAttrib(S, R_DimSymbol))[1];
    pq   = p * q;
    
    if (INTEGER(getAttrib(S, R_DimSymbol))[0] != nobs)
        error("score matrix does not have %d rows", nobs);
    if (LENGTH(cw) != nobs) 
        error("vector of case weights does not have %d elements", nobs);

    /*  compute the sum of the case weights */
        
    scw = REAL(VECTOR_ELT(evSans, 2))[0];
    ES = VECTOR_ELT(evSans, 0);
    VS = VECTOR_ELT(evSans, 1);
    
    /*  allocate storage: the list of return values */

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, expL = allocVector(REALSXP, pq));
    SET_VECTOR_ELT(ans, 1, varL = allocVector(REALSXP, pq));
    
    PROTECT(helpers = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(helpers, 0, wi = allocVector(REALSXP, p));
    SET_VECTOR_ELT(helpers, 1, wii = allocVector(REALSXP, p));

    setAllZero(helpers);

    for (i = 0; i < nobs; i++) {

        if (REAL(cw)[i] == 0.0) continue;

        for (k = 0; k < p; k++) {
            REAL(wi)[k] = REAL(wi)[k] 
                          + REAL(cw)[i] * REAL(W)[aindx(k, i, p)];
        }

        for (k = 0; k < p; k++) {
            REAL(wii)[k] = REAL(wii)[k] 
                         + REAL(cw)[i] * REAL(W)[aindx(k, i, p)] 
                           * REAL(W)[aindx(k, i, p)];
        }
    }

    f1 = scw/(scw - 1);
    f2 = (1/(scw - 1));
    for (k = 0; k < p; k++) {
        for (j = 0; j < q; j++) {
            REAL(expL)[j*p + k] = REAL(ES)[j] * REAL(wi)[k];
            REAL(varL)[j*p + k] = f1*REAL(VS)[j] * REAL(wii)[k] 
                                - f2*REAL(VS)[j] * REAL(wi)[k]*REAL(wi)[k];
        }
    }
    UNPROTECT(5);
    return(ans);
}

SEXP ev(SEXP Weights, SEXP Scores, SEXP cweights) {

    /*
     *
     *   Conditional Expectation and Variance of
     *   Linear Statistics of the form
     *
     *       L = vec(W %*% diag(cw) %*% S)
     *
     */
                                

    SEXP evSans, ans;
    
    evSans = evS(Scores, cweights);
    ans = evL(Weights, Scores, cweights, evSans);
    return(ans);

}
