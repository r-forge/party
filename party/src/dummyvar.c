                                                                   
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "party.h"

SEXP dummyvar(SEXP x, SEXP transposed) {

   int nlevels, n, i, j;
   
   SEXP ans;
   SEXP levels;
   SEXP dimnames;
   
   if (!isFactor(x)) error("x is not a factor\n");
   if (!isLogical(transposed)) error("transposed is not a logical\n");
   
   n = LENGTH(x);
   levels = GET_LEVELS(x);
   nlevels = LENGTH(levels);

   if (LOGICAL(transposed)[0]) { /* == TRUE */
       PROTECT(ans = allocMatrix(INTSXP, nlevels, n));
   
       for (i = 0; i < n; i++) {
           for (j = 0; j < nlevels; j++) {
               if (INTEGER(x)[i] == (j+1)) 
                   INTEGER(ans)[aindx(j, i, nlevels)] = 1;
               else
                   INTEGER(ans)[aindx(j, i, nlevels)] = 0;
           }
       }
       PROTECT(dimnames = allocVector(VECSXP, 2));
       SET_VECTOR_ELT(dimnames, 0, levels);
       SET_VECTOR_ELT(dimnames, 1, R_NilValue);
       setAttrib(ans, R_DimNamesSymbol, dimnames);
   } else {
       PROTECT(ans = allocMatrix(INTSXP, n, nlevels));
   
       for (i = 0; i < n; i++) {
           for (j = 0; j < nlevels; j++) {
               if (INTEGER(x)[i] == (j+1)) 
                   INTEGER(ans)[aindx(i, j, n)] = 1;
               else
                   INTEGER(ans)[aindx(i, j, n)] = 0;
           } 
       }
    
       PROTECT(dimnames = allocVector(VECSXP, 2));
       SET_VECTOR_ELT(dimnames, 0, R_NilValue);
       SET_VECTOR_ELT(dimnames, 1, levels);
       setAttrib(ans, R_DimNamesSymbol, dimnames);
   }
   UNPROTECT(2);
   return(ans);
}

SEXP contcat (SEXP x, SEXP cutpoints) {

    int i, j, nx, nc;
    SEXP ans, rx, rc;
    
    PROTECT(rx = coerceVector(x, REALSXP));
    PROTECT(rc = coerceVector(cutpoints, REALSXP));
    
    nx = LENGTH(rx);
    nc = LENGTH(rc);
    
    PROTECT(ans = allocMatrix(INTSXP, nc, nx));
    
    for (i = 0; i < nx; i++) {
        for (j = 0; j < nc; j++) {
            if (REAL(rx)[i] <= REAL(rc)[j]) {
                INTEGER(ans)[aindx(j, i, nc)] = 1;
            } else {
                INTEGER(ans)[aindx(j, i, nc)] = 0;
            }
        }
    }
    
    UNPROTECT(3);
    return(ans);
}

