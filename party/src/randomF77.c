
/* $Id$
*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include "party.h"

void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(unifrnd)(void) { return unif_rand(); }
