
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

int aindx(int i, int j, int m) {
  /*
      m x n matrix. (i+1,j+1)
  */
  return(j * m + i);
}
            