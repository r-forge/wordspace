/*
 *  Compute different row norms for dense or sparse matrix
 */

/* internal codes for norms:
 *  0 = euclidean
 *  1 = maximum
 *  2 = manhattan
 *  3 = minkowski
 */

#include "function_defs.h"

void row_norms_dense(double *norms, int *nr, int *nc, double *x, int *norm_code, double *p_norm) {
  int row, col, i;
  
  if (*norm_code < 0 || *norm_code > 3)
    error("norm #%d is not defined -- internal error", *norm_code);
  if (*norm_code == 3 && (*p_norm < 1 || !R_FINITE(*p_norm)))
    error("Minkowski p-norm is not defined for p = %g", *p_norm);

  for (row = 0; row < *nr; row++)
    norms[row] = 0; /* initialise norm accumulators */
  
  i = 0;
  for (col = 0; col < *nc; col++) for (row = 0; row < *nr; row++) {
    if      (*norm_code == 0) norms[row] += x[i] * x[i];
    else if (*norm_code == 1) { if (fabs(x[i]) > norms[row]) norms[row] = fabs(x[i]); }
    else if (*norm_code == 2) norms[row] += fabs(x[i]);
    else if (*norm_code == 3) norms[row] += pow(fabs(x[i]), *p_norm);
    i++;
  }
  
  if (*norm_code == 0) {
    for (row = 0; row < *nr; row++) 
      norms[row] = sqrt(norms[row]);
  }
  else if (*norm_code == 3) {
    for (row = 0; row < *nr; row++) 
      norms[row] = pow(norms[row], 1 / *p_norm);
  }
  else {
    /* no adjustment needed for Maximum and Manhattan norms */
  }
}

void row_norms_sparse(double *norms, int *nr, int *nc, int *p, int *row_of, double *x, int *norm_code, double *p_norm) {
  int row, col, i;
  
  if (*norm_code < 0 || *norm_code > 3)
    error("norm #%d is not defined -- internal error", *norm_code);
  if (*norm_code == 3 && (*p_norm < 1 || !R_FINITE(*p_norm)))
    error("Minkowski p-norm is not defined for p = %g", *p_norm);

  for (row = 0; row < *nr; row++)
    norms[row] = 0; /* initialise norm accumulators */
  
  for (col = 0; col < *nc; col++) {
    for (i = p[col]; i < p[col+1]; i++) {
      int row = row_of[i];
      if      (*norm_code == 0) norms[row] += x[i] * x[i];
      else if (*norm_code == 1) { if (fabs(x[i]) > norms[row]) norms[row] = fabs(x[i]); }
      else if (*norm_code == 2) norms[row] += fabs(x[i]);
      else if (*norm_code == 3) norms[row] += pow(fabs(x[i]), *p_norm);
    }
  }
  
  if (*norm_code == 0) {
    for (row = 0; row < *nr; row++) 
      norms[row] = sqrt(norms[row]);
  }
  else if (*norm_code == 3) {
    for (row = 0; row < *nr; row++) 
      norms[row] = pow(norms[row], 1 / *p_norm);
  }
  else {
    /* no adjustment needed for Maximum and Manhattan norms */
  }
}
