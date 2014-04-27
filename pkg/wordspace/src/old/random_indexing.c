/*
 *  Memory-friendly random indexing of a sparse matrix
 */

#include "function_defs.h"

/* Memory-friendly implementation of random indexing uses a very simple algorithm
 * that iterates through the dimensions of the original matrix M.  For each dimensions,
 * it generates a cross-section of the random basis vectors for this dimension, i.e.
 * a column of the random basis matrix Q, and updates the projected vectors with the
 * outer product of the column of Q and the corresponding column of M.
 * Since both M and Q are sparse, the outer product requires are relatively small number
 * of updates to the projected vectors that can be generated very efficiently.
 * Note that this algorithm isn't cache-friendly at all and will be relatively slow.
 */

void
random_indexing_sparse(double *res, int *nr, int *nc, int *p, int *row, double *x, int *n_ri, double *rate, int *verbose) {
  int target_fill = (*n_ri) * (*rate);              /* expected number of nonzero entries in each column of Q */
  int max_fill = 2 * target_fill + 1;               /* maximal number of nonzero entries */
  double *Q_x = R_alloc(max_fill, sizeof(double));  /* buffer for nonzero entries of Q (transient memory) */
  int *Q_row = R_alloc(max_fill, sizeof(int));      /* buffer for row offsets of nonzero entries (transient memory) */
  int *Q_norms = R_alloc(*n_ri, sizeof(double));    /* accumulate Euclidean norms of dynamically generated basis vectors */
  int i, j, col, rd, nnzero;
  double n_updates = 0.0;                           /* performance statistics for <verbose> mode */

  for (i = 0; i < *n_ri; i++) Q_norms[i] = 0.0;     /* initialise array of basis vector norms */
  for (i = 0; i < (*nr) * (*n_ri); i++) res[i] = 0.0;  /* initialise matrix of projected vectors */

  for (col = 0; col < *nc; col++) {
    /* generate column of Q as random vector with fill rate <rate> and values +1 / -1 with equal probability */
    GetRNGstate();
    nnzero = 0;
    rd = rgeom(*rate);                              /* generate random dimension offset from geometric distribution */
    while (rd >= *n_ri) rd = rgeom(*rate);          /* make sure there's at least one nonzero entry */
    while (rd < *n_ri && nnzero < max_fill) {
      Q_row[nnzero] = rd;
      Q_x[nnzero] = (unif_rand() >= 0.5) ? 1.0 : -1.0;
      Q_norms[rd] += 1.0;
      nnzero++;
      rd += rgeom(*rate) + 1;
    }
    PutRNGstate();
  
    /* now generate sparse outer product of matching columns of M and Q, and update <res> */
    for (i = p[col]; i < p[col+1]; i++) {
      for (j = 0; j < nnzero; j++) {
        res[ row[i] + (*nr) * Q_row[j] ] += Q_x[j] * x[i];  /* res[k, m] = Q[m, col] * M[k, col] ... sparse outer product */
        n_updates++;
      }
    }
    
    if (*verbose && ((col+1) % 100000) == 0)
      Rprintf("%6.0fk columns processed (%.1fG memory updates)\n", (col+1.0)/1000, n_updates / 1e9);
  }
  if (*verbose) Rprintf("%.1fG memory updates complete, rescaling RI dimensions\n", n_updates / 1e9);
  
  /* rescale columns of <res>, i.e. random dimensions, with Euclidean norm of basis vectors */
  for (rd = 0; rd < *n_ri; rd++) {
    double factor = 1.0 / sqrt(Q_norms[rd]);
    for (i = 0; i < *nr; i++) {
      res[ i + (*nr) * rd ] *= factor;
    }
  }
  
  /* NB: transient buffers allocated at start of function will be freed by R upon return */
}
