/*
 *  Compute association scores for dense or sparse matrix
 */

#include "function_defs.h"
#include "am.h"

void
dsm_score_dense(double *scores, int *nr, int *nc, double *f, double *f1, double *f2, double *N, int *am_code, int *sparse, int *transform_code) {
  int row, col, i;
  
  if (*am_code < 0 || *am_code >= am_table_entries)
    error("association measure #%d is not defined -- internal error", *am_code);
  am_func AM = am_table[*am_code]; /* selected association measure */

  i = 0;
  for (col = 0; col < *nc; col++) for (row = 0; row < *nr; row++) {
    double score = AM(f[i], f1[row], f2[col], *N, *sparse);
    scores[i] = (*transform_code) ? transform(score, *transform_code) : score;
    i++;
  }
}

void
dsm_score_sparse(double *scores, int *nc, int *p, int *row, double *f, double *f1, double *f2, double *N, int *am_code, int *sparse, int *transform_code) {
  int col, i;
  
  if (*am_code < 0 || *am_code >= am_table_entries)
    error("association measure #%d is not defined -- internal error", *am_code);
  am_func AM = am_table[*am_code]; /* selected association measure */

  if (! *sparse) error("only sparse association scores can be used with sparse matrix representation");
  
  for (col = 0; col < *nc; col++) {
    for (i = p[col]; i < p[col+1]; i++) {
      double score = AM(f[i], f1[row[i]], f2[col], *N, 1);
      scores[i] = (*transform_code) ? transform(score, *transform_code) : score;
    }
  }
}
