#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <limits.h>

#ifndef wordspace_defs_h
#define wordspace_defs_h

void dsm_score_dense(double *scores, int *nr, int *nc, double *f, double *f1, double *f2, double *N, int *am_code, int *sparse, int *transform_code);
void dsm_score_sparse(double *scores, int *nc, int *p, int *row, double *f, double *f1, double *f2, double *N, int *am_code, int *sparse, int *transform_code);

void row_norms_dense(double *norms, int *nr, int *nc, double *x, int *norm_code, double *p_norm);
void row_norms_sparse(double *norms, int *nr, int *nc, int *p, int *row_of, double *x, int *norm_code, double *p_norm);

void col_norms_dense(double *norms, int *nr, int *nc, double *x, int *norm_code, double *p_norm);
void col_norms_sparse(double *norms, int *nr, int *nc, int *p, int *row_of, double *x, int *norm_code, double *p_norm);

void col_dist_dense(double *dist, int *nr, int *nc1, int *nc2, double *x, double *y, int *metric_code, double *param1, int *symmetric);
void col_dist_sparse(double *dist, int *nc1, int *nc2, int *xp, int *xrow, double *x, int *yp, int *yrow, double *y, int *metric_code, double *param1, int *symmetric);

void similarity_to_distance(double *M, int *n, int *code, double *tolerance, int *out_of_range);

void random_indexing_sparse(double *res, int *nr, int *nc, int *p, int *row, double *x, int *n_ri, double *rate, int *verbose);

#endif /* wordspace_defs_h */