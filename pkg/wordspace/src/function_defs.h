#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>

#ifndef wordspace_defs_h
#define wordspace_defs_h

void dsm_score_dense(double *scores, int *nr, int *nc, double *f, double *f1, double *f2, double *N, int *am_code, int *sparse, int *transform_code);

void dsm_score_sparse(double *scores, int *nc, int *p, int *row, double *f, double *f1, double *f2, double *N, int *am_code, int *sparse, int *transform_code);

void row_norms_dense(double *norms, int *nr, int *nc, double *x, int *norm_code, double *p_norm);

void row_norms_sparse(double *norms, int *nr, int *nc, int *p, int *row_of, double *x, int *norm_code, double *p_norm);

void do_sqrt(double *x, int *n, double *result);

#endif /* wordspace_defs_h */