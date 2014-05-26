#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <limits.h>

#ifndef wordspace_defs_h
#define wordspace_defs_h

extern int have_openmp;
extern int openmp_threads;

void similarity_to_distance(double *M, int *n, int *code, double *tolerance, int *out_of_range);

void random_indexing_sparse(double *res, int *nr, int *nc, int *p, int *row, double *x, int *n_ri, double *rate, int *verbose);

#endif /* wordspace_defs_h */
