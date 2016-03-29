#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <R.h>
#include <Rinternals.h>

#include "svdlib.h"

SEXP svdLAS2_(SEXP dim, SEXP i, SEXP p, SEXP x, SEXP dimensions, SEXP iterations, SEXP exclude, SEXP kappa) {
  struct smat M;
  SVDRec svd;
  SEXP res, res_d, res_u, res_v, res_names;

  double *u_dbl, *v_dbl, *mark, *point;
  int *i_int = INTEGER(i);
  int *p_int = INTEGER(p);
  
  int nR = INTEGER(dim)[0];
  int nC = INTEGER(dim)[1];
  int n_cells = length(x);
  int rank, k, j, n_row, n_col;
  
  /* copy M to SMat structure (column-compressed format) */
  M.rows = nR;
  M.cols = nC;
  M.vals = n_cells;
  M.value = REAL(x);
  /* need to make copy of i and p because of different data type (long vs. int) */
  M.pointr = (long *) R_alloc(nC + 1, sizeof(long));
  for (k = 0; k <= nR; k++)
    M.pointr[k] = p_int[k];
  M.rowind = (long *) R_alloc(n_cells, sizeof(long));
  for (k = 0; k < n_cells; k++)
    M.rowind[k] = i_int[k];

  /* execute sparse SVD */
  SVDVerbosity = 0;
  svd = svdLAS2A(&M, 4);
  rank = svd->d;

  /* check matrix dimensions */
  n_row = svd->Ut->cols; /* Ut is the transposed matrix, hence swap row/col counts */
  n_col = svd->Ut->rows;
  if ((n_col != rank) || (n_row != nR)) {
    svdFreeSVDRec(svd);
    error("interal error (U is %d x %d matrix, expected %d x %d)", n_row, n_col, nR, rank);
  }
  n_row = svd->Vt->cols;
  n_col = svd->Vt->rows;
  if ((n_col != rank) || (n_row != nC)) {
    svdFreeSVDRec(svd);
    error("interal error (V is %d x %d matrix, expected %d x %d)", n_row, n_col, nC, rank);
  }
  
  /* extract singular values and matrices of singluar vectors into R objects */
  res_d = PROTECT(allocVector(REALSXP, rank));
  for (k = 0; k < rank; k++)
    REAL(res_d)[k] = (svd->S)[k];
  res_u = PROTECT(allocMatrix(REALSXP, nR, rank));
  u_dbl = REAL(res_u);
  for (k = 0; k < rank; k++) {
    mark = (svd->Ut->value)[k];
    point = u_dbl + k * nR;
    for (j = 0; j < nR; j++)
      *point++ = *mark++;
  }
  res_v = PROTECT(allocMatrix(REALSXP, nC, rank));
  v_dbl = REAL(res_v);
  for (k = 0; k < rank; k++) {
    mark = (svd->Vt->value)[k];
    point = v_dbl + k * nC;
    for (j = 0; j < nC; j++)
      *point++ = *mark++;
  }

  /* free SVDRec after copying to R objects */
  svdFreeSVDRec(svd);

  /* construct result list */
  res = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, res_d);
  SET_VECTOR_ELT(res, 1, res_u);
  SET_VECTOR_ELT(res, 2, res_v);
  res_names = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(res_names, 0, mkChar("d"));
  SET_STRING_ELT(res_names, 1, mkChar("u"));
  SET_STRING_ELT(res_names, 2, mkChar("v"));
  setAttrib(res, R_NamesSymbol, res_names);
  
  UNPROTECT(5);
  return res;
}
