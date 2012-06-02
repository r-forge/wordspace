#include "function_defs.h"

#include <R_ext/Rdynload.h>

R_NativePrimitiveArgType sqrt_args[3] = {REALSXP, INTSXP, REALSXP};

R_NativePrimitiveArgType score_dense_args[10] = {REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, INTSXP};
R_NativePrimitiveArgType score_sparse_args[11] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, INTSXP};

R_NativePrimitiveArgType row_norms_dense_args[6] = {REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};
R_NativePrimitiveArgType row_norms_sparse_args[8] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};

R_NativePrimitiveArgType col_dist_dense_args[9] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, LGLSXP};
R_NativePrimitiveArgType col_dist_sparse_args[12] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, LGLSXP};

R_NativePrimitiveArgType random_indexing_sparse_args[8] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};

const R_CMethodDef cMethods[] = {
  {"C_dsm_score_dense",  (DL_FUNC) &dsm_score_dense,  10, score_dense_args},
  {"C_dsm_score_sparse", (DL_FUNC) &dsm_score_sparse, 11, score_sparse_args},
  {"C_row_norms_dense",  (DL_FUNC) &row_norms_dense,   6, row_norms_dense_args},
  {"C_row_norms_sparse", (DL_FUNC) &row_norms_sparse,  8, row_norms_sparse_args},
  {"C_col_norms_dense",  (DL_FUNC) &col_norms_dense,   6, row_norms_dense_args},
  {"C_col_norms_sparse", (DL_FUNC) &col_norms_sparse,  8, row_norms_sparse_args},
  {"C_col_dist_dense",   (DL_FUNC) &col_dist_dense,    9, col_dist_dense_args},
  {"C_col_dist_sparse",  (DL_FUNC) &col_dist_sparse,  12, col_dist_sparse_args},
  {"C_random_indexing_sparse", (DL_FUNC) &random_indexing_sparse, 8, random_indexing_sparse_args},
  {NULL, NULL, 0}
};

void
R_init_wordspace(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}
