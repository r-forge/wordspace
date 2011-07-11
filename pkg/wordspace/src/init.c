#include "function_defs.h"

#include <R_ext/Rdynload.h>

const R_NativePrimitiveArgType sqrt_args[3] = {REALSXP, INTSXP, REALSXP};

const R_NativePrimitiveArgType score_dense_args[10] = {REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, INTSXP};
const R_NativePrimitiveArgType score_sparse_args[11] = {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, INTSXP};

const R_NativePrimitiveArgType row_norms_dense_args[6] = {REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};
const R_NativePrimitiveArgType row_norms_sparse_args[8] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};

const R_CMethodDef cMethods[] = {
  {"C_sqrt", (DL_FUNC) &do_sqrt, 3, sqrt_args},
  {"C_dsm_score_dense", (DL_FUNC) &dsm_score_dense, 10, score_dense_args},
  {"C_dsm_score_sparse", (DL_FUNC) &dsm_score_sparse, 11, score_sparse_args},
  {"C_row_norms_dense", (DL_FUNC) &row_norms_dense, 6, row_norms_dense_args},
  {"C_row_norms_sparse", (DL_FUNC) &row_norms_sparse, 8, row_norms_sparse_args},
  {NULL, NULL, 0}
};

void
R_init_wordspace(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}
