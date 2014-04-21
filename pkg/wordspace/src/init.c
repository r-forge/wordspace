#include "function_defs.h"

#include <R_ext/Rdynload.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int have_openmp;
int openmp_threads;

void
get_openmp_threads(int *num_threads, int *max_threads) {
  *num_threads = openmp_threads;
#ifdef _OPENMP
  *max_threads = omp_get_max_threads();
#else
  *max_threads = 0;
#endif
}

void
set_openmp_threads(int *num_threads) {
  int n = *num_threads;
  int max_threads = 0;
  if (n < 1) error("number of threads must be >= 1");
#ifdef _OPENMP
  max_threads = omp_get_max_threads();
  if (n > max_threads) n = max_threads;
  openmp_threads = n;
#else
  if (n > 1) warning("OpenMP support not available");
  openmp_threads = 1;
#endif    
}

R_NativePrimitiveArgType similarity_to_distance_args[5] = {REALSXP, INTSXP, INTSXP, REALSXP, INTSXP};

R_NativePrimitiveArgType random_indexing_sparse_args[9] = {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, LGLSXP};

R_NativePrimitiveArgType get_openmp_threads_args[2] = {INTSXP, INTSXP};
R_NativePrimitiveArgType set_openmp_threads_args[1] = {INTSXP};

const R_CMethodDef cMethods[] = {
  {"C_similarity_to_distance", (DL_FUNC) &similarity_to_distance, 5, similarity_to_distance_args},
  {"C_random_indexing_sparse", (DL_FUNC) &random_indexing_sparse, 9, random_indexing_sparse_args},
  {"C_get_openmp_threads", (DL_FUNC) &get_openmp_threads, 2, get_openmp_threads_args},
  {"C_set_openmp_threads", (DL_FUNC) &set_openmp_threads, 1, set_openmp_threads_args},
  {NULL, NULL, 0}
};

void
R_init_wordspace(DllInfo *dll) {
#ifdef _OPENMP
  have_openmp = 1;
  openmp_threads = 1; /* don't use multi-core processing unless explicitly enabled */
#else
  have_openmp = 0;
  openmp_threads = 1;
#endif
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}
