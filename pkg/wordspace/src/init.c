#include "function_defs.h"

#include <R_ext/Rdynload.h>

const R_NativePrimitiveArgType sqrt_args[3] = {REALSXP, INTSXP, REALSXP};

const R_CMethodDef cMethods[] = {
  {"C_sqrt", (DL_FUNC) &do_sqrt, 3, sqrt_args},
  {NULL, NULL, 0}
};

void
R_init_wordspace(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
}
