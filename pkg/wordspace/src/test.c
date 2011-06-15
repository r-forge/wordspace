#include "function_defs.h"

#include <math.h>

void
do_sqrt(double *x, int *n, double *result) {
  int i;
  for (i = 0; i < *n; i++) {
    result[i] = log10(1 + x[i]);
  }
}
