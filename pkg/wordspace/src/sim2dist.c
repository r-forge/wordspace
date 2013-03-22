/*
 *  Transform similarity values to distances (in-place operation on vector or dense matrix)
 */

/* internal codes for type of transformation
 *  0 = cosine -> angle
*/

#include "function_defs.h"

void
similarity_to_distance(double *M, int *n, int *code, double *tolerance, int *out_of_range) {
  int n_clamped = 0, opcode = *code, n_items = *n;
  int i;
  double tol = *tolerance, x, *Mptr;
  
  if (opcode < 0 || opcode > 0)
    error("transformation method #%d is not defined -- internal error", opcode);
  
  for (i = 0; i < n_items; i++) {
    x = M[i];
    switch (opcode) {
      case 0:
        if (x < -(1-tol)) {
          if (x < -(1+tol)) n_clamped++;
          x = -1;
        }
        else if (x > (1-tol)) {
          if (x > (1+tol)) n_clamped++;
          x = 1;
        }
        x = acos(x) * 180 / M_PI; 
        break;
    }
    M[i] = x;
  }
  
  *out_of_range = n_clamped;
}
