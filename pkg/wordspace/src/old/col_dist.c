/*
 *  Compute distances between columns of two dense or sparse matrices
 */

/* internal codes for metric / distance:
 *  0 = euclidean
 *  1 = maximum
 *  2 = manhattan
 *  3 = minkowski   (*param1 = exponent p)
 *  4 = canberra
 */

#include "function_defs.h"

/* make symmetric matrix from right upper triangle */
void
mk_symmetric_matrix(double *x, int nr, int nc) {
  int r, c;
  for (c = 0; c < nc; c++)
    for (r = 0; r < c; r++)
      x[nr * r + c] = x[nr * c + r]; /* x[c, r] = x[r, c] */
}

void
col_dist_dense(double *dist, int *nr, int *nc1, int *nc2, double *x, double *y, int *metric_code, double *param1, int *symmetric) {
  int row, col1, col2, vec_len, col1_max;
  double accum, d_xy;
  double *dist_ptr, *x_ptr, *y_ptr;

  if (*metric_code < 0 || *metric_code > 4)
    error("distance metric #%d is not defined -- internal error", *metric_code);
  if (*metric_code == 3 && (*param1 < 1 || !R_FINITE(*param1)))
    error("Minkowski p-norm is not defined for p = %g", *param1);

  vec_len = *nr;
#pragma omp parallel for \
        if (openmp_threads > 1 && (*nc1 + 0.0) * (*nc2 + 0.0) * (*nr + 0.0) > 100e6) \
        num_threads(openmp_threads) \
        shared(dist) \
        private(row, col1, col2, col1_max, accum, d_xy, dist_ptr, x_ptr, y_ptr)
  for (col2 = 0; col2 < *nc2; col2++) {
    dist_ptr = dist + *nc1 * col2;      /* column <col2> of result matrix <dist> */
    col1_max = (*symmetric) ? col2 + 1 : *nc1;
    for (col1 = 0; col1 < col1_max; col1++) {
      x_ptr = x + vec_len * col1;       /* column <col1> of first matrix <x> */
      y_ptr = y + vec_len * col2;       /* column <col2> of second matrix <y> */
      accum = 0;
      switch (*metric_code) {
      case 0:
        for (row = 0; row < vec_len; row++) {
          d_xy = *(x_ptr++) - *(y_ptr++);
          accum += d_xy * d_xy;
        }
        *(dist_ptr++) = sqrt(accum);
        break;
      case 1:
        for (row = 0; row < vec_len; row++) {
          d_xy = fabs(*(x_ptr++) - *(y_ptr++));
          if (d_xy > accum) accum = d_xy;
        }
        *(dist_ptr++) = accum;
        break;
      case 2:
        for (row = 0; row < vec_len; row++) {
          d_xy = fabs(*(x_ptr++) - *(y_ptr++));
          accum += d_xy;
        }
        *(dist_ptr++) = accum;
        break;
      case 3:
        for (row = 0; row < vec_len; row++) {
          d_xy = fabs(*(x_ptr++) - *(y_ptr++));
          accum += pow(d_xy, *param1);
        }
        *(dist_ptr++) = pow(accum, 1 / *param1);
        break;
      case 4:
        for (row = 0; row < vec_len; row++) {
          double x_plus_y = fabs(*x_ptr) + fabs(*y_ptr);
          d_xy = fabs(*(x_ptr++) - *(y_ptr++));
          if (x_plus_y > 0) accum += d_xy / x_plus_y;
        }
        *(dist_ptr++) = accum;
        break;
      }
    }
  }
  
  if (*symmetric) mk_symmetric_matrix(dist, *nc1, *nc2);
}

void
col_dist_sparse(double *dist, int *nc1, int *nc2, int *xp, int *xrow, double *x, int *yp, int *yrow, double *y, int *metric_code, double *param1, int *symmetric) {
  int xrow_curr, yrow_curr, col1, col2, col1_max;
  int xi, yi, xi_max, yi_max;
  double *dist_ptr, *x_ptr, *y_ptr;
  double accum, d_xy, x_plus_y, x_curr, y_curr;
  /* average number of entries scanned when comparing two columns */
  double avg_nr = (xp[*nc1] - xp[0] + 0.0) / *nc1 + (yp[*nc2] - yp[0] + 0.0) / *nc2; 

  if (*metric_code < 0 || *metric_code > 4)
    error("distance metric #%d is not defined -- internal error", *metric_code);
  if (*metric_code == 3 && (*param1 < 1 || !R_FINITE(*param1)))
    error("Minkowski p-norm is not defined for p = %g", *param1);

#pragma omp parallel for \
        if (openmp_threads > 1 && (*nc1 + 0.0) * (*nc2 + 0.0) * avg_nr > 40e6) \
        num_threads(openmp_threads) \
        shared(dist) \
        private(xrow_curr, yrow_curr, col1, col2, col1_max, xi, yi, xi_max, yi_max, accum, d_xy, x_plus_y, x_curr, y_curr, dist_ptr, x_ptr, y_ptr)
  for (col2 = 0; col2 < *nc2; col2++) {
    dist_ptr = dist + *nc1 * col2;      /* column <col2> of result matrix <dist> */
    col1_max = (*symmetric) ? col2 + 1 : *nc1;
    yi_max = yp[col2 + 1];

    for (col1 = 0; col1 < col1_max; col1++) {
      xi_max = xp[col1 + 1];
      xi = xp[col1];
      yi = yp[col2];
      xrow_curr = (xi < xi_max) ? xrow[xi] : INT_MAX;
      yrow_curr = (yi < yi_max) ? yrow[yi] : INT_MAX;
      
      accum = 0;
      while (xi < xi_max || yi < yi_max) {

        if (xrow_curr < yrow_curr) {
          x_curr = x[xi]; y_curr = 0;
          xi++;
          xrow_curr = (xi < xi_max) ? xrow[xi] : INT_MAX;
        }
        else if (xrow_curr == yrow_curr) {
          x_curr = x[xi]; y_curr = y[yi];
          xi++; yi++;
          xrow_curr = (xi < xi_max) ? xrow[xi] : INT_MAX;
          yrow_curr = (yi < yi_max) ? yrow[yi] : INT_MAX;
        }
        else {
          x_curr = 0; y_curr = y[yi];
          yi++;
          yrow_curr = (yi < yi_max) ? yrow[yi] : INT_MAX;          
        }
        
        switch (*metric_code) {
        case 0:
          d_xy = x_curr - y_curr;
          accum += d_xy * d_xy;
          break;
        case 1:
          d_xy = fabs(x_curr - y_curr);
          if (d_xy > accum) accum = d_xy;
          break;
        case 2:
          d_xy = fabs(x_curr - y_curr);
          accum += d_xy;
          break;
        case 3:
          d_xy = fabs(x_curr - y_curr);
          accum += pow(d_xy, *param1);
          break;
        case 4:
          x_plus_y = fabs(x_curr) + fabs(y_curr);
          d_xy = fabs(x_curr - y_curr);
          if (x_plus_y > 0) accum += d_xy / x_plus_y;
          break;
        }
      } /* while (xi, yi) */

      switch (*metric_code) {
      case 0:
        *(dist_ptr++) = sqrt(accum);
        break;
      case 1:
      case 2:
      case 4:
        *(dist_ptr++) = accum;
        break;
      case 3:
        *(dist_ptr++) = pow(accum, 1 / *param1);
        break;
      }

    } /* for (col1) */
  } /* for (col2) */
  
  if (*symmetric) mk_symmetric_matrix(dist, *nc1, *nc2);
}
