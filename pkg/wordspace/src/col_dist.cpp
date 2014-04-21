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

#include <Rcpp.h>
using namespace Rcpp;

#include "globals.h"

/* make symmetric matrix from right upper triangle */
void mk_symmetric_matrix(NumericMatrix x) {
  int nr = x.nrow(), nc = x.ncol();
  for (int c = 0; c < nc; c++)
    for (int r = 0; r < c; r++)
      x(c, r) = x(r, c);
  //      x[nr * r + c] = x[nr * c + r]; /* x[c, r] = x[r, c] */
}

void check_metric(int metric_code, double p1) {
  if (metric_code < 0 || metric_code > 4)
    stop("internal error -- invalid metric code");
  if (metric_code == 3 && (!R_FINITE(p1) || p1 < 1))
    stop("internal error -- Minkowski metric p out of range [1, Inf)");  
}

// [[Rcpp::export]]
NumericVector CPP_col_dist_dense(NumericMatrix x, NumericMatrix y, int metric_code, double param1, bool symmetric) {
  check_metric(metric_code, param1);
  int nr = x.nrow(), nc1 = x.ncol(), nc2 = y.ncol();
  if (nr != y.nrow()) stop("internal error -- matrices are not conformable");

  NumericMatrix dist(nc1, nc2);

#pragma omp parallel for \
        if (openmp_threads > 1 && (nc1 + 0.0) * (nc2 + 0.0) * (nr + 0.0) > 100e6) \
        num_threads(openmp_threads) \
        shared(dist)
  for (int col2 = 0; col2 < nc2; col2++) {
    NumericVector tmp(nr);
    double accum;
    int col1_max = (symmetric) ? col2 + 1 : nc1;
    for (int col1 = 0; col1 < col1_max; col1++) {
      NumericMatrix::Column vx = x(_, col1);  // column <col1> of first matrix
      NumericMatrix::Column vy = y(_, col2);  // column <col2> of second matrix
      switch (metric_code) {
      case 0: 
        accum = sum((vx - vy) * (vx - vy));
        dist(col1, col2) = sqrt(accum);
        break;
      case 1:
        dist(col1, col2) = max(abs(vx - vy));
        break;
      case 2:
        dist(col1, col2) = sum(abs(vx - vy));
        break;
      case 3:
        accum = sum(pow(abs(vx - vy), param1));
        dist(col1, col2) = pow(accum, 1.0 / param1);
        break;
      case 4:
        tmp = abs(vx) + abs(vy); // denominator |x_i| + |y_i|
        dist(col1, col2) = sum(ifelse(tmp > 0, abs(vx - vy) / tmp, 0.0));
        break;
      }
    }
  }
  
  if (symmetric) mk_symmetric_matrix(dist);
  return(dist);
}
