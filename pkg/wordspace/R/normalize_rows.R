normalize.rows <- function (M, ...) {
  M / rowNorms(M, ...) # faster than sweep(M, 1, rowNorms(M, ...), "/") and works for sparse matrices as well
}