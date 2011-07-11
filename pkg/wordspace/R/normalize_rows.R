normalize.rows <- function (M, ...) {
  M / rowNorms(M, ...) # faster than sweep(M, 1, rowNorms(M, ...), "/") and works for sparse matrices as well
}

normalize.rows.old <- function (M, ...) {
  M / rowNorms.old(M, ...) # faster than sweep(M, 1, rowNorms(M, ...), "/") and works for sparse matrices as well
}