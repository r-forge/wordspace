normalize.rows <- function (M, ...) {
  scaleMargins(M, rows=1/rowNorms(M, ...))
}

normalize.rows.old <- function (M, ...) {
  M / rowNorms.old(M, ...) # faster than sweep(M, 1, rowNorms(M, ...), "/") and works for sparse matrices as well
}