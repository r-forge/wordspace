normalize.rows <- function (M, ...) {
  M / rowNorms(M, ...) # faster than sweep(M, 1, rowNorms(M, ...), "/")
  ## TODO: check whether this works for sparse matrices
}