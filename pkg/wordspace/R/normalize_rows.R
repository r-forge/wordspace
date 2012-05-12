normalize.rows <- function (M, ...) {
  if (inherits(M, "Matrix")) {
    # sparse matrix: generic arithmetic method M / x constructs full dense matrix internally (as of R 2.15.0),
    # so we need to implement explicit code for row scaling (could be made even more memory-efficient with C function)
    norms <- rowNorms(M, ...)
    M@x <- M@x / norms[M@i + 1]
    M
  } else if (is.matrix(M)) {
    # dense matrix
    M / rowNorms(M, ...) # faster than sweep(M, 1, rowNorms(M, ...), "/") 
  } else {
    stop("M must be a dense or sparse matrix")
  }
}

normalize.rows.old <- function (M, ...) {
  M / rowNorms.old(M, ...) # faster than sweep(M, 1, rowNorms(M, ...), "/") and works for sparse matrices as well
}