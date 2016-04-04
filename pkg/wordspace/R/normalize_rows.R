normalize.rows <- function (M, ...) {
  scaleMargins(M, rows = 1 / rowNorms(M, ...))
}

normalize.cols <- function (M, ...) {
  scaleMargins(M, cols = 1 / colNorms(M, ...))
}
