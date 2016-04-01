normalize.rows <- function (M, ...) {
  scaleMargins(M, rows = 1 / rowNorms(M, ...))
}
