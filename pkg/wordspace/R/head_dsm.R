head.dsm <- function (x, n=6L, k=n, ...) {
  info <- check.dsm(x)
  n <- min(n, info$nrow)
  k <- min(k, info$ncol)
  M <- if (info$S$ok) x$S else x$M
  M[1:n, 1:k]
}
