rowNorms <- function (M, method = c("euclidean", "maximum", "manhattan", "minkowski"), p = 2) {
  method <- match.arg(method)
  stopifnot(is.matrix(M))
  if (method == "minkowski" && p == Inf) {
    method <- "maximum"
    p <- 2
  }
  stopifnot(p >= 1 && p < Inf)

  switch(method,
    euclidean = sqrt(rowSums(M * M)),
    maximum = apply(M, 1, function (x) max(abs(x))),
    manhattan = rowSums(abs(M)),
    minkowski = (rowSums(abs(M) ^ p)) ^ (1/p)
  )
}
