.sparse.rowMax <- function (M) {
  # compute the maximum of each row of a sparse matrix in a reasonably efficient way
  M <- abs(M)
  idx <- which(M > 0, arr.ind=TRUE)
  tapply(1:nrow(idx), list(idx[,1]), function (k) max( M[ idx[k,,drop=FALSE] ]))
}

rowNorms <- function (M, method = c("euclidean", "maximum", "manhattan", "minkowski"), p = 2) {
  method <- match.arg(method)
  if (method == "minkowski" && p == Inf) {
    method <- "maximum"
    p <- 2
  }
  stopifnot(p >= 1 && p < Inf)

  sparse.M <- inherits(M, "Matrix")
  stopifnot(is.matrix(M) || sparse.M)

  switch(method,
    euclidean = sqrt(rowSums(M * M)),
    maximum = if (sparse.M) .sparse.rowMax(M) else apply(M, 1, function (x) max(abs(x))),
    manhattan = rowSums(abs(M)),
    minkowski = (rowSums(abs(M) ^ p)) ^ (1/p)
  )
}
