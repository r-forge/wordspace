.check.norm <- function (method, p) {
  if (method == "minkowski" && p == Inf) {
    method <- "maximum"
    p <- 2
  }
  if (method == "minkowski" && (p < 1 || !is.finite(p))) stop("Minkowski p-norm only defined for 1 <= p < Inf")
  
  ## internal codes for selected norm (must match C code in <row_norms.c>)
  code <- switch(method, euclidean=0, maximum=1, manhattan=2, minkowski=3)
  if (is.null(code)) stop("unknown norm selected (internal error)")
  
  list(code=as.integer(code), p=as.double(p))
}

rowNorms <- function (M, method = c("euclidean", "maximum", "manhattan", "minkowski"), p = 2) {
  method <- match.arg(method)
  norm <- .check.norm(method, p)
  
  sparse.M <- inherits(M, "Matrix")
  if (!(is.matrix(M) || sparse.M)) stop("first argumnet must be dense or sparse matrix")
  if (sparse.M && !is(M, "dgCMatrix")) stop("sparse matrix must be in normal form (dgCMatrix)")

  result <- double(nrow(M))
  if (sparse.M) {
    .C(
      C_row_norms_sparse,
      result,
      as.integer(nrow(M)),
      as.integer(ncol(M)),
      as.integer(M@p),
      as.integer(M@i),
      as.double(M@x),
      norm$code,
      norm$p,
      DUP=FALSE, NAOK=FALSE
    )
  } else {
    .C(
      C_row_norms_dense,
      result,
      as.integer(nrow(M)),
      as.integer(ncol(M)),
      as.double(M),
      norm$code,
      norm$p,
      DUP=FALSE, NAOK=FALSE
    )
  }

  names(result) <- rownames(M)
  result
}

colNorms <- function (M, method = c("euclidean", "maximum", "manhattan", "minkowski"), p = 2) {
  method <- match.arg(method)
  norm <- .check.norm(method, p)
  
  sparse.M <- inherits(M, "Matrix")
  if (!(is.matrix(M) || sparse.M)) stop("first argumnet must be dense or sparse matrix")
  if (sparse.M && !is(M, "dgCMatrix")) stop("sparse matrix must be in normal form (dgCMatrix)")

  result <- double(ncol(M))
  if (sparse.M) {
    .C(
      C_col_norms_sparse,
      result,
      as.integer(nrow(M)),
      as.integer(ncol(M)),
      as.integer(M@p),
      as.integer(M@i),
      as.double(M@x),
      norm$code,
      norm$p,
      DUP=FALSE, NAOK=FALSE
    )
  } else {
    .C(
      C_col_norms_dense,
      result,
      as.integer(nrow(M)),
      as.integer(ncol(M)),
      as.double(M),
      norm$code,
      norm$p,
      DUP=FALSE, NAOK=FALSE
    )
  }

  names(result) <- colnames(M)
  result
}

.sparse.rowMax <- function (M) {
  # compute the maximum of each row of a sparse matrix in a reasonably efficient way
  M <- abs(M)
  idx <- which(M > 0, arr.ind=TRUE)
  tapply(1:nrow(idx), list(idx[,1]), function (k) max( M[ idx[k,,drop=FALSE] ]))
}

rowNorms.old <- function (M, method = c("euclidean", "maximum", "manhattan", "minkowski"), p = 2) {
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
