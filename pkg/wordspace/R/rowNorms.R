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

  info <- dsm.is.canonical(M)
  if (!info$canonical) M <- dsm.canonical.matrix(M)

  result <- double(nrow(M))
  if (info$sparse) {
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
  
  info <- dsm.is.canonical(M)
  if (!info$canonical) M <- dsm.canonical.matrix(M)

  result <- double(ncol(M))
  if (info$sparse) {
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
