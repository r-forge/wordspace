## TODO: write C code for more memory efficient row & column scaling of a sparse or dense matrix
scaleMargins <- function (M, rows=NULL, cols=NULL) {
  info <- dsm.is.canonical(M)
  if (!info$canonical) M <- dsm.canonical.matrix(M)

  if (!is.null(rows)) {
    if (length(rows) == 1) rows <- rep(rows, nrow(M))
    stopifnot(length(rows) == nrow(M))
  }
  if (!is.null(cols)) {
    if (length(cols) == 1) cols <- rep(cols, ncol(M))
    stopifnot(length(cols) == ncol(M))
  }
  
  if (info$sparse) {
    # sparse matrix: as of R 2.15.0, Matrix package constructs a full dense matrix internally when computing M / x,
    # so we need to implement explicit code for row scaling (could be made even more memory-efficient with C function)
    if (!is.null(rows)) M@x <- M@x * rows[M@i + 1] # M@i == row numbers of nonzero entries (0-based)
    if (!is.null(cols)) {
      col.idx <- rep(1:ncol(M), diff(M@p)) # compute column numbers of nonzero entries
      M@x <- M@x * cols[col.idx]
    }
  } else {
    if (!is.null(rows)) M <- M * rows # matrix * vector scales rows of matrix
    if (!is.null(cols)) M <- scale(M, center=FALSE, scale=1/cols) # scale() divides by scale= argument (NB: 1/0 = Inf, x/Inf = 0)
  }
  
  M
}
