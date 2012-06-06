subset.dsm <- function (x, subset, select, drop.zeroes=FALSE, matrix.only=FALSE, ...) {
  info <- check.dsm(x, validate=TRUE) # make sure that rows/columns are consistent

  # NB: subset and select conditions cannot be disabled by specifying NULL, since is.null()
  #     would force evaluation of the expression in the current environment (!= eval below)
  if (missing(subset)) {
    row.idx <- 1:info$nrow
    update.nz.cols <- FALSE
  } else {
    condition <- substitute(subset)
    row.idx <- eval(condition, c(x$rows, x$globals), parent.frame())
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
    update.nz.cols <- TRUE # rows may have been deleted, so nnzero counts for columns need to be updated
  }
  
  if (missing(select)) {
    col.idx <- 1:info$ncol
    update.nz.rows <- FALSE
  } else {
    condition <- substitute(select)
    col.idx <- eval(condition, c(x$cols, x$globals), parent.frame())
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
    update.nz.rows <- TRUE # columns may have been deleted, so nnzero counts for rows need to be updated
  }

  if (is.logical(row.idx)) row.idx <- which(row.idx) # make sure we have numeric indices for additional subsetting
  if (is.logical(col.idx)) col.idx <- which(col.idx)

  if (drop.zeroes) {
    M <- if (info$have.S) x$S else x$M  # primary data matrix (use scores if available, which may be sparser than frequencies)
    is.nzero <- M[row.idx, col.idx, drop=FALSE] != 0
    nnzero.rows <- rowSums(is.nzero)
    nnzero.cols <- colSums(is.nzero)
    rm(M, is.nzero)
    keep.rows <- nnzero.rows > 0
    row.idx <- row.idx[keep.rows] # drop rows without nonzero entries
    nnzero.rows <- nnzero.rows[keep.rows] # updated nnzero counts for final subset
    keep.cols <- nnzero.cols > 0
    col.idx <- col.idx[keep.cols] # drop columns without nonzero entries
    nnzero.cols <- nnzero.cols[keep.cols]
    rm(keep.rows, keep.cols)
  }

  if (matrix.only) {
    M <- if (info$have.S) x$S else x$M  # matrix.only=TRUE: just return subset of the appropriate matrix
    return(M[row.idx, col.idx, drop=FALSE])
  }

  # for small result sets, it's much more memory-efficient to construct a new DSM object from scratch
  y <- list(rows=x$rows[row.idx, , drop=FALSE],   # mandatory components
            cols=x$cols[col.idx, , drop=FALSE],
            globals=x$globals, locked=x$locked)
  if ("N" %in% info$slots) y$N <- x$N             # optional components
  if (info$have.M) y$M <- x$M[row.idx, col.idx, drop=FALSE]
  if (info$have.S) y$S <- x$S[row.idx, col.idx, drop=FALSE]

  # if rows and/or columns may have been deleted, update the relevant nnzero counts
  if (drop.zeroes) {
    y$rows$nnzero <- nnzero.rows  # we've already computed these above
    y$cols$nnzero <- nnzero.cols
  } else {
    if (update.nz.rows || update.nz.cols) {
      is.nzero <- if (info$have.S) y$S != 0 else y$M != 0
      if (update.nz.rows) y$rows$nnzero <- rowSums(is.nzero)
      if (update.nz.cols) y$cols$nnzero <- colSums(is.nzero)
    }
  }

  class(y) <- c("dsm", "list")
  return(y)
}
