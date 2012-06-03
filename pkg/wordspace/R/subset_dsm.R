subset.dsm <- function (x, subset, select, ...) {
  info <- check.dsm(x, validate=TRUE) # make sure that rows/columns are consistent

  if (missing(subset) || is.null(subset)) {
    row.idx <- TRUE
    update.nz.cols <- FALSE
  } else {
    condition <- substitute(subset)
    row.idx <- eval(condition, c(x$rows, x$globals), parent.frame())
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
    update.nz.cols <- TRUE # rows may have been deleted, so nnzero counts for columns need to be updated
  }
  
  if (missing(select) || is.null(select)) {
    col.idx <- TRUE
    update.nz.rows <- FALSE
  } else {
    condition <- substitute(select)
    col.idx <- eval(condition, c(x$cols, x$globals), parent.frame())
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
    update.nz.rows <- TRUE # columns may have been deleted, so nnzero counts for rows need to be updated
  }

  # for small result sets, it's much more memory-efficient to construct a new DSM object from scratch
  y <- list(rows=x$rows[row.idx, , drop=FALSE],   # mandatory components
            cols=x$cols[col.idx, , drop=FALSE],
            globals=x$globals, locked=x$locked)
  if ("N" %in% info$slots) y$N <- x$N             # optional components
  if (info$have.M) y$M <- x$M[row.idx, col.idx, drop=FALSE]
  if (info$have.S) y$S <- x$S[row.idx, col.idx, drop=FALSE]

  # if rows and/or columns may have been deleted, update the relevant nnzero counts
  if (update.nz.rows || update.nz.cols) {
    is.nzero <- if (info$have.S) y$S != 0 else y$M != 0
    if (update.nz.rows) y$rows$nnzero <- rowSums(is.nzero)
    if (update.nz.cols) y$cols$nnzero <- colSums(is.nzero)
  }

  class(y) <- c("dsm", "list")
  return(y)
}
