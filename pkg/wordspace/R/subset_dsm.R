subset.dsm <- function (x, subset, select, ...) {
  x.info <- check.dsm(x, validate=TRUE) # make sure that rows/columns are consistent

  if (missing(subset)) {
    row.idx <- TRUE
  } else {
    condition <- substitute(subset)
    row.idx <- eval(condition, c(x$rows, x$globals), parent.frame())
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
  }
  
  if (missing(select)) {
    col.idx <- TRUE
  } else {
    condition <- substitute(select)
    col.idx <- eval(condition, c(x$cols, x$globals), parent.frame())
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
  }

  # for small result sets, it's much more memory-efficient to construct a new DSM object
  y <- list(M=x$M[row.idx, col.idx, drop=FALSE],
            rows=x$rows[row.idx, , drop=FALSE],
            cols=x$cols[col.idx, , drop=FALSE],
            N=x$N, globals=x$globals, locked=x$locked)
  if (x.info$have.S) y$S <- x$S[row.idx, col.idx, drop=FALSE]
  class(y) <- c("dsm", "list")
  return(y)
}
