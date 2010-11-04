subset.dsm <- function (x, subset, select, ...) {
  x.info <- check.dsm(x, validate=TRUE) # make sure that rows/columns are consistent

  if (missing(subset)) {
    row.idx <- TRUE
  } else {
    condition <- substitute(subset)
    row.idx <- eval(condition, x$rows, parent.frame())
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
  }
  
  if (missing(select)) {
    col.idx <- TRUE
  } else {
    condition <- substitute(select)
    col.idx <- eval(condition, x$cols, parent.frame())
    ## todo: check validity (either Boolean of correct length or numeric vector with indexes in range)
  }

  x$M <- x$M[row.idx, col.idx, drop=FALSE]
  x$rows <- x$rows[row.idx, , drop=FALSE]
  x$cols <- x$cols[col.idx, , drop=FALSE]
  if (x.info$have.S) x$S <- x$S[row.idx, col.idx, drop=FALSE]

  return(x)
}
