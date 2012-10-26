t.dsm <- function (x) {
  .info <- check.dsm(x)
  tx <- list(rows=x$cols, cols=x$rows, globals=x$globals, locked=x$locked)
  if (!is.na(.info$N)) tx$N <- .info$N
  if (.info$have.M) tx$M <- t(x$M)
  if (.info$have.S) tx$S <- t(x$S)
  class(tx) <- c("dsm", "list")
  return(tx)
}
