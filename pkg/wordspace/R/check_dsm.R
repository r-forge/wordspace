check.dsm <- function (model, validate=FALSE) {
  stopifnot(inherits(model, "dsm"))
  slots <- names(model)
  stopifnot(all(c("M","rows","cols","N") %in% slots))
  stopifnot(all(c("term","f") %in% colnames(model$rows)))
  stopifnot(all(c("term","f") %in% colnames(model$cols)))
  have.S <- "S" %in% slots
  is.locked <- if ("locked" %in% slots) model$locked else FALSE
  main.matrix <- if (have.S) model$S else model$M
  is.sparse <- inherits(main.matrix, "Matrix")

  n.rows <- nrow(model$M)
  n.cols <- ncol(model$M)
  stopifnot(nrow(model$rows) == n.rows)
  stopifnot(nrow(model$cols) == n.cols)
  if (have.S) {
    stopifnot(nrow(model$S) == n.rows)
    stopifnot(ncol(model$S) == n.cols)
  }

  if (validate) {
    stopifnot(all(rownames(model$M) == model$rows$term))
    stopifnot(all(colnames(model$M) == model$cols$term))
    if (have.S) {
      stopifnot(all(rownames(model$S) == model$rows$term))
      stopifnot(all(colnames(model$S) == model$cols$term))
    }
  }
  
  list(nrow=n.rows, ncol=n.cols, N=model$N, slots=slots, have.S=have.S, locked=is.locked, sparse=is.sparse)
}
