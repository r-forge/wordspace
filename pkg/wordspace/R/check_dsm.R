check.dsm <- function (model, validate=FALSE) {
  stopifnot(inherits(model, "dsm"))
  slots <- names(model)
  stopifnot(all(c("M","rows","cols") %in% slots))
  stopifnot(all(c("term","R1","R2") %in% colnames(model$rows)))
  stopifnot(all(c("term","C1","C2") %in% colnames(model$cols)))
  have.S <- "S" %in% slots
  
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
  
  list(nrow=n.rows, ncol=n.cols, slots=slots, have.S=have.S)
}