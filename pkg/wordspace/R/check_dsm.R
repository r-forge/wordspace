check.dsm <- function (model, validate=FALSE) {
  stopifnot(inherits(model, "dsm"))
  slots <- names(model)
  have.M <- "M" %in% slots
  have.S <- "S" %in% slots
  stopifnot(have.M || have.S) # need to have either frequency matrix or score matrix (or both)
  stopifnot(all(c("rows","cols") %in% slots))
  required <- if (have.M) c("term", "f") else c("term") # required columns in $rows and $cols
  stopifnot(all(required %in% colnames(model$rows)))
  stopifnot(all(required %in% colnames(model$cols)))
  if (have.M) {
    stopifnot("N" %in% slots) # frequency matrix also requires sample size information
    N <- model$N
  } else {
    N <- NA # sample size is irrelevant if there is only a score matrix
  }
  is.locked <- if ("locked" %in% slots) model$locked else FALSE
  n.rows <- nrow(model$rows)
  n.cols <- nrow(model$cols)
  
  is.sparse.M <- FALSE
  if (have.M) {
    is.sparse.M <- inherits(model$M, "Matrix")
    if (is.sparse.M) {
      if (!is(model$M, "dgCMatrix")) stop("sparse cooccurrence matrix M must be in normal form (dgCMatrix)")
    } else {
      if (!is.matrix(model$M)) stop("cooccurrence matrix M must be a dense or sparse matrix")
    }
    stopifnot(nrow(model$M) == n.rows)
    stopifnot(ncol(model$M) == n.cols)
    if (validate) {
      stopifnot(all(rownames(model$M) == model$rows$term))
      stopifnot(all(colnames(model$M) == model$cols$term))
    }
  }

  is.sparse.S <- FALSE
  if (have.S) {
    is.sparse.S <- inherits(model$S, "Matrix")
    if (is.sparse.S) {
      if (!is(model$S, "dgCMatrix")) stop("sparse score matrix S must be in normal form (dgCMatrix)")
    } else {
      if (!is.matrix(model$S)) stop("score matrix S must be a dense or sparse matrix")
    }
    stopifnot(nrow(model$S) == n.rows)
    stopifnot(ncol(model$S) == n.cols)
    if (validate) {
      stopifnot(all(rownames(model$S) == model$rows$term))
      stopifnot(all(colnames(model$S) == model$cols$term))
    }
  }

  if (have.M && have.S) {
    if (is.sparse.M != is.sparse.S) stop("cooccurrence matrix M and score matrix S must either both be sparse or both be dense")
  }
  is.sparse <- is.sparse.M || is.sparse.S
  
  list(nrow=n.rows, ncol=n.cols, N=N, slots=slots, have.M=have.M, have.S=have.S, locked=is.locked, sparse=is.sparse)
}
