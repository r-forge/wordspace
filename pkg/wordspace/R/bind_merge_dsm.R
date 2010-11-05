rbind.dsm <- function (..., term.suffix=NULL, deparse.level=1) {
  models <- list(...) # should be one or more objects of class "dsm" if rbind() dispatches here
  if (!missing(term.suffix) && length(term.suffix) != length(models)) stop("term.suffix must provide as many strings as there are DSMs")

  models.info <- lapply(models, check.dsm, validate=TRUE) # validate models and extract dimensions etc.

  have.S.vec <- sapply(models.info, function (i) i$have.S) # are scored matrices available
  have.S <- all(have.S.vec)
  if (any(have.S.vec) && !have.S.vec) warning("some but not all DSM objects contain score matrix S, dropped from result")

  cols.merged <- .combine.marginals(lapply(models, function (m) m$cols), margin="column", mode="same") # check feature dimensions, then combine column marginals
  marginals.inconsistent <- attr(cols.merged, "adjusted") # if marginal differ between DSMs, they're adjusted to the maximum value to ensure consistency
  if (marginals.inconsistent && !have.S) warning("DSM objects have inconsistent column marginals, should calculate scores before combining them")

  rows.merged <- .bind.marginals(lapply(models, function (m) m$rows), margin="row", term.suffix=term.suffix)
  
  res <- list(
    M = do.call(rbind, lapply(models, function (m) m$M)),
    rows = rows.merged,
    cols = cols.merged,
    locked = marginals.inconsistent
  )
  if (have.S) {
    res$S <- do.call(rbind, lapply(models, function (m) m$S))
  }
  class(res) <- c("dsm", "list")

  return(res)
}

cbind.dsm <- function (..., term.suffix=NULL, deparse.level=1) {
  stop("not yet implemented")
}

merge.dsm <- function (x, y, ..., rows=TRUE, all=FALSE, term.suffix=NULL) {
  models <- list(x, y, ...)
  n.models <- length(models)
  if (!missing(term.suffix) && length(term.suffix) != n.models) stop("term.suffix must provide as many strings as there are DSMs")
  models.info <- lapply(models, check.dsm, validate=TRUE) # validate models and extract dimensions etc.

  if (all) stop("all=TRUE is not yet implemented")
  if (!rows) stop("rows=FALSE is not yet implemented")

  have.S.vec <- sapply(models.info, function (i) i$have.S) # are scored matrices available
  have.S <- all(have.S.vec)
  if (any(have.S.vec) && !have.S.vec) warning("some but not all DSM objects contain score matrix S, dropped from result")

  # bind rows of DSM objects, preserving only terms that are shared by all DSM
  cols.merged <- .combine.marginals(lapply(models, function (m) m$cols), margin="column", mode="intersect")
  marginals.inconsistent <- attr(cols.merged, "adjusted") # if marginal differ between DSMs, they're adjusted to the maximum value to ensure consistency
  if (marginals.inconsistent && !have.S) warning("DSM objects have inconsistent column marginals, should calculate scores before combining them")
  
  rows.merged <- .bind.marginals(lapply(models, function (m) m$rows), margin="row", term.suffix=term.suffix)
  n.rows <- sapply(models.info, function (i) i$nrow)
  first.row <- cumsum(c(1, n.rows)) # rows offsets of individual DSMs in combined matrix M
  
  M <- matrix(nrow=nrow(rows.merged), ncol=nrow(cols.merged), dimnames=list(rows.merged$term, cols.merged$term))
  if (have.S) S <- matrix(nrow=nrow(rows.merged), ncol=nrow(cols.merged), dimnames=list(rows.merged$term, cols.merged$term))
  for (i in 1:n.models) {
    model <- models[[i]]
    col.idx <- na.omit( match(cols.merged$term, model$cols$term) ) # extract columns of i-th DSM matrix that match the common terms
    M[ (first.row[i]):(first.row[i]+n.rows[i]-1), ] <- model$M[ , col.idx]
    if (have.S) S[ (first.row[i]):(first.row[i]+n.rows[i]-1), ] <- model$S[ , col.idx]
  }

  # construct and return merged DSM
  res <- list(
    M = M,
    rows = rows.merged,
    cols = cols.merged,
    locked = marginals.inconsistent
  )
  if (have.S) res$S <- S
  class(res) <- c("dsm", "list")

  return(res)
}  
