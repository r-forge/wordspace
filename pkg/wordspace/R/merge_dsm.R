merge.dsm <- function (x, y, ..., rows=TRUE, all=FALSE) {
  models <- list(x, y, ...)
  models.info <- lapply(models, check.dsm, validate=TRUE) # validate models and extract dimensions etc.

  if (all) stop("all=TRUE is not yet implemented")
  if (!rows) stop("rows=FALSE is not yet implemented")
  
  ncol.vec <- sapply(models.info, function (i) i$ncol) # extract column counts
  if (min(ncol.vec) != max(ncol.vec)) stop("all DSM objects must have the same number of columns (i.e. dimensions)")
  have.S.vec <- sapply(models.info, function (i) i$have.S) # are scored matrices available
  have.S <- all(have.S.vec)
  if (any(have.S.vec) && !have.S.vec) warning("some but not all DSM objects contain score matrix S, dropped from result")
  
  col.names <- models[[1]]$cols$term # column labels must be the same for all models
  colnames.ok.vec <- sapply(models, function (m) m$cols$term == col.names)
  if (!all(colnames.ok.vec)) stop("all DSM objects must have the same features (terms / context labels)")
  row.vars <- colnames(models[[1]]$rows) # row information tables must be compatible
  rows.ok.vec <- sapply(models, function (m) colnames(m$rows) == row.vars)
  if (!all(rows.ok.vec)) stop("row information tables of all DSM objects must be compatible")
  
  res <- list(
    M = do.call(rbind, lapply(models, function (m) m$M)),
    rows = do.call(rbind, lapply(models, function (m) m$rows)),
    cols = models[[1]]$cols
  )
  if (have.S) {
    res$S <- do.call(rbind, lapply(models, function (m) m$S))
  }
  class(res) <- c("dsm", "list")

  return(res)
}
