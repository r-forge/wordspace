context.vectors <- function (M, contexts, split="\\s+", drop.missing=TRUE, row.names=NULL) {
  M <- find.canonical.matrix(M) # ensure that M is a suitable matrix, or extract matrix from DSM
  known.terms <- rownames(M)
  nR <- nrow(M)
  nC <- ncol(M)
  if (is.null(row.names)) {
    row.names <- if (is.null(names(contexts))) 1:length(contexts) else names(contexts)
  } else {
    if (length(row.names) != length(contexts)) stop("row.names= must have same length as contexts=")
  }
  if (is.character(contexts)) {
    tokens.list <- strsplit(contexts, split, perl=TRUE)
  } else {
    if (!is.list(contexts)) stop("contexts= must be either a character vector or a list of index vectors")
    tokens.list <- contexts
  }

  CM <- t(vapply(tokens.list, function (tokens) {
    if (is.character(tokens)) {
      idx <- na.omit(match(tokens, known.terms)) # row numbers of known terms in M (possibly repeated))
    } else {
      idx <- tokens # should be a valid index vector into M (either row numbers or logical)
      if (is.logical(idx)) idx <- which(idx)
      if (length(idx) > 0 && (min(idx) < 1 || max(idx) > nR)) stop("invalid index vector in contexts= (row number out of range)")
    }
    if (length(idx) > 0) {
      colMeans(M[idx, , drop=FALSE]) # centroid vector
    } else {
      if (drop.missing) rep(NA, nC) else rep(0, nC) # return null vector for context without known tokens
    }
  }, FUN.VALUE=numeric(nC), USE.NAMES=FALSE))
  rownames(CM) <- row.names
  if (!is.null(colnames(M))) colnames(CM) <- colnames(M)
  if (drop.missing) {
    idx.miss <- is.na(CM[, 1]) # assuming there were no NAs or NaNs in M
    CM[!idx.miss, ]
  } else {
    CM
  }
}
