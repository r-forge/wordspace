context.vectors <- function (M, contexts, split="\\s+", drop.missing=TRUE) {
  tokens.list <- strsplit(contexts, split, perl=TRUE)
  known.terms <- rownames(M)
  nC <- ncol(M)

  CM <- t(vapply(tokens.list, function (tokens) {
    idx <- na.omit(match(tokens, known.terms)) # row numbers of known terms in M (possibly repeated))
    if (length(idx) > 0) {
      colMeans(M[idx, , drop=FALSE]) # centroid vector
    } else {
      if (drop.missing) rep(NA, nC) else rep(0, nC) # return null vector for context without known tokens
    }
  }, FUN.VALUE=numeric(nC), USE.NAMES=FALSE))
  rownames(CM) <- if (is.null(names(contexts))) 1:length(contexts) else names(contexts)
  if (!is.null(colnames(M))) colnames(CM) <- colnames(M)
  if (drop.missing) {
    idx.miss <- is.na(CM[,1]) # assuming there were no NAs or NaNs in M
    CM[!idx.miss, ]
  } else {
    CM
  }
}
