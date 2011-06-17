cosine <- function (M, M2=M, angles=FALSE, normalized=FALSE) {
  # tcrossprod(M, M2) == M %*% t(M2) calculates dot products between rows of M and rows of M2
  sim <- if (missing(M2)) tcrossprod(M) else tcrossprod(M, M2)
  # need to coerce to regular matrix if sparse (for pmin/pmax, but generally more efficient)
  sim <- as.matrix(sim) # ensure this is a dense matrix
  if (!normalized) {
    norms.M <- rowNorms(M, "euclidean") # norms of row vectors (if not normalised yet)
    norms.M2 <- if (missing(M2)) norms.M else rowNorms(M2, "euclidean")
    sim <- sim / outer(norms.M, norms.M2)
  }

  if (angles) {
    stopifnot(all(sim >= -(1+1e-6) & sim <= 1+1e-6)) # check that cosines are in appropriate range
    sim[sim < -1] <- -1           # clamp to range [-1, 1]
    sim[sim > 1] <- 1             # (pmin/pmax eat many GiB of memory for Matrix class??)
    sim <- acos(sim) / pi * 180   # angles are returned in degrees
  }

  rownames(sim) <- rownames(M)
  colnames(sim) <- rownames(M2)
  return(sim)
}

