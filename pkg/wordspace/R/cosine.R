cosine <-
function (M, M2=M, angles=FALSE, normalised=FALSE) {
  sim <- M %*% t(M2) # should calculate dot products between rows of M and rows of M2
  if (!normalised) {
    norms.M <- p.norm(M, 2) # norms of row vectors (if not normalised yet)
    norms.M2 <- p.norm(M2, 2)
    sim <- sim / outer(norms.M, norms.M2)
  }
  if (angles) {
    stopifnot(all(sim >= -(1+1e-6) & sim <= 1+1e-6)) # check that cosines are in appropriate range
    sim <- pmin(pmax(sim, -1), 1) # clamp to range [-1, 1]
    sim <- acos(sim) / pi * 180 # angles are returned in degrees
  }
  rownames(sim) <- rownames(M)
  colnames(sim) <- rownames(M2)
  return(sim)
}

