##
##  Collection of utility functions for DSM matrices
##


## feature weighting with different association scores and optional logarithmic transformation
am.score <- function (AM=c("frequency", "t-score","z-score","Dice"), M, row.freqs, col.freqs, N, sparse=FALSE, log=FALSE) {
  AM <- match.arg(AM)
  f1 <- row.freqs[rownames(M)]    # first marginal frequencies (rows)
  stopifnot(all(!is.na(f1)))
  f2 <- col.freqs[colnames(M)]    # second marginal frequencies (columns)
  stopifnot(all(!is.na(f2)))

  O <- M # observed frequencies
  if (sparse) {
    scores <- matrix(0, nrow(M), ncol(M))
    nonzero <- M > 0
  }

  if (AM == "t-score") {
    E <- outer(f1, f2) / N # expected frequencies
    if (sparse) {
      idx <- nonzero | O > E
      scores[idx] <- (O[idx] - E[idx]) / sqrt(O[idx])
    }
    else {
      scores <- (O - E) / sqrt(O + 1) # "discounted" t-score
    }
  } else if (AM == "z-score") {
    E <- outer(f1, f2) / N # expected frequencies
    if (sparse) {
      idx <- nonzero | O > E
      scores[idx] <- (O[idx] - E[idx]) / sqrt(E[idx])
    }
    else {
      scores <- (O - E) / sqrt(E)
    }
  } else if (AM == "Dice") {
    scores <- 2 * O / outer(f1, f2, "+")
  } else {
    scores <- O # "frequency" measure = observed frequency
  }

  if (log) scores <- sign(scores) * log(abs(scores) + 1) # "signed log" transformation

  return(scores)
}

## p-norms for rows of matrix; use dist() to calculate distance matrix
p.norm <- function (M, p=2) {
  stopifnot(p >= 1)
  if (p == 1) {
    .norm <- function (x) sum(abs(x))
  } else {
    .norm <- function (x) (sum(abs(x)^p))^(1/p)
  }
  return(apply(M, 1, .norm))
}

## cosine similarities between rows of matrix (or different matrices)
cosine <- function (M, M2=M, angles=FALSE, normalised=FALSE) {
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

# standardise all features (columns) with scale() according to p-norm
norm.rows <- function (M, p=2) {
  M / p.norm(M, p)  # faster than sweep(M, 1, p.norm(M, p), "/")
}

# draw arrow with start/end shortened by configurable amount and optional rotated label
draw.arrow <- function (x1, y1, x2, y2, label=NULL, label.pos=0.5, label.offset=5, cut1=0, cut2=0, head1=FALSE, head2=TRUE, length=.15, lwd=1, ...) {
  .dx <- x2 - x1
  .dy <- y2 - y1
  .l <- sqrt(.dx^2 + .dy^2)
  .code <- head2 * 1 + head1 * 2
  arrows(x1 + cut1 * .dx/.l, y1 + cut1 * .dy/.l, x2 - cut2 * .dx/.l, y2 - cut2 * .dy/.l, code=.code, length=length, lwd=lwd, ...)
  if (!missing(label)) {
    .tx <- x1 + label.pos * .dx - label.offset * .dy/.l
    .ty <- y1 + label.pos * .dy + label.offset * .dx/.l
    .angle <- 180 * atan2(.dy, .dx) / pi
    .angle <- ((.angle + 90) %% 180) - 90
    text(.tx, .ty, labels=label, srt=.angle, ...)
  }
}
