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

  rownames(scores) <- rownames(M) # make sure that row and column names are preserved
  colnames(scores) <- colnames(M)
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

# SVD decomposition to specified number of dimensions
svd.decomp <- function (M, n=min(dim(M))) {
	stopifnot(n <= min(dim(M)))
	SVD <- svd(M, nu=n, nv=n)
	M.svd <- SVD$u %*% diag(SVD$d[1:n], n, n)
	M.svd <- norm.rows(M.svd)  # renormalise rows (for cosine calculation later on)
	rownames(M.svd) <- rownames(M)
	R2 <- sum(SVD$d[1:n]^2) / sum(SVD$d^2) # proportion of variance "explained" by SVD dimensions
	list(d=SVD$d, u=SVD$u, v=SVD$v, M=M.svd, R2=R2)
}

# project matrix into SVD subspace (checks bounds & works correctly for 1-dimensional projection)
svd.projection <- function (SVD, dims=length(SVD$d), labels=NULL) {
	n <- length(SVD$d)
	stopifnot(dims >= 1 & dims <= n)
	d <- SVD$d[1:dims]
	M <- SVD$u[, 1:dims, drop=FALSE] %*% diag(d, nrow=dims)
	if (!missing(labels)) rownames(M) <- labels
	return(M)
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

## -- TODO: general arguments for these functions --
# find nearest neighbours of specified noun (lemma), with optional visualisation
library(MASS)
neighbours <- function (M, term, n=10, plot=c("none", "pca", "mds", "sammon"), n.plot=n, aspect=1, iterate=FALSE, normalised=FALSE, edges=FALSE, edges.range=NULL, edges.lwd=6, edges.col="#AABBFF") {
  plot <- match.arg(plot)
  auto.edgerange <- missing(edges.range) # so range is automatically recomputed in the iterative plot
  if (! term %in% rownames(M)) stop("Noun '", term, "' is not in database.")
  term.vec <- M[term,,drop=FALSE] # make sure vector is interpreted as single-row matrix
  d <- cosine(M, term.vec, normalised=normalised, angles=TRUE)
  idx <- order(d)
  idx.show <- idx[2:(n+1)]  # first "neighbour" is term itself at distance 0
  if (plot != "none") {
    if (n.plot < 5) iterate <- FALSE
    iterate.n <- if (iterate) 5:n.plot else n.plot
    coords <- NULL
    for (n.plot in iterate.n) {
      idx.plot <- idx[1:(n.plot+1)] # map term vector and n.plot nearest neighbours into 2D plane
      N <- M[idx.plot, ]
      N.d <- cosine(N, normalised=normalised, angles=TRUE) # distance matrix (for MDS and similarity graph)
      if (plot == "pca") {
        P <- prcomp(N, center=TRUE, scale=FALSE)
        coords <- P$x[, 1:2]
      } else {
        if (iterate && !is.null(coords)) {
          # add dummy entry for new data point, then shake well
          coords <- jitter(rbind(coords, c(0,0)), amount=10) 
          coords <- sammon(N.d, y=coords, k=2, trace=FALSE)$points
        } else {
          if (plot == "sammon" || iterate) {
            coords <- sammon(N.d, k=2, trace=FALSE)$points
          }
          else {
            coords <- isoMDS(N.d, k=2, p=2, trace=FALSE)$points
          }
        }
      }
      labels <- rownames(M)[idx.plot]
      x.range <- extendrange(coords[, 1], f=.05)
      y.range <- extendrange(coords[, 2], f=.05)
      .asp <- diff(x.range) / diff(y.range)
      if (.asp < aspect) {
        x.range <- extendrange(x.range, f=(aspect/.asp-1)/2)
      } else if (.asp > aspect) {
        y.range <- extendrange(y.range, f=(.asp/aspect-1)/2)
      }
      if (iterate) Sys.sleep(.2)
      plot(coords, pch=20, xlim=x.range, ylim=y.range, xlab="", ylab="")
      if (edges) {
        ## draw connecting lines between points
        .midx <- t(combn(1:nrow(N), 2)) # 2d-index for upper triangle of distance matrix N.d
        .P1 <- .midx[,1] # idx of start node of edge
        .P2 <- .midx[,2] # idx of end node of edge
        .angle <- N.d[.midx] # angle between P1 and P2
        if (auto.edgerange) {
          edges.range <- c(mean(.angle) - 3 * sd(.angle), mean(.angle)) # show half of edges
          cat(sprintf("[showing edges for range %.1f ... %.1f degrees]\n", edges.range[1], edges.range[2]))
        }
        .lwd <- edges.lwd * (edges.range[2] - .angle) / diff(edges.range) # go from 0 to max. lwd in specified range
        .lwd <- pmin(.lwd, edges.lwd)
        .draw <- .lwd > 0.1 # draw only edges with minimum lwd of 0.1
        segments(coords[.P1[.draw],1], coords[.P1[.draw],2], coords[.P2[.draw],1], coords[.P2[.draw],2], lwd=.lwd[.draw], col=edges.col)
        points(coords, pch=20) # redraw points covered by edges
      }
      points(coords[1,1], coords[1,2], pch=20, cex=1.2, col="red")
      text(coords[-1, ], labels=labels[-1], pos=3, font=2)
      text(coords[1,1], coords[1,2], labels=labels[1], pos=3, font=2, cex=1.2, col="red")
    }
  }
  round(d[idx.show, ], 2) 
}

