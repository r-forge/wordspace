## Validate dist.matrix() implementation against pure R code using dist() and matrix operations
## based on small "hieroglyphs" example matrix

library(wordspace)
library(Matrix)

matrix.equal <- function(x, y, name="matrix comparison", tol=1e-10, verbose=TRUE) {
  if (any(is.na(x))) {
    cat(sprintf("%s: missing values in first matrix\n", name))
    invisible(FALSE)
  } else if (any(is.na(y))) {
    cat(sprintf("%s: missing values in second matrix\n", name))
    invisible(FALSE)
  } else if (nrow(x) != nrow(y) || ncol(x) != ncol(y)) {
    if (verbose) cat(sprintf("%s: different matrix formats %d x %d != %d x %d\n", name, nrow(x), ncol(x), nrow(y), ncol(y)))
    invisible(FALSE)
  } else {
    max.diff <- max(abs(x - y))
    if (max.diff < tol) {
      invisible(TRUE)
    } else {
      if (verbose) cat(sprintf("%s: largest difference between matrices = %g exceeds tolerance limit\n", name, max.diff))
      invisible(FALSE)
    }
  }
}

## compute distance matrix based on R function
## d.fnc must compute distances between corresponding rows of two matrices
dist.matrix.R <- function (M, d.fnc, ...) {
  M <- as.matrix(M)
  idx <- seq_len(nrow(M))
  outer(idx, idx, function (i, j) d.fnc(M[i, , drop=FALSE], M[j, , drop=FALSE], ...))
}


## test data: sparse and dense co-occurrence matrix
O <- DSM_HieroglyphsMatrix
E <- outer(rowSums(O), colSums(O)) / sum(O)
M1 <- (O - E) / sqrt(E)               # z-score (signed), dense matrix
M2 <- as(pmax(M1, 0), "sparseMatrix") # sparse z-score, sparse matrix


## Euclidean distance
M1.l2.ws <- dist.matrix(M1, method="euclidean")
M2.l2.ws <- dist.matrix(M2, method="euclidean")

M1.l2.R <- dist(M1, method="euclidean") 
M2.l2.R <- dist(M2, method="euclidean") # converts to dense matrix

stopifnot(matrix.equal(M1.l2.ws, as.matrix(M1.l2.R)))
stopifnot(matrix.equal(M2.l2.ws, as.matrix(M2.l2.R)))
stopifnot(all(diag(M1.l2.ws) == 0)) # must be exact zeroes
stopifnot(all(diag(M2.l2.ws) == 0))


## Manhattan distance
M1.l1.ws <- dist.matrix(M1, method="manhattan")
M2.l1.ws <- dist.matrix(M2, method="manhattan")

M1.l1.R <- dist(M1, method="manhattan") 
M2.l1.R <- dist(M2, method="manhattan") # converts to dense matrix

stopifnot(matrix.equal(M1.l1.ws, as.matrix(M1.l1.R)))
stopifnot(matrix.equal(M2.l1.ws, as.matrix(M2.l1.R)))
stopifnot(all(diag(M1.l1.ws) == 0)) # must be exact zeroes
stopifnot(all(diag(M2.l1.ws) == 0))


## maximum distance
M1.linf.ws <- dist.matrix(M1, method="maximum")
M2.linf.ws <- dist.matrix(M2, method="maximum")

M1.linf.R <- dist(M1, method="maximum") 
M2.linf.R <- dist(M2, method="maximum") # converts to dense matrix

stopifnot(matrix.equal(M1.linf.ws, as.matrix(M1.linf.R)))
stopifnot(matrix.equal(M2.linf.ws, as.matrix(M2.linf.R)))
stopifnot(all(diag(M1.linf.ws) == 0)) # must be exact zeroes
stopifnot(all(diag(M2.linf.ws) == 0))


## Minkowski distance
for (p in c(.1, .25, .5, 1, 1.5, 2, 2.5, 3, 5, 10)) {
  M1.lp.ws <- dist.matrix(M1, method="minkowski", p=p)
  M2.lp.ws <- dist.matrix(M2, method="minkowski", p=p)

  if (p < 1) {
    ## dist.matrix() uses non-standard definition of Minkowski distance for p < 1
    M1.lp.ws <- M1.lp.ws ^ (1/p)
    M2.lp.ws <- M2.lp.ws ^ (1/p)
  }
  
  M1.lp.R <- dist(M1, method="minkowski", p=p) 
  M2.lp.R <- dist(M2, method="minkowski", p=p) # converts to dense matrix

  stopifnot(matrix.equal(M1.lp.ws, as.matrix(M1.lp.R), tol=1e-6)) # dist() seems to use different alogrithm
  stopifnot(matrix.equal(M2.lp.ws, as.matrix(M2.lp.R), tol=1e-6)) # resulting in larger errors on 32-bit i386 
  stopifnot(all(diag(M1.lp.ws) == 0)) # must be exact zeroes
  stopifnot(all(diag(M2.lp.ws) == 0))

  if (p == 1) {
    stopifnot(matrix.equal(M1.lp.ws, as.matrix(M1.l1.ws)))
    stopifnot(matrix.equal(M2.lp.ws, as.matrix(M2.l1.ws)))
  }
  if (p == 2) {
    stopifnot(matrix.equal(M1.lp.ws, as.matrix(M1.l2.ws)))
    stopifnot(matrix.equal(M2.lp.ws, as.matrix(M2.l2.ws)))
  }
}


## Canberra distance
##  - wordspace package implements different form of Canberra distance than dist()
##  - version in dist.matrix() is more sensible for vectors that have both positive and negative entries
##  - dist() and dist.matrix() also disagree on matching zeroes (which dist() imputes as missing values)

M1a <- abs(M1)  # so all entries are positive
M2a <- M2 + 1   # ensure that no two row vectors have matching zeroes
diag(M2a) <- 0
M2a <- as(M2a, "dgCMatrix") # convert back to sparse matrix

M1a.cb.ws <- dist.matrix(M1a, method="canberra")
M2a.cb.ws <- dist.matrix(M2a, method="canberra")

M1a.cb.R <- dist(M1a, method="canberra") 
M2a.cb.R <- dist(M2a, method="canberra") # converts to dense matrix

stopifnot(matrix.equal(M1a.cb.ws, as.matrix(M1a.cb.R)))
stopifnot(matrix.equal(M2a.cb.ws, as.matrix(M2a.cb.R)))
stopifnot(all(diag(M1a.cb.ws) == 0)) # must be exact zeroes
stopifnot(all(diag(M2a.cb.ws) == 0))


## Jaccard metric and overlap measure
Jaccard.fnc <- function (x, y, convert=FALSE) {
  stopifnot(all(x >= 0) && all(y >= 0))
  res <- rowSums(pmin(x, y)) / rowSums(pmax(x, y))
  if (convert) 1 - res else res
}
overlap.fnc <- function (x, y, normalized=FALSE, convert=FALSE) {
  stopifnot(all(x >= 0) && all(y >= 0))
  res <- rowSums(pmin(x, y))
  if (!normalized) res <- res / rowSums(x)
  if (convert) 1 - res else res
}

M1a.jc.ws <- dist.matrix(M1a, method="jaccard", convert=FALSE)
M2.jc.ws <- dist.matrix(M2, method="jaccard", convert=FALSE)
M1a.ol.ws <- dist.matrix(M1a, method="overlap", convert=FALSE)
M2.ol.ws <- dist.matrix(M2, method="overlap", convert=FALSE)

M1a.jc.R <- dist.matrix.R(M1a, Jaccard.fnc, convert=FALSE)
M2.jc.R <- dist.matrix.R(M2, Jaccard.fnc, convert=FALSE)
M1a.ol.R <- dist.matrix.R(M1a, overlap.fnc, convert=FALSE)
M2.ol.R <- dist.matrix.R(M2, overlap.fnc, convert=FALSE)

stopifnot(matrix.equal(M1a.jc.ws, M1a.jc.R))
stopifnot(matrix.equal(M2.jc.ws, M2.jc.R))
stopifnot(matrix.equal(M1a.ol.ws, M1a.ol.R))
stopifnot(matrix.equal(M2.ol.ws, M2.ol.R))

M1a.jc.ws <- dist.matrix(M1a, method="jaccard", convert=TRUE) # Jaccard metric
M2.jc.ws <- dist.matrix(M2, method="jaccard", convert=TRUE)

M1a.jc.R <- dist.matrix.R(M1a, Jaccard.fnc, convert=TRUE)
M2.jc.R <- dist.matrix.R(M2, Jaccard.fnc, convert=TRUE)

stopifnot(matrix.equal(M1a.jc.ws, M1a.jc.R))
stopifnot(matrix.equal(M2.jc.ws, M2.jc.R))

M3a <- rbind(c(0, 2, 1, 0, 1), # test against known results (similarity measures only)
             c(4, 1, 3, 0, 0))
M3a.jc.gold <- rbind(c(   1, 2/10), # gold = expected similarity matrix
                     c(2/10,    1))
M3a.ol.gold <- rbind(c(  1, 2/4),
                     c(2/8,   1))
M3a.olraw.gold <- rbind(c(4, 2),    # unnormalized overlap
                        c(2, 8))

M3b <- rbind(c(0, 0, 0),
             c(0, 0, 0),
             c(1, 1, 0))
M3b.jc.gold <- rbind(c(1, 1, 0),    # J(0, 0) = 1, but J(0, x) = 0
                     c(1, 1, 0),
                     c(0, 0, 1))
M3b.ol.gold <- rbind(c(1, 1, 1),    # o(0, x) = 1, but o(x, 0) = 0
                     c(1, 1, 1),
                     c(0, 0, 1))
M3b.olraw.gold <- rbind(c(0, 0, 0), # raw overlap or(0, x) = or(x, 0) = 0
                        c(0, 0, 0),
                        c(0, 0, 2))

M3a.jc.ws <- dist.matrix(M3a, method="jaccard", convert=FALSE)
M3a.ol.ws <- dist.matrix(M3a, method="overlap", convert=FALSE)
M3a.olraw.ws <- dist.matrix(M3a, method="overlap", convert=FALSE, normalized=TRUE)

stopifnot(matrix.equal(M3a.jc.ws,    M3a.jc.gold))
stopifnot(matrix.equal(M3a.ol.ws,    M3a.ol.gold))
stopifnot(matrix.equal(M3a.olraw.ws, M3a.olraw.gold))

M3b.jc.ws <- dist.matrix(M3b, method="jaccard", convert=FALSE)
M3b.ol.ws <- dist.matrix(M3b, method="overlap", convert=FALSE)
M3b.olraw.ws <- dist.matrix(M3b, method="overlap", convert=FALSE, normalized=TRUE)

stopifnot(matrix.equal(M3b.jc.ws,    M3b.jc.gold))
stopifnot(matrix.equal(M3b.ol.ws,    M3b.ol.gold))
stopifnot(matrix.equal(M3b.olraw.ws, M3b.olraw.gold))


## cosine similarity and angular distance
M1n <- normalize.rows(M1)
M2n <- normalize.rows(M2)

M1.cos.ws <- dist.matrix(M1, method="cosine", convert=FALSE, normalized=FALSE)
M1n.cos.ws <- dist.matrix(M1n, method="cosine", convert=FALSE, normalized=TRUE)

M2.cos.ws <- dist.matrix(M2, method="cosine", convert=FALSE, normalized=FALSE)
M2n.cos.ws <- dist.matrix(M2n, method="cosine", convert=FALSE, normalized=TRUE)

M1n.cos.R <- tcrossprod(M1n)
M2n.cos.R <- as.matrix(tcrossprod(M2n))

stopifnot(matrix.equal(M1.cos.ws, M1n.cos.R))
stopifnot(matrix.equal(M1n.cos.ws, M1n.cos.R))
stopifnot(matrix.equal(M2.cos.ws, M2n.cos.R))
stopifnot(matrix.equal(M2n.cos.ws, M2n.cos.R))

tmp <- pmin(M1n.cos.R, 1)
tmp[tmp > 1-(1e-12)] <- 1 # clamping function of dist.matrix()
M1.deg.R <- acos(tmp) * 180 / pi
tmp <- pmin(M2n.cos.R, 1)
tmp[tmp > 1-(1e-12)] <- 1 # clamping function of dist.matrix()
M2.deg.R <- acos(tmp) * 180 / pi

M1.deg.ws <- dist.matrix(M1, method="cosine", convert=TRUE, normalized=FALSE)
M2.deg.ws <- dist.matrix(M2, method="cosine", convert=TRUE, normalized=FALSE)

stopifnot(matrix.equal(M1.deg.ws, M1.deg.R))
stopifnot(matrix.equal(M2.deg.ws, M2.deg.R))
stopifnot(all(diag(M1.deg.ws) == 0)) # clamping effect guarantees exact zeroes
stopifnot(all(diag(M2.deg.ws) == 0))
