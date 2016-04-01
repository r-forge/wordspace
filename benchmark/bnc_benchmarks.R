##
## Benchmark some "wordspace" functions (and alternative implementations) based on BNC data set
##

library(Matrix)
library(wordspace)
source("benchmark_utils.R")

## data set 1: 4592 x 10000 co-occurrence matrix (BNC dependencies)
load("data/bnc_5k_10k_dsm.rda")
Sparse <- bnc_5k_10k_dsm
Dense <- Sparse
Dense$M <- as.matrix(Dense$M)

## data set 2: 4592 x 100 latent SVD dimensions
load("data/bnc_5k_10k_svd100.rda")
Latent100 <- bnc_5k_10k_svd100


## basic information about the data sets
nr <- nrow(Sparse)
stopifnot(nr == nrow(Latent100))
nc <- ncol(Sparse)
nc100 <- ncol(Latent100)
n.cells <- nr * nc
n.sparse <- nnzero(Sparse$M)
n.cells100 <- nr * nc100

cat("BENCHMARK DATA:\n")
cat(sprintf(" - %d x %d matrix = %.1f M cells\n", nr, nc, n.cells / 1e6))
cat(sprintf(" - %.1f M nonzero entries, fill rate = %.2f%%\n", n.sparse / 1e6, 100 * n.sparse / n.cells))
cat(sprintf(" - matrix storage: %.1f MiB sparse, %.1f MiB dense\n", object.size(Sparse$M) / 2^20, object.size(Dense$M) / 2^20))
cat(sprintf(" - SVD-projected %d x %d matrix with %.1f M cells = %.1f MiB\n", nr, nc100, n.cells100 / 1e6, object.size(Latent100) / 2^20))
cat("\n")


## collect benchmark results in list L
L <- list()


## scoring: t-score with sqrt transformation
t.score.root <- function (O, f1, f2, N) {
  E <- outer(f1, f2) / N
  sqrt(pmax( (O - E) / sqrt(O), 0 ))
}
ops <- nr * nc * 7   # number of FP operations for dense matrix

L <- append.list(L, benchmark( M.scored.dense <<- t.score.root(Dense$M, Dense$rows$f, Dense$cols$f, Dense$globals$N), "D sqrt(t-score) R", ops ))
L <- append.list(L, benchmark( Dense.scored <<- dsm.score(Dense, score="t-score", transform="root"), "D sqrt(t-score) wordspace", ops ))
L <- append.list(L, benchmark( Sparse.scored <<- dsm.score(Sparse, score="t-score", transform="root"), "S sqrt(t-score) wordspace", ops ))

matrix.equal(M.scored.dense, Dense.scored$S)
matrix.equal(M.scored.dense, as.matrix(Sparse.scored$S))


## row and column norms
euclidean.norm <- function (M, byrow=TRUE) {
  if (byrow) {
    n2 <- rowSums(M * M)
  } else {
    n2 <- colSums(M * M)
  }
  sqrt(n2)
}
ops <- nr * nc * 2 + max(nr, nc) # approx. number of FP operations for dense matrix

L <- append.list(L, benchmark( row.norms.R <<- euclidean.norm(Dense.scored$S), "D L2-norm rows R", ops ))
L <- append.list(L, benchmark( col.norms.R <<- euclidean.norm(Dense.scored$S, byrow=FALSE), "D L2-norm columns R", ops ))
L <- append.list(L, benchmark( row.norms.R.sparse <<- euclidean.norm(Sparse.scored$S), "S L2-norm rows R", ops ))
L <- append.list(L, benchmark( col.norms.R.sparse <<- euclidean.norm(Sparse.scored$S, byrow=FALSE), "S L2-norm columns R", ops ))

L <- append.list(L, benchmark( row.norms.ws <<- rowNorms(Dense.scored$S, method="euclidean"), "D L2-norm rows wordspace", ops ))
L <- append.list(L, benchmark( col.norms.ws <<- colNorms(Dense.scored$S, method="euclidean"), "D L2-norm columns wordspace", ops ))
L <- append.list(L, benchmark( row.norms.ws.sparse <<- rowNorms(Sparse.scored$S, method="euclidean"), "S L2-norm rows wordspace", ops ))
L <- append.list(L, benchmark( col.norms.ws.sparse <<- colNorms(Sparse.scored$S, method="euclidean"), "S L2-norm columns wordspace", ops ))

vector.equal(row.norms.R, row.norms.R.sparse)
vector.equal(col.norms.R, col.norms.R.sparse)
vector.equal(row.norms.R, row.norms.ws)
vector.equal(col.norms.R, col.norms.ws)
vector.equal(row.norms.R.sparse, row.norms.ws.sparse)
vector.equal(col.norms.R.sparse, col.norms.ws.sparse)


## normalizing rows
normalize.rows.R <- function (M) {
  M / euclidean.norm(M, byrow=TRUE) # clever use of vector recycling
}
ops <- nr * nc * 2 + nr + nr * nc # FP ops = calculation of norms + scaling

L <- append.list(L, benchmark( M.normalized.R <<- normalize.rows.R(Dense.scored$S), "D normalize R", ops ))
L <- append.list(L, benchmark( M.normalized.R.sparse <<- normalize.rows.R(Sparse.scored$S), "S normalize R", ops ))

L <- append.list(L, benchmark( M.normalized.ws <<- normalize.rows(Dense.scored$S), "D normalize wordspace", ops ))
L <- append.list(L, benchmark( M.normalized.ws.sparse <<- normalize.rows(Sparse.scored$S), "S normalize wordspace", ops ))

matrix.equal(M.normalized.R, M.normalized.R.sparse)
matrix.equal(M.normalized.R, M.normalized.ws)
matrix.equal(M.normalized.R.sparse, M.normalized.ws.sparse)


## Manhattan distances
nr1 <- 400
M1.dense <- Dense.scored$S[1:nr1, ]
M1.sparse <- Sparse.scored$S[1:nr1, ]
ops <- (nr1 * (nr1 - 1) / 2) * (nc * 3) # approx. FP ops = diff + abs() + add per vector element over all pairs of row vectors

L <- append.list(L, benchmark( distances.L1.R <<- dist(M1.dense, method="manhattan"), "D L1 distances R", ops ))
L <- append.list(L, benchmark( distances.L1.ws <<- dist.matrix(M1.dense, method="manhattan"), "D L1 distances wordspace", ops ))
L <- append.list(L, benchmark( distances.L1.ws.sparse <<- dist.matrix(M1.sparse, method="manhattan"), "S L1 distances wordspace", ops ))

matrix.equal(as.matrix(distances.L1.R), distances.L1.ws)
matrix.equal(distances.L1.ws, distances.L1.ws.sparse)

if (wordspace.openmp()$available) {
  for (n.threads in c(2, 4, 8)) {
    if (n.threads <= wordspace.openmp()$max) {
      wordspace.openmp(n.threads)
      L <- append.list(L, benchmark( distances.L1.wsOMP <<- dist.matrix(M1.dense, method="manhattan"), sprintf("D L1 distances wordspace (%d threads)", n.threads), ops ))
      L <- append.list(L, benchmark( distances.L1.wsOMP.sparse <<- dist.matrix(M1.sparse, method="manhattan"), sprintf("S L1 distances wordspace (%d threads)", n.threads), ops ))
      matrix.equal(distances.L1.ws, distances.L1.wsOMP)
      matrix.equal(distances.L1.ws, distances.L1.wsOMP.sparse)
    }
  }
  wordspace.openmp(1)
}


## Minkowski distances (p = 3)
ops <- (nr1 * (nr1 - 1) / 2) * (nc * 4) # approx. FP ops = (diff + abs + exp + sum) per vector element over all pairs of row vectors

L <- append.list(L, benchmark( distances.L3.R <<- dist(M1.dense, method="minkowski", p=3), "D L3 distances R", ops ))
L <- append.list(L, benchmark( distances.L3.ws <<- dist.matrix(M1.dense, method="minkowski", p=3), "D L3 distances wordspace", ops ))
L <- append.list(L, benchmark( distances.L3.ws.sparse <<- dist.matrix(M1.sparse, method="minkowski", p=3), "S L3 distances wordspace", ops ))

matrix.equal(as.matrix(distances.L3.R), distances.L3.ws)
matrix.equal(distances.L3.ws, distances.L3.ws.sparse)

if (wordspace.openmp()$available) {
  for (n.threads in c(2, 4, 8)) {
    if (n.threads <= wordspace.openmp()$max) {
      wordspace.openmp(n.threads)
      L <- append.list(L, benchmark( distances.L3.wsOMP.sparse <<- dist.matrix(M1.sparse, method="minkowski", p=3), sprintf("S L3 distances wordspace (%d threads)", n.threads), ops ))
      matrix.equal(distances.L3.ws, distances.L3.wsOMP.sparse)
    }
  }
  wordspace.openmp(1)
}


## cosine similarities (conversion to angular distance has numerical problems in naive implementation)
cosine <- function (M) {
  norms <- euclidean.norm(M, byrow=TRUE)
  tcrossprod(M) / outer(norms, norms)
}
ops <- (nr1 * (nr1 - 1) / 2) * (nc * 2) + (nr1 * nc * 2 + nr1) + (nr1 * nr1) # cross-product, row norms, scale dot products

L <- append.list(L, benchmark( distances.cosine.R <<- cosine(M1.dense), "D cosine distances R", ops ))
L <- append.list(L, benchmark( distances.cosine.R.sparse <<- cosine(M1.sparse), "S cosine distances R", ops ))
L <- append.list(L, benchmark( distances.cosine.ws <<- dist.matrix(M1.dense, method="cosine", normalized=FALSE, convert=FALSE), "D cosine distances wordspace", ops ))
L <- append.list(L, benchmark( distances.cosine.ws.sparse <<- dist.matrix(M1.sparse, method="cosine", normalized=FALSE, convert=FALSE), "S cosine distances wordspace", ops ))

matrix.equal(distances.cosine.R, distances.cosine.R.sparse)
matrix.equal(distances.cosine.R, distances.cosine.ws)
matrix.equal(distances.cosine.R, distances.cosine.ws.sparse)


## Manhattan and cosine in latent space (100 dimensions)
ops <- ops <- (nr * (nr - 1) / 2) * (nc100 * 3)
L <- append.list(L, benchmark( ldistances.L1.R <<- dist(Latent100, method="manhattan"), "SVD-100 L1 distances R", ops ))
L <- append.list(L, benchmark( ldistances.L1.ws <<- dist.matrix(Latent100, method="manhattan"), "SVD-100 L1 distances wordspace", ops ))

matrix.equal(as.matrix(ldistances.L1.R), ldistances.L1.ws)

ops <- (nr * (nr - 1) / 2) * (nc100 * 2) + (nr * nc100 * 2 + nr) + (nr * nr) # cross-product, row norms, scale dot products
L <- append.list(L, benchmark( ldistances.cosine.R <<- cosine(Latent100), "SVD-100 cosine distances R", ops ))
L <- append.list(L, benchmark( ldistances.cosine.ws <<- dist.matrix(Latent100, method="cosine", normalized=FALSE, convert=FALSE), "SVD-100 cosine distances wordspace", ops ))

matrix.equal(ldistances.cosine.R, ldistances.cosine.ws)


## nearest neighbours (on pre-normalized matrix)
nr.nn <- 500
ops <- (nr.nn * nr) * (nc * 2) + nr.nn * nr * 3 # cross-product + scan each distance vector for 20 NN (rough estimate)
Dense.nn <- normalize.rows(Dense.scored$S)
Sparse.nn <- normalize.rows(Sparse.scored$S)
Latent.nn <- normalize.rows(Latent100)
targets.nn <- rownames(Dense.nn)[1:nr.nn]

L <- append.list(L, benchmark( nn.dense <- nearest.neighbours(Dense.nn, term=targets.nn, n=20, normalized=TRUE, convert=TRUE), "D 20 nearest neighbours wordspace", ops ))
L <- append.list(L, benchmark( nn.sparse <- nearest.neighbours(Sparse.nn, term=targets.nn, n=20, normalized=TRUE, convert=TRUE), "S 20 nearest neighbours wordspace", ops ))

ops <- (nr.nn * nr) * (nc100 * 2) + nr.nn * nr * 3 # cross-product + scan each distance vector for 20 NN (rough estimate)
L <- append.list(L, benchmark( nn.latent <- nearest.neighbours(Latent.nn, term=targets.nn, n=20, normalized=TRUE, convert=TRUE), "SVD-100 20 nearest neighbours wordspace", ops ))

ok <- all.equal(nn.dense, nn.sparse)
if (!isTRUE(ok)) print(ok)


## dimensionality reduction with SVD vs. randomized SVD
nr.svd <- 2000
nc.svd <- 5000
svd.dim <- 100
Dense.svd <- Dense.scored$S[1:nr.svd, 1:nc.svd]
Sparse.svd <- Sparse.scored$S[1:nr.svd, 1:nc.svd]

svd.proj <- function (M, n) {
  res <- svd(M, nu=n, nv=0) # don't compute right singular vectors for projection
  res$u %*% diag(res$d[1:n])
}
ops <- 4 * nc.svd * nc.svd * nr.svd + 8 * nc.svd * nr.svd * nr.svd # complexity of full SVD according to Wikipedia talk page (which assumes nr >> nc, so this is calculated for t(M))
ops <- ops + nr.svd * svd.dim # U %*% diag(d) = rescaling of columns

L <- append.list(L, benchmark( svd100.dense.R <<- svd.proj(Dense.svd, 100), "D SVD to 100 dim. R", ops ))
L <- append.list(L, benchmark( svd100.dense <<- dsm.projection(Dense.svd, n=100, method="svd", verbose=FALSE), "D SVD to 100 dim. wordspace", ops ))
set.seed(42) # for reproducible rSVD
L <- append.list(L, benchmark( rsvd100.dense <<- dsm.projection(Dense.svd, n=100, method="rsvd", oversampling=4, verbose=TRUE), "D rSVD to 100 dim. wordspace", ops ))
L <- append.list(L, benchmark( svd100.sparse <<- dsm.projection(Sparse.svd, n=100, method="svd", verbose=FALSE), "S SVD to 100 dim. wordspace", ops ))
set.seed(42) # for reproducible rSVD
L <- append.list(L, benchmark( rsvd100.sparse <<- dsm.projection(Sparse.svd, n=100, method="rsvd", oversampling=4, verbose=TRUE), "S rSVD to 100 dim. wordspace", ops ))

matrix.equal(svd100.dense.R, svd100.dense)
matrix.equal(svd100.dense.R, svd100.sparse, ignore.sign=TRUE, tol=1e-6) # some numerical differences expected
matrix.equal(rsvd100.dense, rsvd100.sparse, ignore.sign=TRUE, tol=1e-6) # dense and sparse rSVD reasonably close with same random matrices
if (FALSE) {
  ## randomized SVD results in some fairly big differences from "correct" solution 
  matrix.equal(svd100.dense.R, rsvd100.dense, ignore.sign=TRUE)
}


## dimensionality reduction with random indexing
ri.dim <- 2000
ri.rate <- .05
## FP op count assumes semi-sparse implementation using only non-zero entries of random projection matrix Q 
ops <- 2 * (ri.dim * nc * ri.rate) + (nr * ri.dim * nc * ri.rate) + (nr * ri.dim) # generate Q, matrix multiplication, rescale random dimensions

set.seed(42)
L <- append.list(L, benchmark( ri2k.dense <<- dsm.projection(Dense.scored$S, n=ri.dim, method="ri", rate=ri.rate, verbose=TRUE), "D RI to 2000 dim. wordspace", ops ))
set.seed(42)
L <- append.list(L, benchmark( ri2k.sparse <<- dsm.projection(Sparse.scored$S, n=ri.dim, method="ri", rate=ri.rate, verbose=TRUE), "S RI to 2000 dim. wordspace", ops ))
## cannot compare results because random numbers are generated in different ways by R code and C code



## print summary report
cat("\n")
cat(sprintf("===== Benchmark results for wordspace v%s (%s) =====\n", packageVersion("wordspace"), version$version.string))

print(do.call(rbind, L))

