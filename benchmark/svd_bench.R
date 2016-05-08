##
## Benchmark different SVD variants
##

library(Matrix)
library(wordspace)
source("benchmark_utils.R")

## data set 1: 4592 x 10000 co-occurrence matrix (BNC dependencies)
load("data/bnc_5k_10k_dsm.rda")
Sparse <- bnc_5k_10k_dsm
Dense <- Sparse
Dense$M <- as.matrix(Dense$M)

## basic information about the data sets
nr <- nrow(Sparse)
nc <- ncol(Sparse)
n.cells <- nr * nc
n.sparse <- nnzero(Sparse$M)

cat("BENCHMARK DATA:\n")
cat(sprintf(" - %d x %d matrix = %.1f M cells\n", nr, nc, n.cells / 1e6))
cat(sprintf(" - %.1f M nonzero entries, fill rate = %.2f%%\n", n.sparse / 1e6, 100 * n.sparse / n.cells))
cat(sprintf(" - matrix storage: %.1f MiB sparse, %.1f MiB dense\n", object.size(Sparse$M) / 2^20, object.size(Dense$M) / 2^20))
cat("\n")


## collect benchmark results in list L
L <- list()
Dense.scored <- dsm.score(Dense, score="t-score", transform="root")
Sparse.scored <- dsm.score(Sparse, score="t-score", transform="root")

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
set.seed(42)
L <- append.list(L, benchmark( rsvd100.dense <<- dsm.projection(Dense.svd, n=100, method="rsvd", oversampling=4, verbose=TRUE), "D rSVD to 100 dim. wordspace", ops ))
L <- append.list(L, benchmark( svd100.sparse <<- dsm.projection(Sparse.svd, n=100, method="svd", verbose=FALSE), "S SVD to 100 dim. wordspace", ops ))
set.seed(42)
L <- append.list(L, benchmark( rsvd100.sparse <<- dsm.projection(Sparse.svd, n=100, method="rsvd", oversampling=4, verbose=TRUE), "S rSVD to 100 dim. wordspace", ops ))

matrix.equal(svd100.dense.R, svd100.dense)
matrix.equal(svd100.dense.R, svd100.sparse, ignore.sign=TRUE, tol=1e-6) # some numerical differences expected
if (FALSE) {
  ## randomized SVD results in some fairly big differences from "correct" solution 
  matrix.equal(svd100.dense.R, rsvd100.dense, ignore.sign=TRUE)
}
matrix.equal(rsvd100.dense, rsvd100.sparse, ignore.sign=TRUE, tol=1e-6) # dense and sparse rSVD reasonably close with same random matrices


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

