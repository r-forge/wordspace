##
##  Benchmark dense and sparse matrix operations for vector distance / smilarity calculations
##

library(Matrix) # sparse matrix library

load("data/bnc_vobj_basic.rda") # sample data set: 1678 x 4687 with 349378 nonzero cells (4.4% fill rate)

M <- bnc.vobj.basic$M # operate on log-transformed frequency counts
M <- log10(M + 1)

## M <- M[seq(1,nrow(M),5), ] # -- reduce M for test runs

nr <- as.double(nrow(M))
nc <- as.double(ncol(M))

M1 <- M[seq(1, nr, 2), seq(1, nc, 2)] # smaller matrix for slow dist() and svd() benchmarks: 839 x 2343 cells
nr1 <- as.double(nrow(M1))
nc1 <- as.double(ncol(M1))

cat("BENCHMARK DATA:\n")
print(data.frame(rows=c(nr,nr1), cols=c(nc,nc1), cells=c(nr*nc,nr1*nc1), filled=c(sum(M > 0), sum(M1 > 0)), row.names=c("full (M)", "small (M1)")))
cat("\n")

benchmark <- function (expr, name, n.ops=NULL, silent=FALSE) {
  gc()
  .time <- system.time(expr)
  .elapsed <- .time[3]
  .user <- .time[1]
  .mops <- if (missing(n.ops)) NA else (n.ops / 1e6) / .elapsed
  .res <- data.frame(Time=.elapsed, MOPS=round(.mops, 2), CPUTime=.user, row.names=name)
  if (!silent) print(.res)
  .res
}

matrix.equal <- function (M1, M2, name="matrix comparison", tol=1e-12) {
  .error <- max(abs(M1 - M2))
  if (.error > tol) cat(sprintf("%s: max error of %g exceeds tolerance limit\n", name, .error))
}

## normalise vectors to unit length (Euclidean norm)
normalise <- function (M) {
  .len <- sqrt(rowSums(M * M))
  M / .len  # exploits recycling of vectors and column-oriented storage of matrices
}

## cosine measure with/without pre-normalisation (using faster tcrossprod)
cosine <- function (M, normalise=TRUE) {
  if (normalise) {
    M <- normalise(M)
    tcrossprod(M)
  } else {
    .len <- sqrt(rowSums(M * M))
    tcrossprod(M) / outer(.len, .len)
  }
}

## calculate euclidean distances using matrix operations (even though this may be numerically unstable)
## based on || x - y ||^2 = ||x||^2 + ||y||^2 - 2 * <x, y>
euclid <- function (M) {
  .len2 <- rowSums(M * M)
  .dist2 <- outer(.len2, .len2, FUN="+") - 2 * as.matrix(tcrossprod(M)) # convert sparse Matrix to regular dense matrix if necessary
  sqrt(ifelse(.dist2 > 0, .dist2, 0))
}

## full or truncated singular value decomposition (dense matrix only)
svd.wrapper <- function (M, dim=NA, projection=FALSE) {
  n <- min(dim(M))
  if (!missing(dim)) n <- max(dim, n)
  nu <- n
  nv <- if (projection) 0 else n
  svd(M, nu=nu, nv=nv)
}

## ---------- dense matrix operations ----------
cat("===== Running benchmarks for DENSE matrices =====\n")
results.dense <-
  list(benchmark(Mt <<- t(M), "transpose D", nr * nc),
       benchmark(M.norm <<- normalise(M), "normalise D", nr * nc),
       benchmark(d1.dist <<- dist(M1, method="euclidean"), "dist() D", nr1 * nr1 * nc1),
       benchmark(d.inner <<- M %*% Mt, "inner M %*% t(M) D", nr * nr * nc),
       benchmark(d.tcrossprod <<- tcrossprod(M), "inner tcrossprod D", nr * nr * nc),
       benchmark(d.crossprod <<- crossprod(Mt), "inner crossprod t(M) D", nr * nr * nc),
       benchmark(d.cosine1 <<- cosine(M, normalise=TRUE), "cosine normalised D", nr * nr * nc + 2 * nr * nc),
       benchmark(d.cosine2 <<- cosine(M, normalise=FALSE), "cosine general D", nr * nr * nc + nr * nc + nr * nr),
       benchmark(d.euclid <<- euclid(M), "euclid() D", nr * nr * nc + nr * nc + nr * nr),
       benchmark(d1.euclid <<- euclid(M1), "euclid() small D", nr1 * nr1 * nc1 + nr1 * nc1 + nr1 * nr1), # for validation against dist()
       benchmark(svd.full <<- svd.wrapper(M1), "SVD full D", 4 * nc1*nc1*nr1 + 8 * nc1*nr1*nr1 + 9 * nr1*nr1*nr1), # complexity according to a Wikipedia talk page ... (actually for t(M), since Wikipedia assumes nr >> nc)
       benchmark(svd.trunc <<- svd.wrapper(M1, 42), "SVD truncated D", 4 * nc1*nc1*nr1 + 8 * nc1*nr1*nr1 + 9 * nr1*nr1*nr1), # measure speed up wrt. full SVD, i.e. assume operation count of full SVD
       benchmark(svd.proj <<- svd.wrapper(M1, 42, TRUE), "SVD projection D", 4 * nc1*nc1*nr1 + 8 * nr1*nr1*nr1) # projection omits computation of V
       )

cat("----- validating results -----\n")
matrix.equal(d.inner, d.tcrossprod, "M %*% t(M) == tcrossprod(M)")
matrix.equal(d.crossprod, d.tcrossprod, "crossprod(M) == tcrossprod(M)")
matrix.equal(d.cosine1, d.cosine2, "cosine normalised == cosine general")
matrix.equal(as.matrix(d1.dist), d1.euclid, "dist() == euclid()", tol=1e-6) # fast matrix algorithm is very inaccurate
matrix.equal(svd.full$u %*% (svd.full$d * t(svd.full$v)), M1, "U * D * t(V) == M1")


## ---------- sparse matrix operations ----------

cat("===== Running benchmarks for SPARSE matrices =====\n")
results.sparse <-
  list(benchmark(SM <<- Matrix(M, sparse=TRUE), "construct S", nr * nc),
       benchmark(SM1 <<- Matrix(M1, sparse=TRUE), "construct small S", nr1 * nc1),
       benchmark(SMt <<- t(SM), "transpose S", nr * nc),
       benchmark(SM.norm <<- normalise(SM), "normalise S", nr * nc),
       benchmark(Sd1.dist <<- dist(SM1, method="euclidean"), "dist() S", nr1 * nr1 * nc1),
       benchmark(Sd.inner <<- SM %*% SMt, "inner M %*% t(M) S", nr * nr * nc),
       benchmark(Sd.tcrossprod <<- tcrossprod(SM), "inner tcrossprod S", nr * nr * nc),
       benchmark(Sd.crossprod <<- crossprod(SMt), "inner crossprod t(M) S", nr * nr * nc),
       benchmark(Sd.cosine1 <<- cosine(SM, normalise=TRUE), "cosine normalised S", nr * nr * nc + 2 * nr * nc),
       benchmark(Sd.cosine2 <<- cosine(SM, normalise=FALSE), "cosine general S", nr * nr * nc + nr * nc + nr * nr),
       benchmark(Sd.euclid <<- euclid(SM), "euclid() S", nr * nr * nc + nr * nc + nr * nr),
       benchmark(Sd1.euclid <<- euclid(SM1), "euclid() small S", nr1 * nr1 * nc1 + nr1 * nc1 + nr1 * nr1) # for validation against dist()
       )

cat("----- validating results -----\n")
matrix.equal(Sd.inner, Sd.tcrossprod, "M %*% t(M) == tcrossprod(M)")
matrix.equal(Sd.crossprod, Sd.tcrossprod, "crossprod(M) == tcrossprod(M)")
matrix.equal(Sd.cosine1, Sd.cosine2, "cosine normalised == cosine general")
matrix.equal(as.matrix(Sd1.dist), Sd1.euclid, "dist() == euclid()", tol=1e-6) # fast matrix algorithm is very inaccurate


## ---------- result table ----------

cat("===== Benchmark results =====\n")

print(do.call(rbind, c(results.dense, results.sparse)))

