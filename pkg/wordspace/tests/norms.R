## Validate various matrix norms and normalization functions.

library(wordspace)
library(Matrix)

almost.equal <- function (x, y, name="vector comparison", tol=1e-12) {
  if (length(y) == 1) y <- rep(y, length(x)) # comparsion of vector against scalar
  if (length(x) != length(y)) {
    cat(sprintf("%s: vector lengths differ\n", name))
    invisible(FALSE)
  } else if (anyNA(x)) {
    cat(sprintf("%s: missing values in first vector\n", name))
    invisible(FALSE)
  } else if (anyNA(y)) {
    cat(sprintf("%s: missing values in second vector\n", name))
    invisible(FALSE)
  } else {
    d <- max(abs(x - y))
    if (d > tol) {
      cat(sprintf("%s: largest difference = %g exceeds tolerance\n", name, d))
      invisible(FALSE)
    } else {
      invisible(TRUE)
    }
  }
}

## test data: small sparse and dense example matrix (unscaled)
M1 <- DSM_TermTermMatrix    # dense
M2 <- DSM_TermContextMatrix # sparse

## Euclidean, Manhattan and maximum norm, as well as Hamming length
tmp <- M1^2
stopifnot(almost.equal(rowNorms(M1, "euclidean"), sqrt(rowSums(tmp)), name="rowNorms, dense, Euclidean"))
stopifnot(almost.equal(colNorms(M1, "euclidean"), sqrt(colSums(tmp)), name="colNorms, dense, Euclidean"))

tmp <- abs(M1)
stopifnot(almost.equal(rowNorms(M1, "manhattan"), rowSums(tmp), name="rowNorms, dense, Manhattan"))
stopifnot(almost.equal(colNorms(M1, "manhattan"), colSums(tmp), name="colNorms, dense, Manhattan"))

stopifnot(almost.equal(rowNorms(M1, "maximum"), apply(tmp, 1, max), name="rowNorms, dense, maximum"))
stopifnot(almost.equal(colNorms(M1, "maximum"), apply(tmp, 2, max), name="colNorms, dense, maximum"))
stopifnot(almost.equal(rowNorms(M1, "minkowski", p=Inf), apply(tmp, 1, max), name="rowNorms, dense, p = Inf"))
stopifnot(almost.equal(colNorms(M1, "minkowski", p=Inf), apply(tmp, 2, max), name="colNorms, dense, p = Inf"))

tmp <- M1 != 0
stopifnot(almost.equal(rowNorms(M1, "minkowski", p=0), rowSums(tmp), name="rowNorms, dense, Hamming (p = 0)"))
stopifnot(almost.equal(colNorms(M1, "minkowski", p=0), colSums(tmp), name="colNorms, dense, Hamming (p = 0)"))

tmp <- M2^2
stopifnot(almost.equal(rowNorms(M2, "euclidean"), sqrt(rowSums(tmp)), name="rowNorms, sparse, Euclidean"))
stopifnot(almost.equal(colNorms(M2, "euclidean"), sqrt(colSums(tmp)), name="colNorms, sparse, Euclidean"))

tmp <- abs(M2)
stopifnot(almost.equal(rowNorms(M2, "manhattan"), rowSums(tmp), name="rowNorms, sparse, Manhattan"))
stopifnot(almost.equal(colNorms(M2, "manhattan"), colSums(tmp), name="colNorms, sparse, Manhattan"))

stopifnot(almost.equal(rowNorms(M2, "maximum"), apply(tmp, 1, max), name="rowNorms, sparse, maximum"))
stopifnot(almost.equal(colNorms(M2, "maximum"), apply(tmp, 2, max), name="colNorms, sparse, maximum"))
stopifnot(almost.equal(rowNorms(M2, "minkowski", p=Inf), apply(tmp, 1, max), name="rowNorms, sparse, p = Inf"))
stopifnot(almost.equal(colNorms(M2, "minkowski", p=Inf), apply(tmp, 2, max), name="colNorms, sparse, p = Inf"))

tmp <- M2 != 0
stopifnot(almost.equal(rowNorms(M2, "minkowski", p=0), rowSums(tmp), name="rowNorms, sparse, Hamming (p = 0)"))
stopifnot(almost.equal(colNorms(M2, "minkowski", p=0), colSums(tmp), name="colNorms, sparse, Hamming (p = 0)"))

## various Minkowski norms
for (p in c(.1, .2, .5, 1, 1.5, 2, 3, 5, 10)) {
  q <- min(1 / p, 1)

  tmp <- abs(M1) ^ p
  stopifnot(almost.equal(rowNorms(M1, "minkowski", p=p), rowSums(tmp) ^ q, name=paste("rowNorms, dense, Minkowski p =", p)))
  stopifnot(almost.equal(colNorms(M1, "minkowski", p=p), colSums(tmp) ^ q, name=paste("colNorms, dense, Minkowski p =", p)))

  tmp <- abs(M2) ^ p
  stopifnot(almost.equal(rowNorms(M2, "minkowski", p=p), rowSums(tmp) ^ q, name=paste("rowNorms, sparse, Minkowski p =", p)))
  stopifnot(almost.equal(colNorms(M2, "minkowski", p=p), colSums(tmp) ^ q, name=paste("colNorms, sparse, Minkowski p =", p)))
}

## validate row/column normalization (norms must be == 1)
for (norm in c("euclidean", "manhattan", "maximum")) {
  stopifnot(almost.equal(rowNorms(normalize.rows(M1, method=norm), method=norm), 1, name=paste("normalize.rows, dense,", norm)))
  stopifnot(almost.equal(colNorms(normalize.cols(M1, method=norm), method=norm), 1, name=paste("normalize.cols, dense,", norm)))
  stopifnot(almost.equal(rowNorms(normalize.rows(M2, method=norm), method=norm), 1, name=paste("normalize.rows, sparse,", norm)))
  stopifnot(almost.equal(colNorms(normalize.cols(M2, method=norm), method=norm), 1, name=paste("normalize.cols, sparse,", norm)))
}

for (p in c(.1, .2, .5, 1, 1.5, 2, 3, 5, 10)) {
  stopifnot(almost.equal(colNorms(normalize.cols(M1, method="minkowski", p=p), method="minkowski", p=p), 1, name=paste("normalize.cols, dense, Minkowski p =", p)))
  stopifnot(almost.equal(rowNorms(normalize.rows(M1, method="minkowski", p=p), method="minkowski", p=p), 1, name=paste("normalize.rows, dense, Minkowski p =", p)))
  stopifnot(almost.equal(rowNorms(normalize.rows(M2, method="minkowski", p=p), method="minkowski", p=p), 1, name=paste("normalize.rows, sparse, Minkowski p =", p)))
  stopifnot(almost.equal(colNorms(normalize.cols(M2, method="minkowski", p=p), method="minkowski", p=p), 1, name=paste("normalize.cols, sparse, Minkowski p =", p)))
}

## error conditions: normalization not possible / reliable
stopifnot(is(try(normalize.rows(M1, method="minkowski", p=0), silent=TRUE), "try-error"))
stopifnot(is(try(normalize.cols(M1, method="minkowski", p=0), silent=TRUE), "try-error"))
stopifnot(is(try(normalize.rows(M1, method="minkowski", p=.01), silent=TRUE), "try-error"))
stopifnot(is(try(normalize.cols(M1, method="minkowski", p=.01), silent=TRUE), "try-error"))

## near-zero rows/columns are not normalized
fac <- c(1, 0, 1, 1, 1e-9, 0, 1e-12)
M1a <- scaleMargins(M1, rows=fac)
M1b <- scaleMargins(M1, cols=fac)
M2a <- scaleMargins(M2, rows=fac)
M2b <- scaleMargins(M2, cols=fac)

for (p in c(.5, .7, 1, 2, 5)) {
  tol <- if (p < 1) 0.1 else 1e-6 # scaling behaviour is very different for p < 1s
  stopifnot(almost.equal(rowNorms(normalize.rows(M1a, method="minkowski", p=p, tol=tol), method="minkowski", p=p), fac > 1e-6, 
                         name=paste("zeroes in normalize.rows, dense, Minkowski p =", p)))  
  stopifnot(almost.equal(colNorms(normalize.cols(M1b, method="minkowski", p=p, tol=tol), method="minkowski", p=p), fac > 1e-6, 
                         name=paste("zeroes in normalize.cols, dense, Minkowski p =", p)))  
  stopifnot(almost.equal(rowNorms(normalize.rows(M2a, method="minkowski", p=p, tol=tol), method="minkowski", p=p), fac > 1e-6, 
                         name=paste("zeroes in normalize.rows, sparse, Minkowski p =", p)))  
  stopifnot(almost.equal(colNorms(normalize.cols(M2b, method="minkowski", p=p, tol=tol), method="minkowski", p=p), fac > 1e-6, 
                         name=paste("zeroes in normalize.cols, sparse, Minkowski p =", p)))  
}
