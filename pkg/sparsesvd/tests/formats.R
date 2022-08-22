## any sparse matrix format that inherits from dMatrix should work
library(sparsesvd)
library(Matrix)

M <- rbind(
  c(20, 10, 15,  0,  2),
  c(10,  5,  8,  1,  0),
  c( 0,  1,  2,  6,  3),
  c( 1,  0,  0, 10,  5))

M1 <- Matrix(M, sparse=TRUE)          # standard format (column-compressed)
M1 <- as(M1, "CsparseMatrix")
stopifnot(is(M1, "dgCMatrix"))
res1 <- sparsesvd(M1)

M2 <- Matrix(M, sparse=FALSE)         # dense matrix
stopifnot(is(M2, "dgeMatrix"))
res2 <- sparsesvd(M2)

M3 <- Matrix(M, sparse=TRUE)          # triplet format
M3 <- as(M3, "TsparseMatrix")
stopifnot(is(M3, "dgTMatrix"))
res3 <- sparsesvd(M3)

M4 <- Matrix(M, sparse=TRUE)          # row-compressed
M4 <- as(M4, "RsparseMatrix")
stopifnot(is(M4, "dgRMatrix"))
res4 <- sparsesvd(M4)
## TODO -- row-compressed form cannot be converted to dgCMatrix

## check that eigenvalues are the same
stopifnot(all.equal(res1$d, res2$d, tolerance=1e-12))
stopifnot(all.equal(res1$d, res3$d, tolerance=1e-12))
stopifnot(all.equal(res1$d, res4$d, tolerance=1e-12))

## special classes for symmetric matrices
A1 <- crossprod(M1)
A1 <- as(A1, "CsparseMatrix")
stopifnot(is(A1, "dsCMatrix"))
res1a <- sparsesvd(A1)

A2 <- crossprod(M2)
stopifnot(is(A2, "dpoMatrix")) # dense positive-semidefinite symmetric matrix
res2a <- sparsesvd(A2)
## -- fails

A3 <- as(A1, "TsparseMatrix")
stopifnot(is(A3, "dsTMatrix"))
res3a <- sparsesvd(A3)
## -- symmetric matrix can only be converted if already in row-compressed form

A4 <- as(A1, "RsparseMatrix")
stopifnot(is(A4, "dsRMatrix"))
res4a <- sparsesvd(A4)
## -- fails 

## check that eigenvalues are the same and consistent with M
stopifnot(all.equal(res1a$d, (res1$d)^2, tolerance=1e-12))
stopifnot(all.equal(res1a$d, res2a$d, tolerance=1e-12))
stopifnot(all.equal(res1a$d, res3a$d, tolerance=1e-12))
stopifnot(all.equal(res1a$d, res4a$d, tolerance=1e-12))
