## Test validation of and conversion to canonical matrix formats
library(wordspace)

M.dsC <- tcrossprod(DSM_TermContextMatrix) # create different sparse and dense formats for a symmetric matrix
expect_s4_class(M.dsC, "dsCMatrix")

M.dsT <- as(M.dsC, "TsparseMatrix")
expect_s4_class(M.dsT, "dsTMatrix")

M.dgC <- as(M.dsC, "generalMatrix")
expect_s4_class(M.dgC, "dgCMatrix")

M.dgR <- as(M.dgC, "RsparseMatrix")
expect_s4_class(M.dgR, "dgRMatrix")

M.dgT <- as(M.dgC, "TsparseMatrix")
expect_s4_class(M.dgT, "dgTMatrix")

M.dsy <- as(M.dsC, "denseMatrix")
expect_s4_class(M.dsy, "dsyMatrix")

M.dge <- as(M.dsy, "generalMatrix")
expect_s4_class(M.dge, "dgeMatrix")

M.mat <- as.matrix(M.dsC)
expect_true(is.matrix(M.mat))

## canonical format checks
check.canonical <- function (x, canonical=FALSE, sparse=TRUE, nonneg=TRUE, check.nonneg=TRUE) {
  label <- deparse(substitute(x))
  res <- dsm.is.canonical(x, nonneg.check=check.nonneg)
  expect_equal(res$canonical, canonical, label=paste0("dsm.is.canonical(", label, ")$canonical"))
  expect_equal(res$sparse, sparse, label=paste0("dsm.is.canonical(", label, ")$sparse"))
  if (check.nonneg) expect_equal(res$nonneg, nonneg, label=paste0("dsm.is.canonical(", label, ")$nonneg"))
}

check.canonical(M.dsC)
check.canonical(M.dsT)
check.canonical(M.dgC, canonical=TRUE)
check.canonical(M.dgR)
check.canonical(M.dgT)
check.canonical(M.dsy, sparse=FALSE, check.nonneg=FALSE)
check.canonical(M.dge, sparse=FALSE)
check.canonical(M.mat, sparse=FALSE, canonical=TRUE)

mkneg <- function (x) { x[4, 2] <- -42; x }
check.canonical(mkneg(M.dsT), nonneg=FALSE)
check.canonical(mkneg(M.dgC), canonical=TRUE, nonneg=FALSE)
check.canonical(mkneg(M.dgR), nonneg=FALSE)
check.canonical(mkneg(M.dgT), nonneg=FALSE)
check.canonical(mkneg(M.dge), sparse=FALSE, nonneg=FALSE)
check.canonical(mkneg(M.mat), sparse=FALSE, canonical=TRUE, nonneg=FALSE)

## conversion to canonical format
check.conversion <- function (x, sparse=TRUE, triplet=FALSE, reference=M.mat) {
  label <- deparse(substitute(x))
  y <- dsm.canonical.matrix(x, triplet=triplet)
  if (sparse) {
    expected.class <- if (triplet) "dgTMatrix" else "dgCMatrix"
    expect_s4_class(y, expected.class)
  }
  else {
    expect_true(is.matrix(y), label=paste0("is.matrix(", label, ")"))
  }
  y <- as.matrix(y)
  expect_equal(y, reference, label=label)
}

check.conversion(M.dsC)
check.conversion(M.dsT)
check.conversion(M.dgC)
check.conversion(M.dgR)
check.conversion(M.dgT)
check.conversion(M.dsy, sparse=FALSE)
check.conversion(M.dge, sparse=FALSE)
check.conversion(M.mat, sparse=FALSE)

check.conversion(M.dsC, triplet=TRUE)
check.conversion(M.dsT, triplet=TRUE)
check.conversion(M.dgC, triplet=TRUE)
check.conversion(M.dgR, triplet=TRUE)
check.conversion(M.dgT, triplet=TRUE)
