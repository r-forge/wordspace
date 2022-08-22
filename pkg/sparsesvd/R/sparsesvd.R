sparsesvd <- function (M, rank=0L, tol=1e-15, kappa=1e-6) {
  if (is.matrix(M)) stop("argument must be a sparse real matrix")
  if (!is(M, "dMatrix")) stop("only sparse real dMatrix format (from Matrix package) is currently supported")
  if (!is(M, "dgCMatrix")) {
    ## coerce to dgCMatrix format via virtual base classes because
    ##  - direct conversion has been deprecated since Matrix 1.4-2
    ##  - direct path is not available for all pairs of classes anyway
    if (!is(M, "CsparseMatrix")) M <- as(M, "CsparseMatrix")
    if (!is(M, "generalMatrix")) M <- as(M, "generalMatrix")
    if (!is(M, "dgCMatrix")) stop("unable to convert M to dgCMatrix format")
  }
  .Call(svdLAS2_, dim(M), M@i, M@p, M@x, as.integer(rank), as.double(tol * c(-1, 1)), as.double(kappa), PACKAGE="sparsesvd")
}
