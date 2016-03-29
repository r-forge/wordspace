sparsesvd <- function (M) {
  if (is.matrix(M)) stop("dense matrix not allowed")
  if (!is(M, "dMatrix")) stop("wrong format")
  if (!is(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  .Call(svdLAS2_, dim(M), M@i, M@p, M@x)
}
