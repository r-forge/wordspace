as.distmat <- function (x, ...)  UseMethod("as.distmat")

## default method for matrix 
as.distmat.matrix <- function (x, similarity=FALSE, symmetric=FALSE, ...) {
  if (!("dist.matrix" %in% class(x))) class(x) <- c("dist.matrix", class(x))
  attr(x, "dist.matrix") <- TRUE
  attr(x, "similarity") <- similarity
  if (symmetric) {
    if (nrow(x) != ncol(x)) stop("non-square matrix cannot be symmetric")
    if (!is.null(rownames(x)) && !is.null(colnames(x)) && !all(rownames(x) == colnames(x))) stop("matrix cannot be symmetric because row and column labels differ")
    attr(x, "symmetric") <- symmetric
  }
  x
}

as.distmat.sparseMatrix <- function (x, similarity=FALSE, symmetric=FALSE, force.dense=FALSE, ...) {
  if (force.dense) {
    x <- as.matrix(x)
    class(x) <- c("dist.matrix", class(x))
  }
  attr(x, "dist.matrix") <- TRUE
  attr(x, "similarity") <- similarity
  if (symmetric) {
    if (nrow(x) != ncol(x)) stop("non-square matrix cannot be symmetric")
    if (!is.null(rownames(x)) && !is.null(colnames(x)) && !all(rownames(x) == colnames(x))) stop("matrix cannot be symmetric because row and column labels differ")
    attr(x, "symmetric") <- symmetric
  }
  x
}
