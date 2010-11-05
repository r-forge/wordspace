dsm.projection <- function (model, method=c("svd"), n=NA) {
  method <- match.arg(method)
  model.info <- check.dsm(model)
  if (!model.info$have.S) stop("use dsm.score() to compute scored and weighted matrix first")
  if (missing(n)) n <- model.info$nrow
  if (n > min(model.info$nrow, model.info$ncol)) stop("number of target dimensions exceeds dimensionality of DSM matrix")
  
  SVD <- svd(model$S, nu=n, nv=0) # we don't need right singular vectors for the dimensionality reduction
  S <- SVD$u %*% diag(SVD$d[1:n], n, n) # dimensionality-reduced matrix
  rownames(S) <- model$rows$term
  colnames(S) <- paste("svd", 1:n, sep="")
  R2 <- SVD$d[1:n] * SVD$d[1:n] / sum(SVD$d * SVD$d) # proportion of variance "explained" by SVD dimensions

  attr(S, "R2") <- R2
  return(S)
}
