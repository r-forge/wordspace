dsm.projection <- function (model, method=c("svd", "rsvd", "asvd", "ri", "ri+svd"), n=NA, oversampling=NA, q=2, verbose=FALSE) {
  method <- match.arg(method)
  if (is.matrix(model) || inherits(model, "Matrix")) {
    M <- model
  } else if (inherits(model, "dsm")) {
    model.info <- check.dsm(model)
    if (!model.info$have.S) stop("use dsm.score() to compute scored and weighted matrix first")
    M <- model$S
  } else {
    stop("first argument of dsm.projection() must be DSM object or numeric matrix")
  }

  nR <- nrow(M)
  nC <- ncol(M)
  
  if (is.na(n)) n <- min(nR, nC)
  if (n > min(nR, nC)) stop("number of target dimensions exceeds dimensionality of DSM matrix")

  if (is.na(oversampling)) {
    oversampling <- switch(method, svd=1, rsvd=2, asvd=5, ri=1, "ri+svd"=5)
  }

  if (method == "svd") {
    ## -- standard SVD algorithm (dense matrix)

    if (verbose) cat(sprintf("SVD reduction to %d dimensions:\n", n))
    if (verbose) cat(" - SVD decomposition\n")
    SVD <- svd(M, nu=n, nv=0) # we don't need right singular vectors for the dimensionality reduction
    if (verbose) cat(" - composing final matrix\n")
    S <- SVD$u %*% diag(SVD$d[1:n], n, n) # dimensionality-reduced matrix
    R2 <- SVD$d[1:n] * SVD$d[1:n] / sum(SVD$d * SVD$d) # proportion of variance "explained" by SVD dimensions

  } else if (method == "asvd") {
    ## -- approximated SVD based on random sample of rows

    sub.size <- min(n * oversampling, nR)
    if (verbose) cat(sprintf("Approximate SVD reduction to %d dimensions, based on %d rows:\n", n, sub.size))
    sub.idx <- sort(sample(1:nR, sub.size, replace=FALSE))
    M.sub <- M[sub.idx, ]
    if (verbose) cat(" - SVD decomposition\n")
    SVD <- svd(M.sub, nu=0, nv=n) # here we only need the right singular vectors
    if (verbose) cat(" - composing final matrix\n")
    S <- M %*% SVD$v  # V projects columns to first n latent dimensions
    rm(SVD)
    S <- as.matrix(S) # make sure result is an ordinary dense matrix (for efficient further processing)
    R2 <- norm(S, "F")^2 / norm(M, "F")^2 # this should be the proportion of "explained" variance

  } else if (method == "rsvd") {

    ## -- randomized SVD according to Halko, Martinsson & Tropp (2009, p. 9)
    ## we apply this algorithm to A = t(M), which is suitable for M with many rows and limited number of columns
    ## TODO: should automatically apply either to A = M or to A = t(M) depending on format of the matrix
    k2 <- min(oversampling * n, nR, nC)   # = 2*k in the paper
    if (verbose) cat(sprintf("Randomized SVD reduction to %d => %d dimensions:\n", k2, n))
    if (verbose) cat(" - sampling range of A\n")
    Omega <- matrix(rnorm(k2*nR), k2, nR) # = t(Omega)
    Y <- Omega %*% M                      # = t(A * Omega)
    rm(Omega)
    if (q >= 1) for (i in 1:q) {
      if (verbose) cat(sprintf(" - power iteration #%d\n", i))
      Y <- tcrossprod(Y, M) %*% M         # = t( (A * t(A))^i * A * Omega) ) = t(Y)
    }
    if (verbose) cat(" - QR decomposition\n")
    Q <- qr.Q( qr(t(Y)) )                 # = orthonormal basis of rg(Y), from QR decomposition
    rm(Y)
    B <- M %*% Q                          # = t( t(Q) * A ) = t(B)
    if (verbose) cat(" - SVD decomposition\n")
    SVD <- svd(B, nu=n, nv=n)             # t(B) = V * Sigma * t(Uhat), truncated to k target dimensions
    rm(B)
    if (verbose) cat(" - composing final matrix\n")
    S <- M %*% Q %*% SVD$v                # t(t(U) * A) = t(A) * U = t(A) * Q * Uhat
    rm(SVD)
    S <- as.matrix(S) # make sure result is an ordinary dense matrix
    R2 <- norm(S, "F")^2 / norm(M, "F")^2

  } else {
    stop("dimensionality reduction method '", method, "' has not been implemented yet")
  }

  rownames(S) <- rownames(M)
  colnames(S) <- paste(method, 1:n, sep="")
  attr(S, "R2") <- R2
  return(S)
 }
