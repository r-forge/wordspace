dsm.projection <- function (model, method=c("svd", "rsvd", "asvd", "ri", "ri+svd"), n=NA, oversampling=NA, q=2, rate=.01, verbose=FALSE, with.basis=FALSE, blocksize=100000, use.C=FALSE) {
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

  is.sparse <- is(M, "Matrix") # whether to use sparse or dense algorithms
  
  nR <- nrow(M)
  nC <- ncol(M)
  
  if (is.na(n)) n <- min(nR, nC)
  if (method == "ri") {
    if (2*n > nC) stop(sprintf("random indexing from %d to %d dimensions makes no sense", nC, n))
  } else {
    if (n > min(nR, nC)) stop("number of target dimensions exceeds dimensionality of DSM matrix")
  }

  if (is.na(oversampling)) {
    oversampling <- switch(method, svd=1, rsvd=2, asvd=10, ri=1, "ri+svd"=20)
  }
  if (with.basis && method %in% c("ri", "ri+svd")) stop("with.basis=TRUE is not available for RI-based models")
  
  if (method == "svd") {
    ## --- standard SVD algorithm (on dense matrix) ---

    if (verbose) cat(sprintf("SVD reduction to %d dimensions:\n", n))
    if (verbose) cat(" - SVD decomposition\n")
    SVD <- svd(M, nu=n, nv=(if (with.basis) n else 0)) # we don't need right singular vectors for the dimensionality reduction
    if (verbose) cat(" - composing final matrix\n")
    S <- SVD$u %*% diag(SVD$d[1:n], n, n) # dimensionality-reduced matrix
    if (with.basis) B <- SVD$v
    R2 <- SVD$d[1:n] * SVD$d[1:n] / sum(SVD$d * SVD$d) # proportion of variance "explained" by SVD dimensions
    rm(SVD)

    ## *** TODO: use SVDLIBC for efficient sparse SVD if available ***
    
  } else if (method == "asvd") {
    ## --- approximated SVD based on random sample of rows ---

    sub.size <- min(n * oversampling, nR)
    if (verbose) cat(sprintf("Approximate SVD reduction to %d dimensions, based on %d rows:\n", n, sub.size))
    sub.idx <- sort(sample(1:nR, sub.size, replace=FALSE))
    M.sub <- M[sub.idx, ]
    if (verbose) cat(" - SVD decomposition\n")
    SVD <- svd(M.sub, nu=0, nv=n) # here we only need the right singular vectors
    if (verbose) cat(" - composing final matrix\n")
    S <- M %*% SVD$v  # V projects columns to first n latent dimensions
    if (with.basis) B <- SVD$v
    rm(SVD)
    S <- as.matrix(S) # make sure result is an ordinary dense matrix (for efficient further processing)
    R2 <- cumsum(colSums(S*S)) / norm(M, "F")^2 # this should be the proportion of "explained" variance

  } else if (method == "rsvd") {
    ## --- randomized SVD according to Halko, Martinsson & Tropp (2009, p. 9) ---

    ## We can apply the rSVD algorithm either to A = M or to A = t(M), depending on the format of M.
    ## Preliminary testing suggested that the original algorithm (A = M) is suitable for a matrix with many columns,
    ## while the transpose algorithm (A = t(M)) works better if the matrix has many rows and a limited number of columns.
    ## With the current implementation, which uses SVD rather than QR decomposition to obtain an orthonormal basis,
    ## there does not seem to be a substantial difference, so we currently always use the original algorithm.
    k2 <- min(oversampling * n, nR, nC)   # = 2*k in the paper
    transpose <- FALSE # previously: !(nR <= nC)
    if (verbose) cat(sprintf("Randomized SVD reduction%s to %d => %d dimensions:\n", if (transpose) " (transposed)" else "", k2, n))
    if (!transpose) {
      if (verbose) cat(" - sampling range of A\n") # -- original algorithm applied to A = M
      Omega <- matrix(rnorm(nC*k2), nC, k2)
      Y <- M %*% Omega
      rm(Omega)
      if (q >= 1) for (i in 1:q) {
        if (verbose) cat(sprintf(" - power iteration #%d\n", i))
        Y <- M %*% crossprod(M, Y)
      }
      if (verbose) cat(sprintf(" - orthonormal basis of %d x %d matrix\n", nrow(Y), ncol(Y)))
      Q <- svd(Y, nu=k2, nv=0)$u  # orthonormal basis of rg(Y); SVD is faster than and as accurate as QR decomposition
      rm(Y)
      B <- crossprod(Q, M)
      if (verbose) cat(sprintf(" - SVD decomposition of %d x %d matrix\n", nrow(B), ncol(B)))
      SVD <- svd(B, nu=0, nv=n)   # we can either project the original matrix using V, or construct the result from U and Sigma
      rm(B)
      if (verbose) cat(" - composing final matrix\n")
      ## We can either construct the reduced matrix from SVD$u and SVD$d, or project original matrix using SVD$v
      ## U <- Q %*% SVD$u # need to set nu=n above if we want to run this code (Uhat = SVD$u)
      ## --> now construct projected matrix from U and SVD$d (check paper)
      S <- M %*% SVD$v
      if (with.basis) B <- SVD$v
      R2 <- cumsum(SVD$d^2) / norm(M, "F")^2
      rm(SVD)
      S <- as.matrix(S) # make sure result is an ordinary dense matrix
    } else {
      if (verbose) cat(" - sampling range of A\n") # -- transposed algorithm for A = t(M)
      Omega <- matrix(rnorm(k2*nR), k2, nR) # = t(Omega)
      Y <- Omega %*% M                      # = t(A * Omega)
      rm(Omega)
      if (q >= 1) for (i in 1:q) {
        if (verbose) cat(sprintf(" - power iteration #%d\n", i))
        Y <- tcrossprod(Y, M) %*% M         # = t( (A * t(A))^i * A * Omega) ) = t(Y)
      }
      if (verbose) cat(sprintf(" - orthonormal basis of %d x %d matrix\n", ncol(Y), nrow(Y)))
      Q <- svd(Y, nu=0, nv=k2)$v  # orthonormal basis of rg(t(Y)); SVD is faster than and as accurate as QR decomposition
      rm(Y)
      B <- M %*% Q                          # = t( t(Q) * A ) = t(B)
      if (verbose) cat(sprintf(" - SVD decomposition of %d x %d matrix\n", nrow(B), ncol(B)))
      SVD <- svd(B, nu=0, nv=n)             # t(B) = V * Sigma * t(Uhat), truncated to k target dimensions
      rm(B)
      if (verbose) cat(" - composing final matrix\n")
      S <- M %*% Q %*% SVD$v                # t(t(U) * A) = t(A) * U = t(A) * Q * Uhat
      if (with.basis) B <- Q %*% SVD$v
      R2 <- cumsum(SVD$d^2) / norm(M, "F")^2
      rm(SVD)
      S <- as.matrix(S) # make sure result is an ordinary dense matrix
    }

  } else if (method %in% c("ri", "ri+svd")) {
    ## --- straightforward random indexing (with specified fill rate), optionally followed by rSVD ---

    ## *** TODO: check references on statistical guarantees, appropriate fill rates, etc. ***
    if (method == "ri+svd") {
      nRI <- n * oversampling       # number of intermediate random dimensions
      nRI <- max(2*n, min(nRI, floor(nC / 2)))
    } else {
      nRI <- n                      # number of random dimensions
    }
    if (verbose) cat(sprintf("Random Indexing in %d dimensions:\n", nRI))
    
    if (!is.sparse) {
      ## if original matrix can be stored in dense representation, RI should not pose any memory problems
      n.fill <- as.integer(nC * rate) # number of nonzero entries in each random vector
      scale <- 1 / sqrt(n.fill)       # scale random vectors so they are normalised

      if (verbose) cat(sprintf(" - generating %d random vectors with %d of %d nonzero elements\n", nRI, n.fill, nC))
      Q <- matrix(0, nRI, nC)
      for (.d in 1:nRI) {
        .idx <- sort(sample.int(nC, n.fill))
        Q[.d, .idx] <- scale * (1 - 2*rbinom(n.fill, 1, .5))
      }

      if (verbose) cat(" - projecting into random subspace\n")
      S <- tcrossprod(M, Q)
      rm(Q)
    } else {
      
      if (use.C) {
        ## experimental C implementation of RI for sparse matrix
        if (!is(M, "dgCMatrix")) stop("sparse matrix must be in normal form (dgCMatrix) for random indexing")
        S <- matrix(0.0, nrow=nR, ncol=nRI) # pre-allocate projected matrix for C code
        .C(
          C_random_indexing_sparse,
          S,  # will be modified inplace
          as.integer(nR),
          as.integer(nC),
          as.integer(M@p),
          as.integer(M@i),
          as.double(M@x),
          as.integer(nRI),
          as.double(rate),
          DUP=FALSE, NAOK=FALSE
        )
      } else {
        ## standard R implementation of RI for sparse matrix
        ## for sparse matrix, perform RI in blocks of <blocksize> dimensions, iteratively updating the projected vectors    
        S <- matrix(0.0, nrow=nR, ncol=nRI) # pre-allocate projected matrix for block updates
        cumulative.fill <- 0 # total number of nonzero entries (+1 / -1) in random projection vectors
        for (start in seq(1, nC, blocksize)) {
          end <- min(start + blocksize - 1, nC)
          n.cols <- end - start + 1
          n.fill <- as.integer(n.cols * rate) # number of nonzero entries in each random vector
          cumulative.fill <- cumulative.fill + n.fill

          if (verbose) cat(sprintf(" - generating %d sparse random vectors with %d of %d nonzero elements\n", nRI, n.fill, n.cols))
          .i <- integer(n.fill * nRI) # 0-based row offsets of nonzero entries
          for (.d in 0:(nRI-1)) {
            .i[ .d * n.fill + (1:n.fill) ] <- sort(sample.int(n.cols, n.fill)) - 1L
          }
          .p <- n.fill * (0:nRI)    # 0-based offset pointers into .i
          .x <- 1 - 2*rbinom(n.fill * nRI, 1, .5)
          Qt <- new("dgCMatrix", Dim=as.integer(c(n.cols, nRI)), p=.p, i=.i, x=.x)

          if (verbose) cat(" - updating projection into random subspace\n")
          S <- S + as.matrix(M[,start:end] %*% Qt)  # make sure result is a dense matrix
          rm(Qt, .i, .p, .x)
        }

        # all random vectors have Euclidean length sqrt(cumulative.fill), so adjust RI coordinates
        S <- S / sqrt(cumulative.fill)  
      }

    }
    
    if (method == "ri+svd") {
      S <- dsm.projection(S, "rsvd", n, q=q, oversampling=2, verbose=verbose, with.basis=FALSE)
    }
    R2 <- norm(S, "F")^2 / norm(M, "F")^2 
   
  } else {
    stop("dimensionality reduction method '", method, "' has not been implemented yet")
  }

  rownames(S) <- rownames(M)
  colnames(S) <- paste(method, 1:n, sep="")
  if (with.basis) {
    rownames(B) <- colnames(M)
    colnames(B) <- colnames(S)
    attr(S, "basis") <- B
  }
  attr(S, "R2") <- R2
  return(S)
 }
