dsm.projection <- function (model, method=c("svd", "rsvd", "asvd", "ri", "ri+svd"), n=NA, oversampling=NA, q=2, rate=.01, with.basis=FALSE, verbose=FALSE) {
  method <- match.arg(method)
  M <- find.canonical.matrix(model)
  info <- dsm.is.canonical(M) # to find out whether M is sparse or dense
  
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

    ## randomized SVD according to Halko, Martinsson & Tropp (2009, p. 9) 
    ##  - preliminary testing suggests there is no substantial difference between the original and transposed rSVD algorithm
    ##  - so we currently always use the original versions
    SVD <- rsvd(M, n=n, q=q, oversampling=oversampling, transpose=FALSE, verbose=verbose)

    S <- scaleMargins(SVD$u, cols=SVD$d)
    ## -- still necessary? --
    ## S <- as.matrix(S) # make sure result is an ordinary dense matrix
    if (with.basis) B <- SVD$v
    R2 <- cumsum(SVD$d^2) / norm(M, "F")^2
    rm(SVD)

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
    
    if (!info$sparse) {
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
      ## efficient C implementation of RI for sparse matrix (which is guaranteed to be in canonical format)
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
        as.logical(verbose),
        DUP=FALSE, NAOK=FALSE
        )
    }
    
    if (method == "ri+svd") {
      S <- dsm.projection(S, "svd", n, verbose=verbose, with.basis=FALSE) # use plain SVD since matrix is already dense
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
