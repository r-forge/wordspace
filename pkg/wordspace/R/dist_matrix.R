dist.matrix <- function (M, M2=NULL, method=c("cosine", "euclidean", "maximum", "manhattan", "minkowski", "canberra"), p=2, normalized=FALSE, byrow=TRUE, convert=TRUE, as.dist=FALSE, terms=NULL, terms2=terms, skip.missing=FALSE) {
  method <- match.arg(method)
  similarity <- (method %in% c("cosine")) && !convert
  symmetric <- !(method %in% c()) # FALSE if distance/similarity measure is asymmetric
  cross.distance <- !is.null(M2)  # TRUE if calculating (rectangular) cross-distance matrix

  if (method == "minkowski" && (p < 1 || !is.finite(p))) stop("Minkowski p-norm only defined for 1 <= p < Inf")
  if (as.dist && similarity) stop("cannot create 'dist' object from similarity matrix")
  if (as.dist && cross.distance) stop("cannot create 'dist' object from cross-distance matrix")

  sparse.M <- inherits(M, "Matrix") # check that M (and optional M2) are appropriate matrices
  if (!(is.matrix(M) || sparse.M)) stop("M must be a dense or sparse matrix")
  if (cross.distance) {
    sparse.M2 <- inherits(M2, "Matrix")
    if (!(is.matrix(M2) || sparse.M2)) stop("M2 must be a dense or sparse matrix")
    if (byrow) {
      if (ncol(M) != ncol(M2)) stop("M and M2 are not conformable (must have same number of columns)")
    } else {
      if (nrow(M) != nrow(M2)) stop("M and M2 are not conformable (must have same number of rows)")        
    }
  }

  if (!is.null(terms) || !is.null(terms2)) {
    targets.M <- if (byrow) rownames(M) else colnames(M)
    targets.M2 <- if (is.null(M2)) targets.M else if (byrow) rownames(M2) else colnames(M2)

    if (!missing(terms2)) cross.distance <- TRUE # if different filters are applied, we're always dealing with a cross-distance calculation

    # if cross.distance is FALSE, both M2 and terms2 must be missing (and hence M2=M and terms2=terms), so leave M2 set to NULL
    if (!is.null(terms2) && cross.distance) { 
      terms2 <- as.character(terms2) # in case terms2 is a factor
      found <- terms2 %in% targets.M2
      if (!all(found) && !skip.missing) stop("second term(s) not found in M2: ", paste(terms2[!found], collapse=", "))
      terms2 <- terms2[found]
      if (is.null(M2)) {
        M2 <- if (byrow) M[terms2, , drop=FALSE] else M[ , terms2, drop=FALSE] # need to process terms2 first before overwriting M below
      } else {
        M2 <- if (byrow) M2[terms2, , drop=FALSE] else M2[ , terms2, drop=FALSE]
      }
    }

    if (!is.null(terms)) {
      terms <- as.character(terms) # in case terms is a factor
      found <- terms %in% targets.M
      if (!all(found) && !skip.missing) stop("first term(s) not found in M: ", paste(terms[!found], collapse=", "))
      terms <- terms[found]
      M <- if (byrow) M[terms, , drop=FALSE] else M[ , term2, drop=FALSE]
    }
  
  }

  if (method == "cosine") {
    ## cosine / angular measure is computed as very efficient matrix crossproduct
    
    if (byrow) {
      result <- if (is.null(M2)) tcrossprod(M) else tcrossprod(M, M2)
    } else {
      result <- if (is.null(M2)) crossprod(M) else crossprod(M, M2)
    }
    result <- as.matrix(result) # ensure that cosine similarity matrix is dense matrix
    if (!normalized) {
      norms.M <- if (byrow) rowNorms(M, "euclidean") else colNorms(M, "euclidean")
      if (is.null(M2)) {
        norms.M2 <- norms.M
      } else {
        norms.M2 <- if (byrow) rowNorms(M2, "euclidean") else colNorms(M2, "euclidean")        
      }
      result <- scaleMargins(result, rows=1/norms.M, cols=1/norms.M2)
    }

    if (convert) {
      tol <- 1e-12 # rounding errors tolerated for very small angles (cosine approx. 1 or -1)
      if(!all(result >= -(1+tol) & result <= 1+tol)) warning("angular distance may be inaccurate (some cosine values out of range)")
      ## TODO: rewrite angle computation as inplace operation in C to avoid memory overhead
      result[result < -(1-tol)] <- -1     # clamp to range [-1, 1] and snap to endpoints -1 / 1
      result[result > 1-tol] <- 1         # (pmin/pmax eat many GiB of memory for Matrix class??)
      result <- acos(result) * (180 / pi) # angles are returned in degrees
    }    
    rownames(result) <- if (byrow) rownames(M) else colnames(M)
    colnames(result) <- if (is.null(M2)) rownames(result) else if (byrow) rownames(M2) else colnames(M2)

  } else {
    ## other distance measures are implemented in C code, working on columns (transposed matrix) for efficiency

    if (sparse.M && !is(M, "dgCMatrix")) stop("sparse matrix M must be in normal form (dgCMatrix)")
    .M <- if (byrow) t(M) else M

    if (cross.distance) {
      if (sparse.M2 && !is(M2, "dgCMatrix")) stop("sparse matrix M2 must be in normal form (dgCMatrix)")
      if (sparse.M != sparse.M2) stop("M and M2 must be both in dense format or both in sparse format")
      .M2 <- if (byrow) t(M2) else M2
    } else {
      .M2 <- .M
    }

    method.code <- switch(method, euclidean=0, maximum=1, manhattan=2, minkowski=3, canberra=4) # must be kept in sync with C code
    param1 <- switch(method, euclidean=0, maximum=0, manhattan=0, minkowski=p, canberra=0)
  
    result <- matrix(0.0, nrow=ncol(.M), ncol=ncol(.M2))
    if (sparse.M) {
      .C(
        C_col_dist_sparse,
        result,
        as.integer(ncol(.M)),
        as.integer(ncol(.M2)),
        as.integer(.M@p),
        as.integer(.M@i),
        as.double(.M@x),
        as.integer(.M2@p),
        as.integer(.M2@i),
        as.double(.M2@x),
        as.integer(method.code),
        as.double(param1),
        as.logical(symmetric && !cross.distance),
        DUP=FALSE, NAOK=FALSE
      )    
    } else {
      .C(
        C_col_dist_dense,
        result,
        as.integer(nrow(.M)),
        as.integer(ncol(.M)),
        as.integer(ncol(.M2)),
        as.double(.M),
        as.double(.M2),
        as.integer(method.code),
        as.double(param1),
        as.logical(symmetric && !cross.distance),
        DUP=FALSE, NAOK=FALSE
      )
    }
    rownames(result) <- colnames(.M)
    colnames(result) <- colnames(.M2)

  }

  if (as.dist) {
    as.dist(result)
  } else {
    class(result) <- c("dist.matrix", "matrix")
    if (similarity) attr(result, "similarity") <- TRUE
    result
  }
}

cosine <- function (M, M2=M, angles=FALSE, normalized=FALSE) {
  # tcrossprod(M, M2) == M %*% t(M2) calculates dot products between rows of M and rows of M2
  sim <- if (missing(M2)) tcrossprod(M) else tcrossprod(M, M2)
  # need to coerce to regular matrix if sparse (for pmin/pmax, but generally more efficient)
  sim <- as.matrix(sim) # ensure this is a dense matrix
  if (!normalized) {
    norms.M <- rowNorms(M, "euclidean") # norms of row vectors (if not normalised yet)
    norms.M2 <- if (missing(M2)) norms.M else rowNorms(M2, "euclidean")
    sim <- sim / outer(norms.M, norms.M2)
  }

  if (angles) {
    stopifnot(all(sim >= -(1+1e-6) & sim <= 1+1e-6)) # check that cosines are in appropriate range
    sim[sim < -1] <- -1           # clamp to range [-1, 1]
    sim[sim > 1] <- 1             # (pmin/pmax eat many GiB of memory for Matrix class??)
    sim <- acos(sim) / pi * 180   # angles are returned in degrees
  }

  rownames(sim) <- rownames(M)
  colnames(sim) <- rownames(M2)
  return(sim)
}

