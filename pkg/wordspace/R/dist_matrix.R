dist.matrix <- function (M, M2=NULL, method=c("cosine", "euclidean", "maximum", "manhattan", "minkowski", "canberra"), p=2, normalized=FALSE, byrow=TRUE, convert=TRUE, as.dist=FALSE, terms=NULL, terms2=terms, skip.missing=FALSE) {
  method <- match.arg(method)
  similarity <- (method %in% c("cosine")) && !convert
  symmetric <- !(method %in% c()) # FALSE if distance/similarity measure is asymmetric
  cross.distance <- !is.null(M2)  # TRUE if calculating (rectangular) cross-distance matrix

  if (method == "minkowski" && (p < 1 || !is.finite(p))) stop("Minkowski p-norm only defined for 1 <= p < Inf")
  if (as.dist && similarity) stop("cannot create 'dist' object from similarity matrix")
  if (as.dist && cross.distance) stop("cannot create 'dist' object from cross-distance matrix")

  M <- find.canonical.matrix(M) # extract co-occurence matrix from DSM object, ensure canonical format
  sparse.M <- dsm.is.canonical(M)$sparse
  if (cross.distance) {
    M2 <- find.canonical.matrix(M2)
    sparse.M2 <- dsm.is.canonical(M2)$sparse
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

    ## if cross.distance is FALSE, M2 must be NULL and terms2 unspecified (i.e. M2=M and terms2=terms), so leave M2 set to NULL
    if (cross.distance) {
      if (!is.null(terms2)) {
        terms2 <- as.character(terms2) # in case terms2 is a factor
        found <- terms2 %in% targets.M2
        if (!all(found) && !skip.missing) stop("second term(s) not found in M2: ", paste(terms2[!found], collapse=", "))
        terms2 <- terms2[found]
        if (is.null(M2)) {
          M2 <- if (byrow) M[terms2, , drop=FALSE] else M[ , terms2, drop=FALSE] # need to process terms2 first before overwriting M below
        } else {
          M2 <- if (byrow) M2[terms2, , drop=FALSE] else M2[ , terms2, drop=FALSE]
        }
      } else {
        if (is.null(M2)) M2 <- M # cross-distances with terms2=NULL -> M2 = copy of M before subsetting
      }
      sparse.M2 <- dsm.is.canonical(M2)$sparse # define/update sparse.M2
    }

    if (!is.null(terms)) {
      terms <- as.character(terms) # in case terms is a factor
      found <- terms %in% targets.M
      if (!all(found) && !skip.missing) stop("first term(s) not found in M: ", paste(terms[!found], collapse=", "))
      terms <- terms[found]
      M <- if (byrow) M[terms, , drop=FALSE] else M[ , terms, drop=FALSE]
    }
  
  }
  
  if (method == "cosine") {
    ## cosine / angular measure is computed as very efficient matrix crossproduct
    
    if (byrow) {
      result <- if (is.null(M2)) tcrossprod(M) else tcrossprod(M, M2)
    } else {
      result <- if (is.null(M2)) crossprod(M) else crossprod(M, M2)
    }
    result <- as.matrix(result) # ensure that cosine similarity matrix is in dense representation
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
      n <- prod(dim(result))
      transform_code <- 0L # cosine -> angle transformation
      n.clamped <- 0L
      tol <- 1e-12
      .C(
        C_similarity_to_distance,
        result,
        as.integer(n),
        transform_code,
        tol,
        n.clamped,
        DUP=FALSE, NAOK=FALSE
      )
      if (n.clamped > 0) warning("angular distance may be inaccurate (some cosine values out of range)")
      
      # tol <- 1e-12 # rounding errors tolerated for very small angles (cosine approx. 1 or -1)
      # if(!all(result >= -(1+tol) & result <= 1+tol)) warning("angular distance may be inaccurate (some cosine values out of range)")
      # ## TODO: rewrite angle computation as inplace operation in C to avoid memory overhead
      # result[result < -(1-tol)] <- -1     # clamp to range [-1, 1] and snap to endpoints -1 / 1
      # result[result > 1-tol] <- 1         # (pmin/pmax eat many GiB of memory for Matrix class??)
      # result <- acos(result) * (180 / pi) # angles are returned in degrees
    }    
    rownames(result) <- if (byrow) rownames(M) else colnames(M)
    colnames(result) <- if (is.null(M2)) rownames(result) else if (byrow) rownames(M2) else colnames(M2)

  } else {
    ## other distance measures are implemented in C code, working on columns (transposed matrix) for efficiency

    .M <- if (byrow) t(M) else M
    if (cross.distance) {
      if (sparse.M != sparse.M2) stop("M and M2 must either be both in dense format or both in sparse format")
      .M2 <- if (byrow) t(M2) else M2
    } else {
      .M2 <- .M
    }

    method.code <- switch(method, euclidean=0, maximum=1, manhattan=2, minkowski=3, canberra=4) # must be kept in sync with C code
    param1 <- switch(method, euclidean=0, maximum=0, manhattan=0, minkowski=p, canberra=0)
  
    if (sparse.M) {
      result <- matrix(0.0, nrow=ncol(.M), ncol=ncol(.M2))
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
      result <- CPP_col_dist_dense(.M, .M2, method.code, param1, symmetric && !cross.distance)
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
