nearest.neighbours <-function (M, term, n=10, drop=TRUE, skip.missing=FALSE, min.distance=1e-6, byrow=TRUE, ...) {
  is.dist <- inherits(M, "dist.matrix")
  if (!is.dist) {
    if (!(is.matrix(M) || is(M, "Matrix"))) stop("M must be either a dense or sparse DSM matrix, or a pre-computed dist.matrix object")
  }

  term.labels <- if (byrow) rownames(M) else colnames(M)
  found <- term %in% term.labels
  if (any(!found) && !skip.missing) stop("target term(s) not found in M: ", paste(term[!found], collapse=", "))
  term <- term[found]
  if (length(term) == 0) return(NULL)
  
  if (!is.dist) {
    M.term <- if (byrow) M[term,, drop=FALSE] else M[,term, drop=FALSE]
    M <- dist.matrix(M=M.term, M2=M, byrow=byrow, ...)
    byrow <- TRUE  # dist.matrix computed on the fly is always accessed by row
  }
  similarity <- attr(M, "similarity")
  if (is.null(similarity)) similarity <- FALSE
  
  result <- lapply(term, function (.t) {
    neighbours <- if (byrow) M[.t,] else M[,.t]
    if (!similarity && min.distance > 0) neighbours <- neighbours[neighbours >= min.distance]
    neighbours <- sort(neighbours, decreasing=similarity)
    head(neighbours, n)
  })
  names(result) <- term
  
  if (drop && length(result) == 1) result[[1]] else result
}
