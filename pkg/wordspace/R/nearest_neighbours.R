nearest.neighbours <-function (M, term, n=10, drop=TRUE, skip.missing=FALSE, min.distance=1e-6, byrow=TRUE, ..., batchsize=50e6, verbose=FALSE) {
  is.dist <- inherits(M, "dist.matrix")
  if (!is.dist) M <- find.canonical.matrix(M) # ensure that M is a suitable matrix, or extract from DSM object

  term.labels <- if (byrow) rownames(M) else colnames(M)
  found <- term %in% term.labels
  if (any(!found) && !skip.missing) stop("target term(s) not found in M: ", paste(term[!found], collapse=", "))
  term <- term[found]
  n.terms <- length(term)
  if (n.terms == 0) return(NULL)
  
  ## unless input is a pre-computed dist.matrix, process vector of lookup terms in moderately sized batches
  if (!is.dist) {
    n.cand <- if (byrow) nrow(M) else ncol(M)
    if (n.cand > 1 && as.double(n.terms) * n.cand > batchsize) {
      items.per.batch <- ceiling(batchsize / n.cand)
      res.list <- lapply(seq(1, n.terms, items.per.batch), function (i.start) {
        i.end <- min(i.start + items.per.batch - 1, n.terms)
        if (verbose) cat(sprintf(" - nearest.neighbours(): terms #%d .. #%d of %d (size = %.1fM)\n", i.start, i.end, n.terms, (i.end-i.start+1) * n.cand / 1e6))
        nearest.neighbours(M, term[i.start:i.end], n=n, drop=FALSE, skip.missing=FALSE, min.distance=min.distance, byrow=byrow, ..., batchsize=Inf, verbose=verbose)
      })
      return(do.call(c, res.list))
    }
  }
  
  if (!is.dist) {
    M.term <- if (byrow) M[term, , drop=FALSE] else M[, term, drop=FALSE]
    M <- dist.matrix(M=M.term, M2=M, byrow=byrow, ...)
    byrow <- TRUE  # dist.matrix computed on the fly is always accessed by row
  }
  similarity <- attr(M, "similarity")
  if (is.null(similarity)) similarity <- FALSE
  
  result <- lapply(term, function (.t) {
    neighbours <- if (byrow) M[.t, ] else M[, .t]
    if (!similarity && min.distance > 0) neighbours <- neighbours[neighbours >= min.distance]
    neighbours <- sort(neighbours, decreasing=similarity)
    head(neighbours, n)
  })
  names(result) <- term
  
  if (drop && length(result) == 1) result[[1]] else result
}
