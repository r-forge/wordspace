pair.distances <- function (w1, w2, M, ..., rank=c("none", "fwd", "bwd", "avg"), transform=NULL, avg.method=c("arithmetic", "geometric", "harmonic"), batchsize=10e6, verbose=FALSE) {
  rank <- match.arg(rank)
  avg.method <- match.arg(avg.method)
  if (!is.null(transform) && !is.function(transform)) stop("transform= must be a vectorized function expecting a single argument")
  w1 <- as.character(w1)
  w2 <- as.character(w2)
  stopifnot(length(w1) == length(w2))
  n <- length(w1)
  types1 <- unique(w1)
  types2 <- if (rank == "none") unique(w2) else NULL
  n.types1 <- as.double(length(types1))
  n.types2 <- as.double(length(types2))
  M <- find.canonical.matrix(M) # ensure matrix is in canonical format, or extract from DSM object
  
  if (rank != "none") {
    ## dispatch to specialised function for computing neighbour ranks
    if (rank == "fwd") {
      r <- pair.ranks(w1, w2, M, ..., batchsize=batchsize, verbose=verbose)
    } else if (rank == "bwd") {
      r <- pair.ranks(w2, w1, M, ..., batchsize=batchsize, verbose=verbose)
    } else {
      r1 <- pair.ranks(w1, w2, M, ..., batchsize=batchsize, verbose=verbose)
      r2 <- pair.ranks(w2, w1, M, ..., batchsize=batchsize, verbose=verbose)
      r <- switch(avg.method,
                  arithmetic = (r1 + r2) / 2,
                  geometric = sqrt(r1 * r2),
                  harmonic = ifelse(is.finite(pmax(r1, r2)), 2 * r1 * r2 / (r1 + r2), Inf))
    }
    if (!is.null(transform)) transform(r) else r
  } else {
    ## if there are too many distinct types (such that intermediate distance matrix would have > chunksize elements), 
    ## partition input recursively and combine result vectors
    n.elements <- n.types1 * n.types2 # size of distance matrix to be computed
    split.batch <- n.elements > batchsize && n >= 4
    if (verbose) cat(sprintf("%s- pair.distances(): %d pairs, %d x %d types = %.1fM elements %s\n", paste(rep(" ", verbose), collapse=""), n, length(types1), length(types2), n.elements/1e6, if (split.batch) "" else "***"))
    if (split.batch) {
      pivot <- floor(n/2)
      verbose.val <- if (verbose) verbose + 2 else FALSE
      res <- c(
        pair.distances(w1[1:pivot], w2[1:pivot], M, ..., batchsize=batchsize, verbose=verbose.val),
        pair.distances(w1[(pivot+1):n], w2[(pivot+1):n], M, ..., batchsize=batchsize, verbose=verbose.val)
      )
    }
    else {  
      distances <- dist.matrix(M, byrow=TRUE, terms=types1, terms2=types2, skip.missing=TRUE, ...)
      found1 <- rownames(distances)
      found2 <- colnames(distances)
      is.similarity <- attr(distances, "similarity")
      miss.val <- if (!is.null(is.similarity) && is.similarity) -Inf else Inf
      res <- rep(miss.val, n)
      is.known <- w1 %in% found1 & w2 %in% found2 
      res[is.known] <- distances[cbind(w1[is.known], w2[is.known])]
      names(res) <- paste(w1, w2, sep="/")
    }
    if (!is.null(transform)) transform(res) else res
  }
}

## need specialised implemenation for neighbour ranks, which loops over w1 types instead of word pairs,
## in order to avoid repeated expensive computation of all-neighbour rankings
pair.ranks <- 
function (w1, w2, M, ..., batchsize=10e6, verbose=FALSE) {
  w1 <- as.character(w1)
  w2 <- as.character(w2)
  stopifnot(length(w1) == length(w2))
  res <- rep(Inf, length(w1))

  known.words <- rownames(M)
  is.known <- w1 %in% known.words & w2 %in% known.words # word pairs found in DSM
  w1.types <- sort(unique(w1[is.known])) # list of w1 types found in DSM
  n.types <- length(w1.types)
  items.per.batch <- ceiling(batchsize / nrow(M)) # number of w1 types we can process in a single batch
  
  for (i.start in seq(1, n.types, items.per.batch)) {
    i.end <- min(i.start + items.per.batch - 1, n.types)
    if (verbose) cat(sprintf(" - pair.ranks(): types #%d .. #%d of %d (%s .. %s)\n", i.start, i.end, n.types, w1.types[i.start], w1.types[i.end]))
    batch.types <- w1.types[i.start:i.end]
    distances <- dist.matrix(M, byrow=TRUE, terms=batch.types, terms2=NULL, skip.missing=FALSE, ...)
    is.similarity <- attr(distances, "similarity")
    if (!is.null(is.similarity) && is.similarity) distances <- -distances
    ranks <- t(apply(distances, 1, function (x) rank(x, ties.method="min")))
    idx <- (w1 %in% batch.types) & is.known
    res[idx] <- ranks[cbind(w1[idx], w2[idx])] - 1 # adjust for w1 as its own neighbour
  }
  names(res) <- paste(w1, w2, sep="/")
  res
}
