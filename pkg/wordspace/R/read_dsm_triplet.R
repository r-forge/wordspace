read.dsm.triplet <- function (filename, freq=FALSE, value.first=FALSE, tokens=FALSE, sep="\t", quote="", nmax=-1, sort=FALSE, encoding=getOption("encoding"), verbose=FALSE) {
  if (verbose) cat(sprintf("Loading triplets from '%s' ... ", filename))
  is.pipe <- grepl("\\|\\s*$", filename, perl=TRUE)
  if (is.pipe) {
    filename <- sub("\\s*\\|\\s*$", "", filename, perl=TRUE)
    fh <- pipe(filename, open="rt", encoding=encoding)
  } else {
    fh <- file(filename, open="rt", encoding=encoding)
  }
  
  if (tokens) {
    if (freq || value.first) warning("freq= and value.first= options are ignored with tokens=TRUE")
    triplets <- scan(fh, what=list(l1=character(), l2=character()), sep=sep, quote=quote, nmax=nmax, quiet=TRUE, multi.line=FALSE, na.strings="")
    triplets$val <- rep(1, length(triplets$l1))
    freq=TRUE
  } else {
    format <- if (value.first) list(val=double(), l1=character(), l2=character()) else list(l1=character(), l2=character(), val=double())
    triplets <- scan(fh, what=format, sep=sep, quote=quote, nmax=nmax, quiet=TRUE, multi.line=FALSE, na.strings="")
  }
  close(fh)
  
  if (verbose) cat(sprintf("%.2fM %s\n", length(triplets$l1) / 1e6, if (tokens) "tokens" else "items"))
  if (verbose) cat("Compiling target & feature terms ... \n")
  
  l1.dict <- unique(triplets$l1) # compile lists of target terms and feature terms
  l2.dict <- unique(triplets$l2)
  
  if (sort) {
    l1.dict <- sort(l1.dict)
    l2.dict <- sort(l2.dict)
  }
  
  row.idx <- match(triplets$l1, l1.dict) # row/column indices of triplets (according to term lists)
  col.idx <- match(triplets$l2, l2.dict)
  
  if (verbose) cat("Building sparse matrix ... ")
  M <- sparseMatrix(i=row.idx, j=col.idx, x=triplets$val) # this may be either M or S
  rm(triplets, row.idx, col.idx) # remove unused large objects
  
  rownames(M) <- l1.dict
  colnames(M) <- l2.dict
  if (verbose) cat(sprintf("%d x %d\n", nrow(M), ncol(M)))
  
  if (freq) {
    if (verbose) cat("Computing marginal frequencies ... \n")
    is.nzero <- M > 0
    # NB: rowNorms(x, "manhattan") is a memory-efficient alternative to rowSums for non-negative matrix
    rows <- data.frame(term=l1.dict, f=rowNorms(M, "manhattan"), nnzero=rowSums(is.nzero), stringsAsFactors=FALSE)
    cols <- data.frame(term=l2.dict, f=colNorms(M, "manhattan"), nnzero=colSums(is.nzero), stringsAsFactors=FALSE)
    
    rm(is.nzero) # remove unused large object
    N <- sum(M)
    globals <- list(N=N)
    dsm <- list(M=M, rows=rows, cols=cols, N=N, globals=globals, locked=FALSE)
  } else {
    if (verbose) cat("Constructing pre-scored DSM object ... \n")
    is.nzero <- M != 0
    rows <- data.frame(term=l1.dict, nnzero=rowSums(is.nzero), stringsAsFactors=FALSE)
    cols <- data.frame(term=l2.dict, nnzero=colSums(is.nzero), stringsAsFactors=FALSE)
    rm(is.nzero) # remove unused large object
    dsm <- list(S=M, rows=rows, cols=cols, globals=list(), locked=FALSE)
  }

  if (verbose) cat(sprintf("%d x %d matrix with %d nonzero entries (fill rate = %.2f%%)\n", nrow(M), ncol(M), nnzero(M), 100 * nnzero(M) / prod(dim(M))))

  class(dsm) <- c("dsm", "list")
  return(dsm)
}