read.dsm.triplet <- function (filename, freq=FALSE, value.first=FALSE, tokens=FALSE, sep="\t", quote="", nmax=-1, sort=FALSE, encoding=getOption("encoding"), verbose=FALSE) {
  if (verbose) cat(sprintf("Loading DSM triplets from '%s' ... ", filename))
  is.pipe <- grepl("\\|\\s*$", filename, perl=TRUE)
  if (is.pipe) {
    filename <- sub("\\s*\\|\\s*$", "", filename, perl=TRUE)
    fh <- pipe(filename, encoding=encoding) # don't open connection yet, so scan() always converts to UTF-8
  } else {
    fh <- file(filename, encoding=encoding)
  }
  
  if (tokens) {
    if (freq || value.first) warning("freq= and value.first= options are ignored with tokens=TRUE")
    triplets <- scan(fh, what=list(l1=character(), l2=character()), sep=sep, quote=quote, nmax=nmax, quiet=TRUE, multi.line=FALSE, na.strings="", comment.char="")
    triplets$val <- 1 # will automatically be replicated by dsm()
    freq <- TRUE
  } else {
    format <- if (value.first) list(val=double(), l1=character(), l2=character()) else list(l1=character(), l2=character(), val=double())
    triplets <- scan(fh, what=format, sep=sep, quote=quote, nmax=nmax, quiet=TRUE, multi.line=FALSE, na.strings="", comment.char="")
  }
  close(fh)
  
  if (verbose) cat(sprintf("%.2fM %s\n", length(triplets$l1) / 1e6, if (tokens) "tokens" else "items"))

  dsm(target=triplets$l1, feature=triplets$l2, score=triplets$val, raw.freq=freq, sort=sort, verbose=verbose)
}
