read.dsm.triplet <- function (filename, freq=FALSE, value.first=FALSE, tokens=FALSE, sep="\t", quote="", nmax=-1, sort=FALSE, encoding=getOption("encoding"), verbose=FALSE) {
  if (verbose) cat(sprintf("Loading DSM triplets from '%s' ... ", filename))
  is.pipe <- grepl("\\|\\s*$", filename, perl=TRUE)
  if (is.pipe) {
    filename <- sub("\\s*\\|\\s*$", "", filename, perl=TRUE)
    fh <- pipe(filename)
  } else {
    fh <- filename
  }
  
  if (tokens) {
    if (freq || value.first) warning("freq= and value.first= options are ignored with tokens=TRUE")
    triplets <- read_delim(fh, "\t", col_names=c("l1", "l2"), col_types="cc", locale=locale(encoding=encoding), n_max=nmax, quote="", comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE)
    triplets$val <- 1
    freq <- TRUE
  } else {
    col.types <- if (value.first) "dcc" else "ccd"
    col.names <- if (value.first) c("val", "l1", "l2") else c("l1", "l2", "val")
    triplets <- read_delim(fh, "\t", col_names=col.names, col_types=col.types, locale=locale(encoding=encoding), n_max=nmax, quote="", comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE)
  }
  if (is.pipe) close(fh)
  
  if (verbose) cat(sprintf("%.2fM %s\n", length(triplets$l1) / 1e6, if (tokens) "tokens" else "items"))

  dsm(target=triplets$l1, feature=triplets$l2, score=triplets$val, raw.freq=freq, sort=sort, verbose=verbose)
}
