read.dsm.triplet <- function (filename, freq=FALSE, value.first=FALSE, tokens=FALSE, sep="\t", quote="", nmax=-1, sort=FALSE, encoding=getOption("encoding"), verbose=FALSE) {
  if (verbose) cat(sprintf("Loading DSM triplets from '%s' ... ", filename))
  is.pipe <- grepl("\\|\\s*$", filename, perl=TRUE)
  if (is.pipe) {
    filename <- sub("\\s*\\|\\s*$", "", filename, perl=TRUE)
    fh <- pipe(filename, encoding=encoding)
  } else {
    fh <- file(filename, encoding=encoding)
  }
  
  if (tokens) {
    if (freq || value.first) warning("freq= and value.first= options are ignored with tokens=TRUE")
    triplets <- structure(
      read.delim.raw(fh, sep="\t", quote="", header=FALSE, colClasses=c("character", "character"), nrows=nmax, strict=TRUE),
      names=c("l1", "l2"))
    ## iotools::read.delim.raw is better than readr::read_delim:
    ##  - readr has many "expensive" dependencies (esp. Boost in package 'BH')
    ##  - readr doesn't support the default "native.enc" encoding and cannot read from all types of connections
    ##  - iotools is slightly faster and leaner (less memory overhead) than readr
    ## triplets <- read_delim(fh, "\t", col_names=c("l1", "l2"), col_types="cc", locale=locale(encoding=encoding), n_max=nmax, quote="", comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE)
    triplets$val <- 1
    freq <- TRUE
  } else {
    col.types <- if (value.first) c("numeric", "character", "character") else c("character", "character", "numeric")
    col.names <- if (value.first) c("val", "l1", "l2") else c("l1", "l2", "val")
    triplets <- structure(
      read.delim.raw(fh, sep="\t", quote="", header=FALSE, colClasses=col.types, nrows=nmax, strict=TRUE),
      names=col.names)
    ## triplets <- read_delim(fh, "\t", col_names=col.names, col_types=col.types, locale=locale(encoding=encoding), n_max=nmax, quote="", comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE)
  }
  ## close(fh) not needed, because read.delim.raw automatically opens and closes the connection
  
  if (verbose) cat(sprintf("%.2fM %s\n", length(triplets$l1) / 1e6, if (tokens) "tokens" else "items"))

  dsm(target=triplets$l1, feature=triplets$l2, score=triplets$val, raw.freq=freq, sort=sort, verbose=verbose)
}
