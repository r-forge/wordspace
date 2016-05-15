## replacement for iotools::read.delim.raw, which re-encodes input file after reading
## - most options are fixed and column names and classes have to be specified by the user
## - fh must be a connection that can be properly opened in binary mode (with automatic decompression if necessary)
## - colClasses should be a named vector to set column names of the table
## - unless encoding is "" or "native.enc", input will be converted to UTF-8 and marked as such
.read.delim.raw <- function (fh, colClasses, nrows=-1, encoding=getOption("encoding")) {
  is.native <- encoding == "" || encoding == "native.enc"
  is.utf8 <- any(grepl("^utf-?8$", encoding, ignore.case=TRUE))
  bindata <- readAsRaw(fh) # read entire file as raw vector (may want to make block size n= larger than default of 64k bytes)
  if (!is.native && !is.utf8) {
    ## use iconv to recode the raw vector into UTF-8
    bindata <- iconv(list(bindata), from=encoding, to="UTF-8", toRaw=TRUE)[[1]]
  }
  ## now parse table in raw vector into a data frame
  res <- dstrsplit(bindata, colClasses, sep="\t", quote="", strict=TRUE, nrows=nrows)
  if (!is.native) {
    ## mark all character variables as UTF-8 (unless read with native encoding)
    for (i in seq_along(colClasses)) {
      if (colClasses[i] == "character") Encoding(res[[i]]) <- "UTF-8"
    }
  }
  res
}

read.dsm.triplet <- function (filename, freq=FALSE, value.first=FALSE, tokens=FALSE, sep="\t", quote="", nmax=-1, sort=FALSE, encoding=getOption("encoding"), verbose=FALSE) {
  if (verbose) cat(sprintf("Loading DSM triplets from '%s' ... ", filename))
  is.pipe <- grepl("\\|\\s*$", filename, perl=TRUE)
  if (is.pipe) {
    filename <- sub("\\s*\\|\\s*$", "", filename, perl=TRUE)
    fh <- pipe(filename) # don't set encoding= parameter because .read.delim.raw will read binary data and convert later
  } else {
    fh <- file(filename) # don't set encoding= parameter
  }
  
  if (tokens) {
    if (freq || value.first) warning("freq= and value.first= options are ignored with tokens=TRUE")
    triplets <- .read.delim.raw(fh, colClasses=c(l1="character", l2="character"), nrows=nmax, encoding=encoding)
    ## I prefer iotools::read.delim.raw over readr::read_delim because:
    ##  - readr has many "expensive" dependencies (esp. Boost in package 'BH')
    ##  - readr doesn't support the default "native.enc" encoding and cannot read from all types of connections
    ##  - iotools is slightly faster and leaner (less memory overhead) than readr
    ##  - unfortunately, read.delim.raw doesn't convert character encodings at all, but we work around this with .read.delim.raw above
    ## Alternative version using readr::read_delim:
    ##   triplets <- read_delim(fh, "\t", col_names=c("l1", "l2"), col_types="cc", locale=locale(encoding=encoding), n_max=nmax, quote="", comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE)
    triplets$val <- 1
    freq <- TRUE
  } else {
    col.types <- if (value.first) c(val="numeric", l1="character", l2="character") else c(l1="character", l2="character", val="numeric")
    triplets <- .read.delim.raw(fh, colClasses=col.types, nrows=nmax, encoding=encoding)
    ## Alternative version using readr::read_delim:
    ##   col.types <- if (value.first) "dcc" else "ccd"
    ##   col.names <- if (value.first) c("val", "l1", "l2") else c("l1", "l2", "val")
    ##   triplets <- read_delim(fh, "\t", col_names=col.names, col_types=col.types, locale=locale(encoding=encoding), n_max=nmax, quote="", comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE)
  }
  ## close(fh) not needed because read.delim.raw automatically opens and closes the connection
  
  if (verbose) cat(sprintf("%.2fM %s\n", length(triplets$l1) / 1e6, if (tokens) "tokens" else "items"))

  dsm(target=triplets$l1, feature=triplets$l2, score=triplets$val, raw.freq=freq, sort=sort, verbose=verbose)
}
