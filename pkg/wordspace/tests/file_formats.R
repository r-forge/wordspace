## Test reading data sets in different file formats, based on the example DSM_TermContext matrix

library(wordspace)

## example data: DSM_TermContext
TC <- DSM_TermContext        # order rows and columns alphabetically for comparison with input files
idxR <- order(TC$rows$term)
idxC <- order(TC$cols$term)
TC$M <- TC$M[idxR, idxC]
TC$rows <- TC$rows[idxR, ]
TC$cols <- TC$cols[idxC, ]
invisible(check.dsm(TC, validate=TRUE))
head(TC, Inf)

## Co-occurrence matrix in triplet representation
triplet.file <- system.file("extdata", "term_context_triplets.gz", package="wordspace", mustWork=TRUE)
TC.triplet <- read.dsm.triplet(triplet.file, freq=TRUE, sort=TRUE, encoding="ascii")
invisible(check.dsm(TC.triplet, validate=TRUE))

res <- all.equal(TC.triplet$M, TC$M) # correct marginal frequencies cannot be obtained from triplet file
if (!isTRUE(res)) {
  cat("Mismatch for matrix loaded from triplet file:\n")
  cat(paste0("\t", res, "\n"))
  stop("test failed")
}

## UCS export format
ucs.file <- system.file("extdata", "term_context_ucs", package="wordspace", mustWork=TRUE)
TC.ucs <- read.dsm.ucs(ucs.file, encoding="ascii")
invisible(check.dsm(TC.ucs, validate=TRUE))

res <- all.equal(TC.ucs, TC, check.attributes=FALSE) # recursive all.equal seems broken for att's in list-like structures
if (!isTRUE(res)) {
  cat("Mismatch for DSM loaded from UCS export format:\n")
  cat(paste0("\t", res, "\n"))
  stop("test failed")
}

## Check that non-ASCII characters in different encodings are read correctly
Encode.ref <- c("Test", "\u{0164}\u{00E9}\u{015F}t") # expected rownames of DSM

test.encoding <- function (filename, encoding, ref=Encode.ref, force=FALSE) {
  if (!force && !(encoding %in% iconvlist())) {
    warning(sprintf("the %s character encoding is not supported on this platform", encoding))
  } else {
    model <- read.dsm.triplet(system.file("extdata", filename, package="wordspace", mustWork=TRUE), encoding=encoding, tokens=TRUE, sort=FALSE)
    if (!all(ref == model$rows$term)) stop(sprintf("failed to load %s triplet file", encoding))
  }
}

test.encoding("tokens_utf8.txt", "UTF-8", force=TRUE) # should work on all platforms
test.encoding("tokens_latin2.txt", "ISO-8859-2")
test.encoding("tokens_utf16.txt", "UTF-16")
