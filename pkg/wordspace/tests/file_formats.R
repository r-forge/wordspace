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
