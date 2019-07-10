## Test whether DSM object has expected structure

library(wordspace)

## Verify that reconstructed built-in DSMs have exactly the same structure and attributes
reconstruct.dsm <- function (Orig, name=sprintf("DSM object '%s'", deparse(substitute(Orig)))) {
  ## extract relevant information in a "minimally similar" way
  if (dsm.is.canonical(Orig$M)$sparse) {
    M <- Orig$M
    attr(M, "nonneg") <- FALSE
  } else {
    M <- matrix(Orig$M, nrow(Orig$M), ncol(Orig$M), dimnames=dimnames(Orig$M))
  }
  row.f <- Orig$rows[sample(nrow(M)), c("term", "f")] # ordering of data frame shouldn't matter
  col.f <- Orig$cols[sample(ncol(M)), c("term", "f")]
  N <- Orig$globals$N
  New <- dsm(M=M, rowinfo=row.f, colinfo=col.f, N=N, raw.freq=TRUE)
  if (!isTRUE(all.equal(Orig, New))) stop(sprintf("structural difference in %s", name))
}

reconstruct.dsm(DSM_TermTerm)
reconstruct.dsm(DSM_TermContext)
