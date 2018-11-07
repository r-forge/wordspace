## Test computation of centroids as context vectors

library(wordspace)

vec.compare <- function (x, y, name="vector comparison", normalize=TRUE, tol=1e-10, verbose=TRUE) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (any(is.na(x))) {
    cat(sprintf("%s: missing values in first vector\n", name))
    invisible(FALSE)
  } else if (any(is.na(y))) {
    cat(sprintf("%s: missing values in second vector\n", name))
    invisible(FALSE)
  } else if (length(x) != length(y)) {
    if (verbose) cat(sprintf("%s: different vector lengths %d != %d\n", name, length(x), length(y)))
    invisible(FALSE)
  } else {
    if (normalize) {
      x <- x / sqrt(sum(x ^ 2))
      y <- y / sqrt(sum(y ^ 2))
    }
    max.diff <- max(abs(x - y))
    if (max.diff < tol) {
      invisible(TRUE)
    } else {
      if (verbose) cat(sprintf("%s: largest difference between vectors = %g exceeds tolerance limit\n", name, max.diff))
      invisible(FALSE)
    }
  }
}


## test case: centroid vector for the document "cat cat dog cause cause"
M <- DSM_TermTermMatrix
x.ref <- 2 * M["cat", ] + M["dog", ] + 2 * M["cause", ]

doc <- "cat cat dog cause cause" # centroid of document as string
x1 <- context.vectors(M, doc)
stopifnot( vec.compare(x1, x.ref, "centroid of document string") )

doc <- c("cat", "cat", "dog", "cause", "cause") # centroid of document as list of tokens
x2 <- context.vectors(M, list(doc))
stopifnot( vec.compare(x2, x.ref, "centroid of vector of tokens") )

doc <- c("cat", "cause", "cause", "dog", "cat") # should be independent of ordering
x3 <- context.vectors(M, list(doc))
stopifnot( vec.compare(x3, x.ref, "centroid of vector of tokens (reordered)") )

doc <- c(cat=1, cat=1, dog=1, cause=1, cause=1) # centroid of document as bag of words
x4 <- context.vectors(M, list(doc))
stopifnot( vec.compare(x4, x.ref, "centroid of bag of words (tokens)") )

doc <- c(cat=1, cat=1, dog=1, cause=1, cause=1) # centroid of document as bag of words
x4 <- context.vectors(M, list(doc))
stopifnot( vec.compare(x4, x.ref, "centroid of bag of words (tokens)") )

doc <- c(cat=2, dog=1, cause=2) # with aggregated frequency counts
x5 <- context.vectors(M, list(doc))
stopifnot( vec.compare(x5, x.ref, "centroid of bag of words (types)") )

doc <- c(cat=.2, dog=.1, cause=.2) # check that non-integer values work correctly
x6 <- context.vectors(M, list(doc))
stopifnot( vec.compare(x6, x.ref, "centroid of bag of words (weighted)") )

stopifnot( vec.compare(x5, x6, "scaling of bag-of-words weights", normalize=FALSE))
