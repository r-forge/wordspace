## Test various convenience methods for DSM objects and related functionality

library(wordspace)

dsm1 <- DSM_TermContext
dsm2 <- subset(DSM_TermTerm, f < 100000) # reduce to 5x7 matrix
dsm2 <- dsm.score(dsm2, function(O, E, ...) round(log2(O/E), 2))


## dim() of DSM object
stopifnot(all.equal(dim(dsm1), c(7, 7)))
stopifnot(all.equal(nrow(dsm2), 5))
stopifnot(all.equal(ncol(dsm2), 7))


## reading and setting dimnames
stopifnot(all.equal(dimnames(dsm1), dimnames(dsm1$M)))
stopifnot(all.equal(dimnames(dsm2), dimnames(dsm2$M)))
rn <- c("cat", "dog", "animal", "reason", "cause") # verify row/column names directly
cn <- c("breed", "tail", "feed", "kill", "important", "explain", "likely")
stopifnot(all.equal(rownames(dsm2), rn))
stopifnot(all.equal(colnames(dsm2), cn))

dsm2n <- dsm2 # modification of row/column names
rownames(dsm2n)[3] <- "pet"
stopifnot(all.equal(dimnames(dsm2), list(rn, cn))) # must not affect original object
stopifnot(all.equal(colnames(dsm2n), cn)) # no changes to column names yet
rn2 <- rn; rn2[3] <- "pet" 
stopifnot(all.equal(rownames(dsm2n), rn2)) # but one element of rownames has been changed
cn2 <- LETTERS[1:7] # replace all column names
colnames(dsm2n) <- cn2
stopifnot(all.equal(dimnames(dsm2), list(rn, cn)))
stopifnot(all.equal(rownames(dsm2n), rn2))
stopifnot(all.equal(colnames(dsm2n), cn2))

invisible(check.dsm(dsm2n, validate=TRUE)) # checks that rows$term and cols$term have been updated


## extraction of co-occurrence or score matrix
stopifnot(all.equal(as.matrix(dsm1), dsm1$M))
stopifnot(all.equal(as.matrix(dsm2), dsm2$S)) # automatic selection
stopifnot(all.equal(as.matrix(dsm2, what="M"), dsm2$M))
stopifnot(all.equal(as.matrix(dsm2, what="S"), dsm2$S))


## transposition
dsm2t <- t(dsm2)
stopifnot(nrow(dsm2t) == ncol(dsm2))
stopifnot(ncol(dsm2t) == nrow(dsm2))
stopifnot(all.equal(dsm2t$M, t(dsm2$M)))
stopifnot(all.equal(dsm2t$S, t(dsm2$S)))


## efficient checks for non-negativity asnd nonzero count
R.signcount <- function (x) {
  if (is(x, "Matrix")) x <- as.matrix(x) # just test on small sparse matrix examples
  c(pos=sum(x > 0), zero=sum(x == 0), neg=sum(x < 0))
}
R.nonneg <- function (x) {
  if (is(x, "Matrix")) x <- as.matrix(x)
  !any(x < 0)
}
R.nnzero <- function (x) {
  if (is(x, "Matrix")) x <- as.matrix(x)
  sum(x != 0)
}

x <- round(rnorm(1e6), .1) # test on a large numeric vector
y <- as.integer(round(x))  # and integer counterpart
stopifnot(all.equal(signcount(x), R.signcount(x)))
stopifnot(sum(signcount(x)) == length(x))
stopifnot(all.equal(signcount(y), R.signcount(y)))
stopifnot(sum(signcount(y)) == length(y))
stopifnot(all.equal(signcount(x, "nonneg"), R.nonneg(x)))
stopifnot(all.equal(signcount(y, "nonneg"), R.nonneg(y)))
stopifnot(all.equal(signcount(x, "nnzero"), R.nnzero(x)))
stopifnot(all.equal(signcount(y, "nnzero"), R.nnzero(y)))

stopifnot(all.equal(signcount(DSM_TermTermMatrix), R.signcount(DSM_TermTermMatrix))) # dense numeric matrix
M <- Matrix(DSM_HieroglyphsMatrix) # and a dense Matrix object
stopifnot(all.equal(signcount(M), R.signcount(M)))
stopifnot(sum(signcount(M)) == prod(dim(M)))
stopifnot(all.equal(signcount(M, "nonneg"), R.nonneg(M)))
stopifnot(all.equal(signcount(M, "nnzero"), R.nnzero(M)))

stopifnot(all.equal(signcount(DSM_TermContextMatrix), R.signcount(DSM_TermContextMatrix))) # sparse matrix (dgCMatrix)
M <- as(DSM_TermContextMatrix, "dgTMatrix") # triplet representation is not supported
stopifnot(is(try(signcount(M), silent=TRUE), "try-error"))
M <- as(as.matrix(M), "dgRMatrix") # sparse matrix (dgRMatrix), has only minimal support in Matrix package
stopifnot(all.equal(signcount(M), R.signcount(M)))
stopifnot(all.equal(signcount(M, "nonneg"), R.nonneg(M)))
stopifnot(all.equal(signcount(M, "nnzero"), R.nnzero(M)))
