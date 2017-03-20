## Validate C code in dsm_score() against pure R implementation
## based on small "hieroglyphs" example matrix

library(wordspace)
library(Matrix)

vec.equal <- function(x, y, tol=1e-12, verbose=TRUE) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) == length(y)) {
    noteq <- x != y # abs(x-y) undefined if both are +Inf or -Inf
    max.diff <- if (any(noteq)) max(abs(x[noteq] - y[noteq])) else 0
    if (verbose && max.diff > max(.Machine$double.eps, .Machine$double.neg.eps)) {
      cat(sprintf("vec.equal: largest difference between vectors = %g\n", max.diff))
    }
    if (max.diff < tol) TRUE else FALSE
  } else {
    if (verbose) cat(sprintf("vec.equal: different vector lengths %d != %d\n", length(x), length(y)))
    FALSE
  }
}

## test data: co-occurrence matrix and corresponding DSM object
M <- DSM_HieroglyphsMatrix
M.dsm <- dsm(M, raw.freq=TRUE)

## marginals and expected frequencies
f1 <- rowSums(M)
f2 <- colSums(M)
N <- sum(M)

M.exp <- outer(f1, f2) / N

stopifnot(vec.equal(f1, M.dsm$rows$f))
stopifnot(vec.equal(f2, M.dsm$cols$f))
stopifnot(vec.equal(N, M.dsm$globals$N))

## validate several score functions and transformations
f <- M # raw co-occurrence frequency
f.dsm <- dsm.score(M.dsm, score="frequency", matrix.only=TRUE)
stopifnot(vec.equal(f, f.dsm))

MI <- pmax(log2(M / M.exp), 0) # sparse pointwise Mutual Information
MI.dsm <- dsm.score(M.dsm, score="MI", matrix.only=TRUE)
stopifnot(vec.equal(MI, MI.dsm))

t <- ifelse(M >= M.exp, (M - M.exp) / sqrt(M), 0) # sparse t-score
t.dsm <- dsm.score(M.dsm, score="t-score", matrix.only=TRUE)
stopifnot(vec.equal(t, t.dsm))

tsqrt.dsm <- dsm.score(M.dsm, score="t-score", transform="root", matrix.only=TRUE) # square root transformation
stopifnot(vec.equal(sqrt(t), tsqrt.dsm))

tlog.dsm <- dsm.score(M.dsm, score="t-score", transform="log", matrix.only=TRUE) # logarithmic transformation
stopifnot(vec.equal(log(t+1), tlog.dsm))

## compare with sparse matrix format
S <- as(M, "sparseMatrix")
stopifnot(dsm.is.canonical(S)$canonical && dsm.is.canonical(S)$sparse)
S.dsm <- dsm(S, raw.freq=TRUE)

f.sparse <- dsm.score(S.dsm, score="frequency", matrix.only=TRUE)
stopifnot(vec.equal(f.dsm, f.sparse))
MI.sparse <- dsm.score(S.dsm, score="MI", matrix.only=TRUE)
stopifnot(vec.equal(MI.dsm, MI.sparse))
t.sparse <- dsm.score(S.dsm, score="t-score", matrix.only=TRUE)
stopifnot(vec.equal(t.dsm, t.sparse))
tsqrt.sparse <- dsm.score(S.dsm, score="t-score", transform="root", matrix.only=TRUE)
stopifnot(vec.equal(tsqrt.dsm, tsqrt.sparse))
tlog.sparse <- dsm.score(S.dsm, score="t-score", transform="log", matrix.only=TRUE)
stopifnot(vec.equal(tsqrt.dsm, tsqrt.sparse))

## column scaling
tsqrt.std <- scale(sqrt(t)) # standardization
tsqrt.std.dsm <- dsm.score(M.dsm, score="t-score", transform="root", scale="standardize", matrix.only=TRUE) 
stopifnot(vec.equal(apply(tsqrt.std.dsm, 2, mean), rep(0, ncol(M))))
stopifnot(vec.equal(apply(tsqrt.std.dsm, 2, sd), rep(1, ncol(M))))
stopifnot(vec.equal(tsqrt.std.dsm, tsqrt.std))

## NB: sparse DSM doesn't allow full standardization unless negative.ok=TRUE (and this is not recommended)

tsqrt.scale <- scale(sqrt(t), center=FALSE) # scale to unit variance without centering
tsqrt.scale.dsm <- dsm.score(M.dsm, score="t-score", transform="root", scale="scale", matrix.only=TRUE)
rms <- function (x) sqrt(sum(x * x) / (length(x) - 1)) # root mean square (cf. ?scale)
stopifnot(vec.equal(apply(tsqrt.scale.dsm, 2, rms), rep(1, ncol(M))))
stopifnot(vec.equal(tsqrt.scale.dsm, tsqrt.scale))

tsqrt.scale.sparse <- dsm.score(S.dsm, score="t-score", transform="root", scale="scale", matrix.only=TRUE) # scale to unit variance without centering
stopifnot(vec.equal(apply(tsqrt.scale.sparse, 2, rms), rep(1, ncol(M))))
stopifnot(vec.equal(tsqrt.scale.sparse, tsqrt.scale))

## row normalization
tsqrt.norm <- sqrt(t) / rowSums(abs(sqrt(t))) # simplest case: Manhattan norm
tsqrt.norm.dsm <- dsm.score(M.dsm, score="t-score", transform="root", normalize=TRUE, method="manhattan", matrix.only=TRUE)
stopifnot(vec.equal(tsqrt.norm.dsm, tsqrt.norm))

tsqrt.norm.sparse <- dsm.score(S.dsm, score="t-score", transform="root", normalize=TRUE, method="manhattan", matrix.only=TRUE)
stopifnot(vec.equal(tsqrt.norm.sparse, tsqrt.norm))

## validate user-defined association measures against built-in C code
my.MI <- function (O11, E11, ...) log2(O11 / E11)
my.tscore <- function (O, E, ...) (O - E) / sqrt(O)
my.tscore.disc <- function (O, E, ...) (O - E) / sqrt(O + 1) # discounted t-score used as non-sparse measure
my.Dice <- function (O11, R1, C1, ...) 2 * O11 / (R1 + C1)
my.ll <- function (O, E, ...) {
  ll <- 2 * (ifelse(O > 0, O * log(O / E), 0) - (O - E))
  ifelse(O >= E, ll, -ll)
}
my.tfidf <- function (O11, cols, ..., K=7) O11 * log(1 / cols$df) # need to annotate df explicitly

TT <- DSM_TermTerm    # dense co-occurrence matrix
TT$cols <- transform(TT$cols, df=((nnzero + 1) / (nrow(TT) + 1)))
TC <- DSM_TermContext # sparse co-occurrence matrix
TC$cols <- transform(TC$cols, df=((nnzero + 1) / (nrow(TC) + 1)))

test.AM <- function(model, AM.name, AM.fun, transform="none", sparse=TRUE) {
  stopifnot(vec.equal(
    dsm.score(model, AM.name, transform=transform, sparse=sparse, negative.ok=TRUE, matrix.only=TRUE),
    dsm.score(model, AM.fun,  transform=transform, sparse=sparse, negative.ok=TRUE, matrix.only=TRUE, batchsize=10)))
}

## dense matrix, sparse AM
test.AM(TT, "MI", my.MI)
test.AM(TT, "t-score", my.tscore)
test.AM(TT, "Dice", my.Dice)
test.AM(TT, "simple-ll", my.ll)
test.AM(TT, "tf.idf", my.tfidf)

test.AM(TT, "simple-ll", my.ll, transform="log")
test.AM(TT, "simple-ll", my.ll, transform="root")
test.AM(TT, "simple-ll", my.ll, transform="sigmoid")

## dense matrix, dense AM
test.AM(TT, "MI", my.MI, sparse=FALSE)
test.AM(TT, "t-score", my.tscore.disc, sparse=FALSE)
test.AM(TT, "Dice", my.Dice, sparse=FALSE)
test.AM(TT, "simple-ll", my.ll, sparse=FALSE)
test.AM(TT, "tf.idf", my.tfidf, sparse=FALSE)

test.AM(TT, "simple-ll", my.ll, transform="log", sparse=FALSE)
test.AM(TT, "simple-ll", my.ll, transform="root", sparse=FALSE)
test.AM(TT, "simple-ll", my.ll, transform="sigmoid", sparse=FALSE)

## sparse matrix, sparse AM
test.AM(TC, "MI", my.MI)
test.AM(TC, "t-score", my.tscore)
test.AM(TC, "Dice", my.Dice)
test.AM(TC, "simple-ll", my.ll)
test.AM(TC, "tf.idf", my.tfidf)

test.AM(TC, "simple-ll", my.ll, transform="log")
test.AM(TC, "simple-ll", my.ll, transform="root")
test.AM(TC, "simple-ll", my.ll, transform="sigmoid")

## sparse matrix, dense AM (casts to dense matrix)
test.AM(TC, "MI", my.MI, sparse=FALSE)
test.AM(TC, "t-score", my.tscore.disc, sparse=FALSE)
test.AM(TC, "Dice", my.Dice, sparse=FALSE)
test.AM(TC, "simple-ll", my.ll, sparse=FALSE)
test.AM(TC, "tf.idf", my.tfidf, sparse=FALSE)

test.AM(TC, "simple-ll", my.ll, transform="log", sparse=FALSE)
test.AM(TC, "simple-ll", my.ll, transform="root", sparse=FALSE)
test.AM(TC, "simple-ll", my.ll, transform="sigmoid", sparse=FALSE)
