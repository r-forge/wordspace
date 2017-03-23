## Validate pair.distances() implementation

library(wordspace)

qw <- function (x) unlist(strsplit(x, "\\s+", perl=TRUE)) # Perl's qw()

M <- DSM_TermContextMatrix                # so Manhattan distances are simple
DM <- dist.matrix(M, method="manhattan")  # test against elements from distance matrix

w1 <-  qw("cat  cat  dog  rule  dog    rule")
w2 <-  qw("dog  pig  cat  cat   cause  pig")
d.exp <- c(24, Inf,  24,  Inf,  32,    Inf) # expected distances,
r.exp <- c(2,  Inf,  1,   Inf,  6,     Inf) # forward ranks,
b.exp <- c(1,  Inf,  2,   Inf,  4,     Inf) # and backward ranks

## compute pair distances from vectors
d1 <- pair.distances(w1, w2, M, method="manhattan")
stopifnot(all(d1 == d.exp))
stopifnot(all(names(d1) == paste(w1, w2, sep="/"))) # check correct labels

r1 <- pair.distances(w1, w2, M, method="manhattan", rank="fwd")
stopifnot(all(r1 == r.exp))
b1 <- pair.distances(w1, w2, M, method="manhattan", rank="bwd")
stopifnot(all(b1 == b.exp))

## compute pair distances by direct lookup in pre-computed matrix
d2 <- pair.distances(w1, w2, DM, method="euclidean") # method will be ignored!
stopifnot(all(d2 == d.exp))
stopifnot(all(names(d2) == paste(w1, w2, sep="/")))

r2 <- pair.distances(w1, w2, DM, method="euclidean", rank="fwd")
stopifnot(all(r2 == r.exp))
b2 <- pair.distances(w1, w2, DM, method="euclidean", rank="bwd")
stopifnot(all(b2 == b.exp))

## using pair.distances() to look up co-occurrence frequencies (with a fake dist.matrix)
DM <- as.distmat(DSM_TermContextMatrix, similarity=TRUE)
w1 <-  qw("cat    cat   time  time       yyz   cause      yyz")
w2 <-  qw("Feral  Rush  Boat  Back_Pain  Kant  Back_Pain  Rush")
d.exp <- c(7,     -Inf, 2,    0,         -Inf, 6,         -Inf)
r.exp <- c(3,     Inf,  1,    Inf,       Inf,  1,         Inf)
b.exp <- c(2,     Inf,  2,    Inf,       Inf,  1,         Inf)

d3 <- pair.distances(w1, w2, DM)
stopifnot(all(d3 == d.exp))
stopifnot(all(names(d3) == paste(w1, w2, sep="/")))

r3 <- pair.distances(w1, w2, DM, rank="fwd")
stopifnot(all(r3 == r.exp))
b3 <- pair.distances(w1, w2, DM, rank="bwd")
stopifnot(all(b3 == b.exp))

