## Test exporting and importing a DSM matrix in different file formats

library(wordspace)

## export word2vec text format
word2vec_hiero <- system.file("extdata", "word2vec_hiero.txt", package="wordspace", mustWork=TRUE)
M <- DSM_HieroglyphsMatrix / 100
gold.lines <- readLines(word2vec_hiero)

fn <- tempfile(fileext=".txt")
write.dsm.matrix(M, fn, format="word2vec", round=1)
if (!isTRUE(all.equal(gold.lines, readLines(fn)))) stop("error in word2vec format (with round=1)")

fn <- tempfile(fileext=".txt.gz")
fh <- gzfile(fn, encoding="UTF-8")
write.dsm.matrix(M, fh, format="word2vec", round=1)
if (!isTRUE(all.equal(gold.lines, readLines(fn)))) stop("error in word2vec.gz format (with round=1)")

## import word2vec text format
M2 <- read.dsm.matrix(word2vec_hiero, format="word2vec")
colnames(M2) <- colnames(M) # not in word2vec format
if (!isTRUE(all.equal(round(M, 1), M2))) stop("error reading word2vec format")

## round trip for larger matrix
BNC <- with(DSM_VerbNounTriples_BNC,
            dsm(target=noun, feature=verb, score=f, raw.freq=TRUE))
BNC <- dsm.score(BNC, score="simple-ll", transform="log", normalize=TRUE, update.nnzero=TRUE)
BNC100 <- dsm.projection(BNC, n=100, method="svd")

fn <- tempfile(fileext=".txt")
write.dsm.matrix(BNC100, fn, format="word2vec", round=3)
Read <- read.dsm.matrix(fn, format="word2vec")

Gold <- round(BNC100[, ], 3) # round the original matrix,
colnames(Gold) <- NULL       # removing extra attributes and column names
stopifnot(all.equal(Gold, Read))
