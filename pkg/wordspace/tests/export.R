## Test exporting a DSM in different file formats

library(wordspace)

## word2vec text format
M <- DSM_HieroglyphsMatrix / 100
gold.lines <- readLines(system.file("extdata", "word2vec_hiero.txt", package="wordspace", mustWork=TRUE))

fn <- tempfile(fileext=".txt")
write.dsm.matrix(M, fn, format="word2vec", round=1)
if (!isTRUE(all.equal(gold.lines, readLines(fn)))) stop("error in word2vec format (with round=1)")

fn <- tempfile(fileext=".txt.gz")
fh <- gzfile(fn, encoding="UTF-8")
write.dsm.matrix(M, fh, format="word2vec", round=1)
if (!isTRUE(all.equal(gold.lines, readLines(gzfile(fn))))) stop("error in word2vec.gz format (with round=1)")
