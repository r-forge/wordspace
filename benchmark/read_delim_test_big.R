##
## Benchmark different ways of reading TAB-delimited file with 60M rows
## (requires external data set)
##
source("benchmark_utils.R")

library(readr)   # popular new package, but heavy dependencies
library(iotools) # claimed to be fast for reading simple CSV formats

## use a .gz-compressed file to test processing speed (.bz2 has considerably more overhead)
filename <- "~/Writings/Conferences/2014/COLING2014/wordspace/dm-demo/typedm_w_w.txt.gz" # 462 MiB compressed, 1559 MiB uncompressed
n.cells <- 60337357 * 3 # number of cells in the table

for (pass in 1:2) {
cat(sprintf("----- pass #%d -----\n", pass))
L1 <- benchmark(
	vn1 <<- read.delim(filename, header=FALSE, colClasses=c("character", "character", "factor"), col.names=c("w1", "w2", "score"), sep="\t", quote="", comment.char="", fileEncoding="utf8"), 
	name="read.delim()", n.ops=n.cells)
## MOPS = number of cells read per second

cat("object size = ", format(object.size(vn1), units="MB"), "\n") # object size = 1247 MiB = 1.2 GiB

L2 <- benchmark(
	vn2 <<- read_delim(filename, "\t", col_names=c("w1", "w2", "score"), col_types="ccd", locale=locale(encoding="utf8"), quote="", comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE),
	name="readr::read_delim()", n.ops=n.cells)
all.equal(vn1, vn2, check.attributes=FALSE)
rm(vn2)

L3 <- benchmark(
	vn3 <<- structure(read.delim.raw(file(filename, encoding="utf8"), header=FALSE, colClasses=c("character", "character", "numeric"), sep="\t", quote=""), names=c("w1", "w2", "score")),
	name="iotools::read.delim.raw()", n.ops=n.cells)
all.equal(vn1, vn3, check.attributes=FALSE)
rm(vn3)
rm(vn1)

L <- do.call(rbind, list(L1, L2, L3))
print(L)
}
