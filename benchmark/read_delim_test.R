##
## Benchmark different ways of reading files
##
source("benchmark_utils.R")

library(readr)   # popular new package, but heavy dependencies
library(iotools) # claimed to be fast for reading simple CSV formats

filename <- "../datasets/orig/bnc_cc_vn_triples.txt.gz" # only 10 MiB compressed, 61 MiB uncompressed
n.cells <- 2236106 * 5 # number of cells in the table

for (pass in 1:2) {
cat(sprintf("----- pass #%d -----\n", pass))
L1 <- benchmark(
	vn1 <<- read.delim(filename, header=FALSE, colClasses=c("numeric", "character", "factor", "character", "factor"), col.names=c("f", "n", "rel", "v", "mode"), sep="\t", quote="", comment.char="", fileEncoding="utf8"), 
	name="read.delim() w/ factors", n.ops=n.cells)
## MOPS = number of cells read per second

cat("object size = ", format(object.size(vn1), units="MB"), "\n") # final object size = 70.5 MiB

L2 <- benchmark(
	vn2 <<- read_delim(filename, "\t", col_names=c("f", "n", "rel", "v", "mode"), col_types=cols(f=col_double(), n=col_character(), rel=col_factor(c("obj", "subj")), v=col_character(), mode=col_factor(c("spoken", "written"))), locale=locale(encoding="utf8"), quote="", comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE),
	name="readr::read_delim() w/ factors", n.ops=n.cells)
all.equal(vn1, vn2, check.attributes=FALSE)
rm(vn2)

L3 <- benchmark({
	vn3 <<- structure(read.delim.raw(file(filename, encoding="utf8"), header=FALSE, colClasses=c("numeric", "character", "character", "character", "character"), sep="\t", quote=""), names=c("f", "n", "rel", "v", "mode"))
	vn3$rel <- as.factor(vn3$rel)
	vn3$mode <- as.factor(vn3$mode)
	}, name="iotools::read.delim.raw() w/ factors", n.ops=n.cells)
all.equal(vn1, vn3, check.attributes=FALSE)
rm(vn3)
rm(vn1)

L4 <- benchmark(
	vn1 <<- read.delim(filename, header=FALSE, colClasses=c("numeric", "character", "character", "character", "character"), col.names=c("f", "n", "rel", "v", "mode"), sep="\t", quote="", comment.char="", fileEncoding="utf8"), 
	name="read.delim()", n.ops=n.cells)

cat("object size = ", format(object.size(vn1), units="MB"), "\n") # final object size = 87.6 MiB

L5 <- benchmark(
	vn2 <<- read_delim(filename, "\t", col_names=c("f", "n", "rel", "v", "mode"), col_types="dcccc", locale=locale(encoding="utf8"), quote="", comment="", na=character(), escape_double=FALSE, escape_backslash=FALSE, progress=FALSE),
	name="readr::read_delim()", n.ops=n.cells)
all.equal(vn1, vn2, check.attributes=FALSE)
rm(vn2)

L6 <- benchmark(
	vn3 <<- structure(read.delim.raw(file(filename, encoding="utf8"), header=FALSE, colClasses=c("numeric", "character", "character", "character", "character"), sep="\t", quote=""), names=c("f", "n", "rel", "v", "mode")),
	name="iotools::read.delim.raw()", n.ops=n.cells)
all.equal(vn1, vn3, check.attributes=FALSE)
rm(vn3)
rm(vn1)

L <- do.call(rbind, list(L1, L2, L3, L4, L5, L6))
print(L)
}
