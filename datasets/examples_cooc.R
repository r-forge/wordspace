##
##  Small co-occurrence tables that may be included as examples in wordspace package
##  .rda files are written to local data/ directory (selected data sets can then be copied to package)
##

library(wordspace)

desc.vn <- read.delim("orig/desc_cc_vn_triples.txt.gz", header=FALSE, quote="", col.names=c("f", "noun", "rel", "verb"), colClasses=c("numeric", "character", "character", "character"), fileEncoding="utf8")
DSM_VerbNounTriples_DESC <- desc.vn[, c("noun", "rel", "verb", "f")]

str(DSM_VerbNounTriples_DESC)
save(DSM_VerbNounTriples_DESC, file="data/DSM_VerbNounTriples_DESC.rda", compress="xz", compression_level=6)


bnc.vn <- read.delim("orig/bnc_cc_vn_triples.txt.gz", header=FALSE, quote="", col.names=c("f", "noun", "rel", "verb", "mode"), colClasses=c("numeric", "character", "character", "character", "factor"), fileEncoding="utf8")
DSM_VerbNounTriples_BNC <- bnc.vn[, c("noun", "rel", "verb", "f", "mode")]

str(DSM_VerbNounTriples_BNC)
save(DSM_VerbNounTriples_BNC, file="data/DSM_VerbNounTriples_BNC.rda", compress="xz", compression_level=6)
