##
##  Small co-occurrence tables that may be included as examples in wordspace package
##  .rda files are written to local data/ directory (selected data sets can then be copied to package)
##

library(wordspace)

desc.vn <- read.delim("orig/desc_cc_vn_triples.txt.gz", header=FALSE, quote="", col.names=c("f", "noun", "rel", "verb"), colClasses=c("numeric", "character", "character", "character"), fileEncoding="utf8")
DSM_VerbNounTriples_DESC <- subset(desc.vn, f >= 2, c("noun", "rel", "verb", "f"))

str(DSM_VerbNounTriples_DESC)
save(DSM_VerbNounTriples_DESC, file="data/DSM_VerbNounTriples_DESC.rda", compress="xz", compression_level=6)


bnc.vn <- read.delim("orig/bnc_cc_vn_triples.txt.gz", header=FALSE, quote="", col.names=c("f", "noun", "rel", "verb", "mode"), colClasses=c("numeric", "character", "character", "character", "factor"), fileEncoding="utf8")
DSM_VerbNounTriples_BNC <- subset(bnc.vn, f >= 5, c("noun", "rel", "verb", "f", "mode")) # use f >= 5 for reduced data size

str(DSM_VerbNounTriples_BNC)
save(DSM_VerbNounTriples_BNC, file="data/DSM_VerbNounTriples_BNC.rda", compress="xz", compression_level=6)

dim(DSM_VerbNounTriples_BNC); object.size(DSM_VerbNounTriples_BNC) / 1e6
TRUE

if (TRUE) {
  ## validate new data extraction (without duplicate entries) against old data sets
  ## NB: new tuples are sorted lexicographically, so need to resort DSM matrices for comparison
 
  bnc.vn.old <- read.delim("orig/bnc_cc_vn_triples_old.txt.gz", header=FALSE, quote="", col.names=c("f", "noun", "rel", "verb", "mode"), colClasses=c("numeric", "character", "character", "character", "factor"), fileEncoding="utf8")
  Mold <- with(bnc.vn.old, dsm(target=noun, feature=paste(rel, verb), score=f, raw.freq=TRUE, sort=TRUE, verbose=TRUE))
  Mnew <- with(bnc.vn, dsm(target=noun, feature=paste(rel, verb), score=f, raw.freq=TRUE, sort=TRUE, verbose=TRUE))
  all.equal(dim(Mold), dim(Mnew))
  all.equal(Mold$M, Mnew$M)

  desc.vn.old <- read.delim("orig/desc_cc_vn_triples_old.txt.gz", header=FALSE, quote="", col.names=c("f", "noun", "rel", "verb"), colClasses=c("numeric", "character", "character", "character"), fileEncoding="utf8")
  Mold <- with(desc.vn.old, dsm(target=noun, feature=paste(rel, verb), score=f, raw.freq=TRUE, sort=TRUE, verbose=TRUE))
  Mnew <- with(desc.vn, dsm(target=noun, feature=paste(rel, verb), score=f, raw.freq=TRUE, sort=TRUE, verbose=TRUE))
  all.equal(dim(Mold), dim(Mnew))
  all.equal(Mold$M, Mnew$M)
}
