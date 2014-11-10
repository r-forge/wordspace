##
##  Pre-compiled reduced DSM vectors for small selected vocabulary.
##

library(wordspace)
library(wordspaceEval)

argv <- commandArgs(trailingOnly=TRUE)
if (length(argv) != 1) stop("Usage:  R --no-save -f precompiled_dsm.R --args <path to Web-based DSMs>")

large.dsm <- argv[1] # path where large Web-based DSM can be found
load(paste(large.dsm, "rda/web_150k_30k_l4r4_svd1000.rda", sep="/"), verbose=TRUE)

## vocabulary should cover basic evaluation tasks
vocab <- ESSLLI08_Nouns$word
vocab <- union(vocab, with(RG65, c(word1, word2)))
vocab <- union(vocab, with(WordSim353, c(word1, word2)))
vocab <- union(vocab, with(TOEFL80, c(target, correct, distractor1, distractor2, distractor3)))
vocab <- union(vocab, AP402$word)

## as well as terms from example matrices
vocab <- union(vocab, paste0(rownames(DSM_HieroglyphsMatrix), "_N"))
vocab <- union(vocab, paste0(colnames(DSM_HieroglyphsMatrix), "_V"))
vocab <- union(vocab, paste0(rownames(DSM_TermContextMatrix), "_N"))
# - we omit columns of DSM_TermTermMatrix because they are from different parts of speech

## as well as nearest neighbours of a few selected words
words <- c("white_J",                   # includes many colours
           "apple_N",                   # various fruit
           "kindness_N",                # always good :-)
           "walk_V") 
vocab <- union(vocab, words)
for (w in words) vocab <- union(vocab, names(nearest.neighbours(web_150k_30k_l4r4_svd1000, w, 40)))


vocab <- vocab[vocab %in% rownames(web_150k_30k_l4r4_svd1000)]
length(vocab) # 1323 target words remaining

DSM_Vectors <- web_150k_30k_l4r4_svd1000[vocab, 1:100]
DSM_Vectors <- normalize.rows(DSM_Vectors) # renormalize row vectors
dim(DSM_Vectors)
cat(sprintf("Approx. size: %.2f MiB\n", object.size(DSM_Vectors) / 2^20))

save(DSM_Vectors, file="data/DSM_Vectors.rda", compress="xz", compression_level=9)

