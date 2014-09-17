##
##  Pre-compiled reduced DSM vectors for small selected vocabulary.
##

library(wordspace)
library(wordspaceEval)

large.dsm <- "." # path where large Web-based DSM can be found
load(paste(large.dsm, "rda/web_150k_30k_l4r4_svd1000.rda", sep="/"), verbose=TRUE)

## vocabulary should cover basic evaluation tasks
vocab <- ESSLLI08_Nouns$word
vocab <- union(vocab, with(WordSim353, c(word1, word2)))
vocab <- union(vocab, with(TOEFL80, c(target, correct, distractor1, distractor2, distractor3)))
vocab <- union(vocab, AP402$word)

## as well as nearest neighbours of a few selected words
words <- c("white_J",                   # includes many colours
           "apple_N",                   # various fruit
           "kindness_N",                # always good :-)
           "walk_V") 
vocab <- union(vocab, words)
for (w in words) vocab <- union(vocab, names(nearest.neighbours(web_150k_30k_l4r4_svd1000, w, 40)))


vocab <- vocab[vocab %in% rownames(web_150k_30k_l4r4_svd1000)]
length(vocab) # 1314 target words remaining

DSM_Vectors <- web_150k_30k_l4r4_svd1000[vocab, 1:100]
DSM_Vectors <- normalize.rows(DSM_Vectors) # renormalize row vectors
dim(DSM_Vectors)
cat(sprintf("Approx. size: %.2f MiB\n", object.size(DSM_Vectors) / 2^20))

save(DSM_Vectors, file="data/DSM_Vectors.rda", compress="xz", compression_level=9)

