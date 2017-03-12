## Exercise for part 3:
## Explore and evaluate distributional semantic models

library(wordspace)
library(wordspaceEval) # if you have the non-public data sets

## Several pre-compiled DSMs based on the English Wikipedia (WP500 corpus)
## using different co-occurrence contexts are available for download from 
##
##     http://www.collocations.de/data/#dsm
##
## The following contexts are available:
##   TermDoc  ...   term-document matrix
##   Win30    ...   30-word span (L30/R30)
##   Win5     ...   5-word span (L5/R5)
##   Win2     ...   2-word span (L2/R2)
##   DepFilter ...  dependency-filtered
##   DepStruct ...  dependency-structured
##
## All models are available as raw co-occurrence counts and in a pre-compiled
## version with SVD dimensionality reduction to 500 dimensions.

## Start by loading one of the pre-compiled models ("*_svd500.rda").
## It is easiest to put the downloaded file in the same directory as this R script.
load("WP500_Win5_Lemma_svd500.rda", verbose=TRUE) # verbose option prints name of the DSM matrix

M <- WP500_Win5_Lemma_svd500 # now assign the matrix to a shorter variable name

## compute neighbours for selected terms (lemma_POS format)
nearest.neighbours(M, "love_V")
nearest.neighbours(M, "love_V", method="maximum")

plot(nearest.neighbours(M, "semantics_N", n=20, dist.matrix=TRUE))

## compute semantic map for selected words (make your own set!)
words <- unique(c(RG65$word1, RG65$word2))
words
plot(dist.matrix(M, terms=words, skip.missing=TRUE)) # ignore missing words

## evaluate model in various tasks
eval.clustering(ESSLLI08_Nouns, M) # also try various distance measures
plot(eval.similarity.correlation(RG65, M, details=TRUE))
eval.similarity.correlation(WordSim353, M) # larger & more difficult data set
eval.multiple.choice(TOEFL80, M)


## Your task:
##  - explore distances, neighbours and semantic maps for different DSMs and distance measures
##  - evaluate each model in the three standard tasks listed above
## Summarize your findings:
##  - How different are the various co-occurrence contexts?
##  - Can you put a finger on the kind of semantic relations in each model?
##  - Are some parts of speech, semantic classes, etc. represented better than others?
##  - How much influence does the distance measure have?  Is one measure better than all others?
