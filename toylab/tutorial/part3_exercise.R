##
## Exercise for part 3:
## Explore and evaluate distributional semantic models
##


library(wordspace)
## evaluation tasks:
##   ESSLLI08_Nouns ... clustering
##   RG65           ... similarity ratings
##   WordSim353     ... similarity ratings
##   SemCorWSD      ... word sense disambiguation (Schuetze-style)

library(wordspaceEval) # if you have the non-public data sets
## additional non-public evaluation tasks:
##   TOEFL80        ... multiple choice (synonyms)
##   SPP_Items      ... multiple choice (various relations)
##   GEK_Items      ... multiple choice (various relations)
##   AP402          ... clustering
##   Battig82       ... clustering


## Several pre-compiled DSMs based on the English Wikipedia (WP500 corpus)
## using different co-occurrence contexts are available for download from 
##
##    http://wordspace.collocations.de/doku.php/course:material#pre-compiled_dsms
##
## The following co-occurrence contexts are available:
##   TermDoc   ...  term-document matrix
##   Win30     ...  30-word span (L30/R30)
##   Win5      ...  5-word span (L5/R5)
##   Win2      ...  2-word span (L2/R2)
##   DepFilter ...  dependency-filtered
##   DepStruct ...  dependency-structured
##   Ctype_L1R1 ... L1+R1 context types (pattern of left & right word)
##   Ctype_L2R2 ... L2+R2 context types (left & right 2 words, very sparse)
##   Ctype_L2R2pos ... L2+R2 part-of-speech context types
##
## You can also try two non-lemmatized models, but you will have to specify the
## option format="HWLC" for all evaluation functions (see ?convert.lemma).
##
## All models are available as raw co-occurrence counts with marginal frequencies,
## as well as in a pre-compiled version (log simple-ll, SVD to 500 dimensions, P = 0).


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

## Check out all data sets in the wordspaceEval package and read their help pages:
data(package="wordspaceEval")
?SPP_Items


## Your task:
##  - explore distances, neighbours and semantic maps for different DSMs and distance measures
##  - remember that you can apply post-hoc power scaling and skip dimensions to the SVD-reduced DSM
##  - evaluate each model in the three standard tasks listed above
## Summarize your findings:
##  - How different are the various co-occurrence contexts?
##  - Can you put a finger on the kind of semantic relations in each model?
##  - Are some parts of speech, semantic classes, etc. represented better than others?
##  - How much influence does the distance measure have?  Is one measure better than all others?
