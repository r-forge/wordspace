######################################################################
## Code examples and hands-on exercise for DSM tutorial Part 3
##

library(wordspace)
library(wordspaceEval)

M <- DSM_Vectors
nearest.neighbours(M, "walk_V")


######################################################################
## TASK 1: MULTIPLE CHOICE

# DATASET: TOEFL
head(TOEFL80)
eval.multiple.choice(TOEFL80, M) 

# With details=TRUE we can look at how well we did on the individual items. 
eval.multiple.choice(TOEFL80, M, details=TRUE) 
eval.multiple.choice(TOEFL80, M)
eval.multiple.choice(TOEFL80, M, method="manhattan") # manhattan is slightly better
eval.multiple.choice(TOEFL80, M, method="maximum")   # maximum a disaster :)

# Comparing distance to rank
eval.multiple.choice(TOEFL80, M) 
eval.multiple.choice(TOEFL80, M, rank="fwd") # unsurprisingly identical to distance, given that the ranking is the same.  
eval.multiple.choice(TOEFL80, M, rank="bwd") # it is bwd rank which can make a difference here - and it does, but not for the better

# An interesting field of the output of details=TRUE is "correct.rank" 
# which tells us, among the 4 choices,
# where did the correct one land (1 means the item was predicted correctly).
# In a way, this is a more "forgiving" way of assessing performance,
# one where getting the right choice in position 2 is not as wrong 
# as getting it as the last one.

# cosine
mean(eval.multiple.choice(TOEFL80, M, details=TRUE)$correct.rank)
# manhattan
mean(eval.multiple.choice(TOEFL80, M, method="manhattan", details=TRUE)$correct.rank)
# maximum
mean(eval.multiple.choice(TOEFL80, M, method="maximum", details=TRUE)$correct.rank)


######################################################################
## TASK 2: PREDICTION OF SIMILARITY RATINGS

# DATASET: RG65
head(RG65)
# Let us just have a look at a handful of items (from 1 to 61 in steps of 5)
RG65[seq(1,61,5), ]

eval.similarity.correlation(RG65, M, convert=FALSE)
# With details=TRUE we can look at the individual items
eval.similarity.correlation(RG65, M, details=TRUE, convert=FALSE)
# ... and plot the correlation between similarity/distance and rating
plot(eval.similarity.correlation(RG65, M, details=TRUE, convert=FALSE))

# DATASET: WordSim353
head(WordSim353)

eval.similarity.correlation(WordSim353, M, convert=FALSE)
plot(eval.similarity.correlation(WordSim353, M, details=TRUE, convert=FALSE))

# WordSim353 contains both similarity and relatedness ratings, indicated by Boolean variables
# relatedness and similarity. We can use this information to conduct separate evaluations.
# E.g, looking at items annotated for relatedness
eval.similarity.correlation(subset(WordSim353, relatedness), M, convert=FALSE)
# versus items annotated for similarity
eval.similarity.correlation(subset(WordSim353, similarity), M, convert=FALSE)
# Yeah, our model likes similarity more, clearly... 

# BONUS: Note that you can use the output from details=TRUE for further analysis, e.g. with a regression model
result <- eval.similarity.correlation(WordSim353, M, convert=FALSE, details=TRUE)
summary(lm(distance ~ score + relatedness, data=result))
# the linear model confirms us that the highest scores have higher similarties
# (even if it is called distance, remember we had convert=FALSE)
# and that when the items are annotated as TRUE for relatedness, 
# their similarity score is lower 


######################################################################
# TASK 3: CONCEPT CLUSTERING

# DATASET: ESSLLI08_Nouns
head(ESSLLI08_Nouns)
ESSLLI08_Nouns[seq(1,40,5), ]

eval.clustering(ESSLLI08_Nouns, M)
eval.clustering(ESSLLI08_Nouns, M, details=TRUE)

# Note that the ESSLLI dataset has multiple annotations available, at different granularities
summary(ESSLLI08_Nouns)

# The eval.clustering function allows you to specify which column of your dataset 
# is to be used as a gold standard class for clustering
# 6 classes
eval.clustering(ESSLLI08_Nouns, M)
# 3 classes
eval.clustering(ESSLLI08_Nouns, M, class.name="class2")
# 2 classes
eval.clustering(ESSLLI08_Nouns, M, class.name="class3")
# clearly the task is easier with fewer classes to be predicted...


######################################################################
# DIGGING DOWN: ERROR ANALYSIS

# **TODO** example code -> rework rough draft from workshop below
# - see if we can do it with DSM_Vectors, but probably need to load WP500_Win5 as below
# - always runs standard eval first, then with details=TRUE for error analysis
# - add step-by-step explanations and interpretations
# - work out error analysis for similarity correlation

res <- eval.multiple.choice(TOEFL80, M, details=TRUE)
head(res)
summary(res$correct)
err <- subset(res, !correct)
err
nearest.neighbours(M, "perseverance_N")
nearest.neighbours(M, "endurance_N")
nearest.neighbours(M, "advent_N")

res <- eval.clustering(ESSLLI08_Nouns, M, details=TRUE)
head(res)
xtabs(~ paste(cluster, label) + gold, data=res)
subset(res, cluster == 3)
nearest.neighbours(M, "snail_N", 16)

res <- eval.similarity.correlation(RG65, M, details=TRUE)
res[seq(1, 65, 8), ]
# need regression to work out what biggest errors are
# try via smoother with loess() + predict()


######################################################################
## Exercise for Part 3:
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
load("models/WP500_Win5_Lemma_svd500.rda", verbose=TRUE) # verbose option prints name of the DSM matrix

M <- WP500_Win5_Lemma_svd500 # now assign the matrix to a shorter variable name

## compute neighbours for selected terms (lemma_POS format)
nearest.neighbours(M, "love_V")
nearest.neighbours(M, "love_V", method="manhattan")

plot(nearest.neighbours(M, "semantics_N", n=20, dist.matrix=TRUE))

## compute semantic map for selected words (make your own set!)
words <- unique(c(RG65$word1, RG65$word2))
words
plot(dist.matrix(M, terms=words, skip.missing=TRUE), show.edges=FALSE) # ignore missing words

## evaluate model in various tasks
eval.clustering(ESSLLI08_Nouns, M) # also try various distance measures
plot(eval.similarity.correlation(RG65, M, details=TRUE))
eval.similarity.correlation(WordSim353, M) # larger & more difficult data set
eval.multiple.choice(TOEFL80, M)

## Check out all data sets in the wordspaceEval package and read their help pages:
data(package="wordspaceEval")
?AP402
?Battig82
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
