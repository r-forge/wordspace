## Exercise for part 4:
## Explore and evaluate DSM parameters

library(wordspace)
## evaluation tasks:
##   ESSLLI08_Nouns ... clustering
##   RG65           ... similarity ratings
##   WordSim353     ... similarity ratings
##   SemCorWSD      ... word sense disambiguation (Schuetze-style)

library(wordspaceEval)
## additional non-public evaluation tasks:
##   TOEFL80        ... multiple choice (synonyms)
##   SPP_Items      ... multiple choice (various relations)
##   GEK_Items      ... multiple choice (various relations)
##   AP402          ... clustering
##   Battig82       ... clustering

## Co-occurrence data for several DSMs based on the English Wikipedia (WP500 corpus)
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
##   Ctype_L1R1 ... type context: left + right word (lemma)
##   Ctype_L2R2pos ... type context: 2 left + 2 right POS pattern around target 
##
## You can also try two non-lemmatized models, but you will have to specify the
## option format="HWLC" for all evaluation functions (see ?convert.lemma).

## Download one of the raw co-occurrence data sets (not a pre-compiled DSM) and
## read it into R.  We will use the Win30 context below, but you may want to choose
## a smaller model depending on how powerful your computer is.

## It is easiest to store the model in your working directory (i.e. the same 
## directory as this script), so you don't have to specify the directory path.  
load("WP500_Win30_Lemma.rda", verbose=TRUE)

## Most models have a co-occurrence matrix of approx. 50,000 x 50,000 rows and columns
WP500_Win30_Lemma


## If you have sufficient amounts of RAM (at least 8 GB) and patience (you're willing to wait
## 30 minutes or more for SVD dimensionality reduction), you can work on the full matrix.
## Otherwise it's probably a good idea to reduce the number of rows and columns with the
## subset() function, which has a special method for DSM objects.  You can specify two conditions
## for the rows and columns, respectively (see ?subset.dsm for more options).

## Use grepl() to filter terms with regular expressions, e.g. for adjectives as targets and 
## nouns as features:
DSM <- subset(WP500_Win30_Lemma, grepl("_J$", term), grepl("_N$", term))
DSM # now ca. 6k x 38k 

## Or you can filter out low-frequency target and feature terms.  Let us look at the distribution
## of marginal frequencies first.
hist(log10(WP500_Win30_Lemma$rows$f)) # row marginals are scaled by span size!
hist(log10(WP500_Win30_Lemma$cols$f)) # 2 = 100, 3 = 1000, 4 = 10,000, 5 = 100,000, 6 = 1,000,000

## Let us keep targets with R >= 10,000 (choose a different suitable number for other DSMs!)
## and features in a mid-frequency range 1000 <= f <= 20,000
DSM <- subset(WP500_Win30_Lemma, f >= 10000, f >= 1000 & f <= 20000)
DSM # approx. 30k x 10k now

## If you're short on RAM, delete the original model now and clean up
rm(WP500_Win30_Lemma)
gc() # run garbage collector to free up RAM


## Experiment with an unreduced model first, which is less time-consuming than SVD
DSM <- dsm.score(DSM, score="Dice") # Dice scores without normalization

## Look at some nearest neighbours and evaluate the model in various task. Try different
## distance metrics. Go back and change parameters, then re-run all evaluation step (this
## is particularly convenient with an R script in RStudio).
nearest.neighbours(DSM, "mouse_N")
nearest.neighbours(DSM, "mouse_N", method="euclidean")
nearest.neighbours(DSM, "mouse_N", method="manhattan")

eval.multiple.choice(TOEFL80, DSM) # note that the missing items count as errors!
## Try many other tasks here!

## You can also check directly whether normalization seems to be necessary by plotting
## the distribution of row vector norms (need to specify the matrix rather than DSM object):
hist(rowNorms(DSM$S, method="euclidean"))
hist(rowNorms(DSM$S, method="manhattan")) # now huge differences here


## Once you're satisfied with your DSM, try to see if dimensionality reduction improves
## the representation.  The more latent dimensions you ask for, the longer this will take,
## but you can then select the first r dimensions from the reduced matrix or skip a few.
DSM300 <- dsm.projection(DSM, n=300, method="svd")

nearest.neighbours(DSM300, "mouse_N") # most people use cosine similarity with SVD-reduced models
eval.multiple.choice(TOEFL80, DSM300)

## Use matrix subsetting to pick fewer dimensions or skip the first dimensions
M <- DSM300[, 1:100]  # first 100 dim's only
M <- DSM300[, 51:150] # skip 50, then take next 100 dim's

nearest.neighbours(M, "mouse_N")
eval.multiple.choice(TOEFL80, M)
