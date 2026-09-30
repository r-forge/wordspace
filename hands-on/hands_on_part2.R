######################################################################
## Code examples from Part 2 of the tutorial (Taxonomy of DSM Parameters),
## including some hands-on exercises
##

## Try these code examples as you follow the lecture and modify them!

library(wordspace)

TC <- DSM_TermContext
head(TC, Inf)

TT <- DSM_TermTerm
head(TT, Inf)

## TT is a 'dsm' object with marginal frequencies
TT$rows
TT$cols
TT$globals$N
TT$M # the actual co-oc matrix

options(digits=3)
dsm.score(TT, score="MI", sparse=FALSE, matrix.only=TRUE)

dsm.score(TT, score="MI", matrix.only=TRUE)
dsm.score(TT, score="simple-ll", transform="log", matrix.only=TRUE)
?dsm.score # try other measures & transformations listed in the help page

## in wordspace v0.2-2 and newer you can define your own association measures
logDice <- function (O, R1, C1, ...) 14 + log2(2 * O / (R1 + C1)) # used by SketchEngine
dsm.score(TT, score=logDice, sparse=FALSE, matrix.only=TRUE)

TT <- dsm.score(TT, score="frequency", transform="log")
TT$S # co-oc matrix after feature scaling

pair.distances(c("cat","cause"), c("animal","effect"), TT, method="euclidean")

dist.matrix(TT, method="euclidean")
dist.matrix(TT, method="minkowski", p=4)

rowNorms(TT$S, method="euclidean")
TT <- dsm.score(TT, score="freq", transform="log", normalize=TRUE, method="euclidean")
rowNorms(TT$S, method="euclidean")
dist.matrix(TT, method="euclidean")

dist.matrix(TT, method="cosine", convert=FALSE) # cosine similarity
dist.matrix(TT, method="cosine") # angular distance metric

?dist.matrix # try other distance / similarity measures listed in the help page
## wordspace v0.2-2 offers Jaccard and overlap for non-negative DSM vectors


######################################################################
## Projection methods for dimensionality reduction
##

TT2 <- dsm.projection(TT, n=2, method="svd")
TT2

x <- TT2[, 1]
y <- TT2[, 2]
plot(x, y, pch=20, col="red", xlim=extendrange(x), ylim=extendrange(y))
text(x, y, rownames(TT2), pos=3)

## wrap plot in ad-hoc function so we can apply it repeatedly (with a few more options)
latent.plot <- function (M, dim=1:2, xlim=extendrange(x), ylim=extendrange(y)) {
  x <- M[, dim[1]]
  y <- M[, dim[2]]
  plot(x, y, pch=20, col="red", xlim=xlim, ylim=ylim, xlab=colnames(M)[dim[1]], ylab=colnames(M)[dim[2]])
  text(x, y, rownames(M), pos=3)
  t(M)[dim, ]
}

## Random Indexing produces different result on each run
latent.plot( dsm.projection(TT, n=2, method="ri") )

## visualise effects of SVD power scaling 
latent.plot(dsm.projection(TT, n=2, method="svd"), xlim=c(-1, 1), ylim=c(-1, 1)) # first SVD dim just points to centroid of data set
latent.plot(dsm.projection(TT, n=3, method="svd"), dim=c(3, 2), xlim=c(-1, 1), ylim=c(-1, 1))
latent.plot(dsm.projection(TT, n=3, method="svd", power=0.5), dim=c(3, 2), xlim=c(-1, 1), ylim=c(-1, 1))
latent.plot(dsm.projection(TT, n=3, method="svd", power=0), dim=c(3, 2), xlim=c(-1, 1), ylim=c(-1, 1))

## bonus for R hackers: post-hoc power scaling (more efficient for experimenting with different values of P)
TT2 <- dsm.projection(TT, n=2, method="svd", power=0) # start from P=0
sigma <- attr(TT2, "sigma") # the singular values
scaleMargins(TT2, cols=sigma ^ 0.5) # compare to P=0.5 below
dsm.projection(TT, n=2, method="svd", power=0.5)


######################################################################
## Building your own DSM from a table of dependency triples
##

subset(DSM_VerbNounTriples_BNC, noun == "dog" & verb == "walk")

tri <- subset(DSM_VerbNounTriples_BNC, rel == "obj")
VObj <- dsm(target=tri$noun, feature=tri$verb, score=tri$f, raw.freq=TRUE)
VObj
head(VObj$rows, 20)

VObj <- dsm.score(VObj, score="MI", normalize=TRUE)
nearest.neighbours(VObj, "dog")
nearest.neighbours(VObj, "dog", method="euclidean") # equivalent to cosine

nearest.neighbours(VObj, "dog", method="manhattan") # NB: incompatible normalization!

## select appropriate norm for row vector normalisation
VObj.L1 <- dsm.score(VObj, score="MI", normalize=TRUE, method="manhattan")
nearest.neighbours(VObj.L1, "dog", method="manhattan")

VObj50 <- dsm.projection(VObj, n=50, method="svd")
nearest.neighbours(VObj50, "dog", method="euclidean") 

## should re-normalise latent vectors for Euclidean distance
VObj50 <- normalize.rows(VObj50, method="euclidean") # applies to matrix, not DSM object
nearest.neighbours(VObj50, "dog", method="euclidean") 

## Exercise:
##  - How many other co-occurrence matrices can you build from DSM_VerbNounTriples_BNC?
##  - Explore different parameter settings for these models.
##  - Compare nearest neighbours of selected words for different co-occurrence matrices and parameters.


######################################################################
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
## All models are available as raw co-occurrence counts with marginal frequencies,
## as well as in a pre-compiled version (log simple-ll, SVD to 500 dimensions, P = 0).

## In this exercise, we will use the raw co-occurrence matrix so that we can test
## different settings for all the DSM parameters.  Download one of the models, e.g.
## 
##    WP500_Win5_Lemma.rda
##
## and save the file to your models/ subdirectory (usually the same as this R script).
## If your computer has limited RAM, it may be better to pick one of the more specific
## co-occurrence contexts resulting in a smaller model (file size < 60 MB).

## Load the co-occurrence matrix and marginal frequency data
load("models/WP500_Win5_Lemma.rda", verbose=TRUE) # verbose option prints name of the DSM object

## It is convenient to assign the model to a shorter variable name
WP <- WP500_Win5_Lemma # R only makes a copy if you modify one of the two variables

## Take a first look at the data
WP # shows number of rows and columns as well as fill rate (-> sparseness)
head(WP, 15) # top left corner of co-occurrence matrix
head(WP$rows, 10) # row marginals (NB: f == 10 * frequency because of span-size adjustment)
head(WP$cols, 10) # column marginals (here f is the corpus frequency of the term)
## Terms in this model are POS-disambiguated lemmas with part-of-speech codes
## _N (noun), _V (verb), _J (adjective), _R (adverb)

## Now apply different scores, transformation normalizations, etc. to the model and
## take a look at nearest neighbours for selected words.  It would be best to always
## look at the same targets as you change parameters, so let's collect them in a vector
## (you might want to select a different set of targets, of course):
words <- c("dog_N", "book_N", "walk_V", "smile_V")

## E.g. PPMI scores with L3 normalization (for Minkowski p=3 distance):
WP <- dsm.score(WP, score="MI", normalize=TRUE, method="minkowski", p=3)

## nearest.neighbours() accepts multiple target words
nearest.neighbours(WP, words, n=10, method="minkowski", p=3)
## NB: method= and p= arguments are documented in ?dist.matrix

## You can also create semantic maps such as the one shown on the first day.
## A simplified version can be obtained with just a few wordspace commands!
nouns <- ESSLLI08_Nouns$word # list of nouns to be visualised in the map
nouns.cat <- ESSLLI08_Nouns$class # corresponding semantic classes for colour-coding

DM <- dist.matrix(WP, terms=nouns, method="manhattan") # distance matrix for these nouns 
plot(DM, col=nouns.cat, show.edges=FALSE, method="sammon") # plot as semantic map
## now repeat these two commands for different DSM parameter settings and distance measures

## Nearest-neighbour (NN) search among the rows of a sparse matrix is fairly slow for 
## technical reasons.  If you run out of patience, here are two things you can do:

## 1) Compute NN for a large number of words you're interested in at the same time:
NN.list <- nearest.neighbours(WP, words, n=20, method="minkowski", p=3)
NN.list$dog_N        # now show individual neighbour sets from the list
NN.list[["smile_V"]] # safer in case a word contains non-standard characters

## 2) If your vesion of R supports multi-threading, try to activate it
wordspace.openmp()
wordspace.openmp(threads=4) # rule of thumb: number of physical CPU cores

nearest.neighbours(WP, words, n=10, method="minkowski", p=3)


## If you also want to experiment with dimensionality reduction, you may need to
## simplify your model by filtering rows and columns.  This will also speed up NN
## search in the unreduced matrix, of course.

## The subset() function allows you to select rows and columns with the arguments
## subset= (rows) and select= (columns).  The filtering expressions are based on
## the variables term, f, nnzero from the row/column marginal tables.
## See the method documentation ?subset.dsm for further information and examples.

## For example, let us only look at verbs as targets and medium-frequency features
## with 400 <= f <= 10000 (based on the histogram below)
hist(log10(WP$cols$f), breaks=40) # 2 = 100, 3 = 1000, 4 = 10000, ...
abline(v=log10(c(400, 10000)), col="red")

WP2 <- subset(WP, subset=(grepl("_V$", term)), select=(f >= 400 & f <= 10000), update.nnzero=TRUE)
WP2 # considerably smaller now

WP2 <- dsm.score(WP2, score="MI", normalize=TRUE, method="minkowski", p=3) # need to renormalize rows (why?)

## Use skip.missing=TRUE to ignore target words not found in the model
nearest.neighbours(WP2, words, n=10, method="minkowski", p=3, skip.missing=TRUE)

## SVD dimensionality reduction usually applied after Euclidean normalization (=> Part 5)
WP2 <- dsm.score(WP2, score="MI", normalize=TRUE) # default normalization is Euclidean

## SVD-based dimensionality reduction to n=100 latent dimension with whitening (P = 0)
WP100 <- dsm.projection(WP2, n=100, method="svd", power=0)
## NB: WP100 is a plain matrix (rows labelled with target terms), not a DSM object
## NB: 100 latent dimensions are often not enough for good semantics => try n=300 or n=500

nearest.neighbours(WP100, words, skip.missing=TRUE) # NN are much faster in reduced matrix

## If you want to experiment with power scaling, it's much faster to adjust post-hoc:
sigma <- attr(WP100, "sigma")
WP100.P1 <- scaleMargins(WP100, cols=(sigma ^ 1)) # because we started with P = 0

## You can use matrix indexing to skip latent dimensions
WP100.P1.skip <- WP100.P1[, 21:100] # skips first 20 latent dims

## Don't forget to re-normalize the vectors appropriately (unless you're using cosine/angle)
WP100.P1.skip <- normalize.rows(WP100.P1.skip, method="manhattan")
nearest.neighbours(WP100.P1.skip, words, n=10, method="manhattan", skip.missing=TRUE)

## Exercise:
##  - experiment with all the parameters you have learned in the tutorial
##  - can you clearly identify "good" and "poor" parameter settings?
##  - are there parameter settings you cannot realize with the wordspace package?
