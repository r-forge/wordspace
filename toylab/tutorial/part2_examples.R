##
## Code examples from tutorial part 2
##

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

TT2 <- dsm.projection(TT, n=2, method="svd")
TT2

x <- TT2[, 1]
y <- TT2[, 2]
plot(x, y, pch=20, col="red", xlim=extendrange(x), ylim=extendrange(y))
text(x, y, rownames(TT2), pos=3)

## wrap plot in ad-hoc function so we can apply it repeatedly
latent.plot <- function (M) {
  x <- M[, 1]
  y <- M[, 2]
  plot(x, y, pch=20, col="red", xlim=extendrange(x), ylim=extendrange(y))
  text(x, y, rownames(M), pos=3)
}

latent.plot( dsm.projection(TT, n=2, method="ri") ) # RI produces different result on each run

dsm.projection(TT, n=2, method="svd") # effects of power scaling
dsm.projection(TT, n=2, method="svd", power=0.5)
dsm.projection(TT, n=2, method="svd", power=0)

TT2 <- dsm.projection(TT, n=2, method="svd", power=0) # post-hoc power scaling: start from P=0
sigma <- attr(TT2, "sigma") # the singular values
scaleMargins(TT2, cols=sigma ^ 0.5) # compare to P=0.5 above


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
