##
## Exercise for part 4:
## Roll your own DSM
##

library(wordspace)

## Load a relatively dense co-occurrence matrix from the WP500 corpus, e.g. L5/R5
load("WP500_Win5_Lemma.rda", verbose=TRUE)

## Your implementation will use standard matrix algebra for dense matrices built into R.
## Therefore you should first reduce the matrix to the most frequent targets and features.
## Keep in mind that a dense representation of the full matrix would require
cat(sprintf("%.1f GiB RAM\n", prod(dim(WP500_Win5_Lemma)) * 8 / 2^30))

## Real-life DSMs need to be processed with sparse matrix algebra (R package: "Matrix"),
## but this can be quite tricky. The wordspace package is carefully designed to carry out
## most operations efficiently with sparse matrices.

## Here's how you select the top-1000 targets and features.  This DSM will describe nouns
## as targets by adjectives as features, so we first apply a POS filter:
Model <- subset(WP500_Win5_Lemma, grepl("_N$", term), grepl("_J$", term))

## Now we can select the most frequent items from the remaining targets and features
## (rank(-f) computes ranks for the marginal frequencies in decreasing order):
Model <- subset(Model, rank(-f) <= 1000, rank(-f) <= 1000)

## Finally, extract the co-occurrence matrix, row/column marginals and sample size:
M <- as.matrix(Model$M) # convert to dense matrix
r <- Model$rows$f
c <- Model$cols$f
N <- Model$globals$N


## Your task:
##  - apply the steps from the lecture slide to build your own DSM from these data
##  - can you compute other association scores, transformations, ... ?
##  - how long does the dense SVD projection take? compare it to dsm.projection(Model)


