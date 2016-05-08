##
##  Artificial example matrices (defined inline or derived from co-occurrence tables)
##  .rda files are directly written to package data/ directory
##  prerequisites: data/DSM_VerbNounTriples_BNC.rda
##

library(wordspace)

##
## Hieroglyphs example (with English labels)
##
hieroglyphs.txt <- "get see use hear eat kill
  knife   51  20  84    0   3    0
  cat     52  58   4    4   6   26
  dog    115  83  10   42  33   17
  boat    59  39  23    4   0    0
  cup     98  14   6    2   1    0
  pig     12  17   3    2   9   27
  banana  11   2   2    0  18    0
"
fh <- textConnection(hieroglyphs.txt)
DSM_HieroglyphsMatrix <- as.matrix(read.table(fh))
mode(DSM_HieroglyphsMatrix) <- "double"
close(fh)

print(DSM_HieroglyphsMatrix)
str(DSM_HieroglyphsMatrix)
save(DSM_HieroglyphsMatrix, file="../pkg/wordspace/data/DSM_HieroglyphsMatrix.rda", compress="gzip")


##
## Wikipedia term-context matrix
##
tc.txt <- "Felidae Pet Feral Boat Philosophy Kant Back_Pain
cat             10  10     7    0          0    0         0
dog              0  10     4   11          0    0         0
animal           2  15    10    0          0    0         0
time             1   0     0    2          2    1         0
reason           0   1     0    0          1    4         1
cause            0   0     0    2          1    2         6
effect           0   0     0    1          0    1         0
"
fh <- textConnection(tc.txt)
DSM_TermContextMatrix <- Matrix(as.matrix(read.table(fh))) # converts to mode double
close(fh)

row_info <- data.frame(
  term = rownames(DSM_TermContextMatrix),
  f = c(6103, 14180, 21538, 316591, 26470, 15913, 36933),
  stringsAsFactors = FALSE)
col_info <- data.frame(
  term = colnames(DSM_TermContextMatrix),
  f = c(270, 305, 308, 313, 273, 263, 310),
  stringsAsFactors = FALSE)

## corresponding DSM object with sparse representation
DSM_TermContext <- dsm(M=Matrix(DSM_TermContextMatrix), rowinfo=row_info, colinfo=col_info, N=108771103, raw.freq=TRUE)

if (FALSE) {
  ## extract from original WP500 data set with marginal frequencies
  load("~/Project/DSM/models/WP500/rda/WP500_TermDoc_Lemma.rda")
  model <- WP500_TermDoc_Lemma
  targets <- c("cat_N", "dog_N", "animal_N", "time_N", "reason_N", "cause_N", "effect_N")
  features <- c("Felidae", "Pet", "Feral", "Boat", "Philosophy", "Immanuel Kant", "Back pain")
  idxR <- match(targets, model$rows$term)
  idxC <- match(features, model$cols$term)
  tc <- subset(model, idxR, idxC)
  head(tc, Inf) # co-occurrence matrix (slightly modified in example for illustrational purposes)
  tc$rows       # row and column marginals
  tc$cols
  tc$globals$N  # underlying sample size
}

print(DSM_TermContextMatrix)
str(DSM_TermContextMatrix)
print(DSM_TermContext)
head(DSM_TermContext, Inf)

save(DSM_TermContextMatrix, file="../pkg/wordspace/data/DSM_TermContextMatrix.rda", compress="gzip")
save(DSM_TermContext, file="../pkg/wordspace/data/DSM_TermContext.rda", compress="gzip")


##
## Export term-context matrix as triplet file and in UCS format
##
MT <- as(DSM_TermContextMatrix, "dgTMatrix")
triplets <- data.frame(
  target = rownames(MT)[MT@i + 1],
  feature = colnames(MT)[MT@j + 1],
  f = MT@x)
write.table(triplets, file=gzfile("../pkg/wordspace/inst/extdata/term_context_triplets.gz"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

ucs.ds <- with(
  DSM_TermContext,
  data.frame(
    l1 = rows$term[MT@i + 1],
    l2 = cols$term[MT@j + 1],
    f = MT@x,
    f1 = rows$f[MT@i + 1],
    f2 = cols$f[MT@j + 1],
    N = globals$N))
write.table(ucs.ds, file=gzfile("data/term_context.ds.gz"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

## Now export to read.ucs() format with UCS command-line tools:
##   ucs-tool export-dsm-matrix -S -gz -ts alpha -fs alpha -f data/term_context.ds.gz ../pkg/wordspace/inst/extdata/term_context_ucs

##
## Wikipedia term-term matrix
##
tt.txt <- "breed    tail   feed   kill   important   explain   likely
cat           84      17      8     38           0         2        0
dog          579      14     32     63           1         2        2
animal        45      11     86    136          13         5        4
time          19       8     29    134          94        44      100
reason         1       0      1     18          71       140       39
cause          0       1      0      3          55        35       51
effect         0       1      1      6          62        37       14
"
## -- previous version of the matrix
## cat           83    17     7    37          0        1       0
## dog          561    13    30    60          1        2       4
## animal        42    10   109   134         13        5       5
## time          19     9    29   117         81       34     109
## reason         1     0     2    14         68      140      47
## cause          0     1     0     4         55       34      55
## effect         0     0     1     6         60       35      17
fh <- textConnection(tt.txt)
DSM_TermTermMatrix <- as.matrix(read.table(fh))
mode(DSM_TermTermMatrix) <- "double"
close(fh)

row_info <- data.frame(
  term = rownames(DSM_TermTermMatrix),
  f = c(22007, 50807, 77053, 1156693, 95047, 54739, 133102),
  stringsAsFactors = FALSE)
col_info <- data.frame(
  term = colnames(DSM_TermTermMatrix),
  f = c(10329, 6511, 8307, 49534, 53824, 17939, 11096),
  stringsAsFactors = FALSE)

## corresponding DSM object with sparse representation
DSM_TermTerm <- dsm(M=Matrix(DSM_TermTermMatrix), rowinfo=row_info, colinfo=col_info, N=199902178, raw.freq=TRUE)

if (FALSE) {
  ## extract from original WP500 data set with marginal frequencies
  load("~/Project/DSM/models/WP500/rda/WP500_Win2_Lemma.rda")
  model <- WP500_Win2_Lemma 
  targets <- c("cat_N", "dog_N", "animal_N", "time_N", "reason_N", "cause_N", "effect_N")
  features <- c("breed_N", "breed_V", "tail_N", "feed_V", "kill_V", "important_J", "explain_V", "likely_J")
  idxR <- match(targets, model$rows$term)
  idxC <- match(features, model$cols$term)
  tc <- subset(model, idxR, idxC)
  head(tc, Inf) # co-occurrence matrix (partially adapted, original data set w/o POS disambiguation)
  tc$rows       # row and column marginals
  tc$cols       # only major difference from headword model is breed_N + breed_V, which should be aggregated
  tc$globals$N  # underlying sample size
}

print(DSM_TermTermMatrix)
str(DSM_TermTermMatrix)
print(DSM_TermTerm)
head(DSM_TermTerm, Inf)

save(DSM_TermTermMatrix, file="../pkg/wordspace/data/DSM_TermTermMatrix.rda", compress="gzip")
save(DSM_TermTerm, file="../pkg/wordspace/data/DSM_TermTerm.rda", compress="gzip")


##
## Illustration of dimensionality reduction: nouns denoting goods dimensions "buy", "sell" and "own"
##
load("data/DSM_VerbNounTriples_BNC.rda")
VObj <- with(subset(DSM_VerbNounTriples_BNC, rel=="obj"), dsm(target=noun, feature=verb, score=f, raw.freq=TRUE))

Goods <- subset(VObj, select=(term %in% c("buy", "sell", "own"))) # updates nnzero counts
Goods <- dsm.score(Goods, score="simple-ll", transform="log", sparse=FALSE, negative.ok=TRUE, update.nnzero=FALSE) # log-transformed NON-sparse G2 scores; make sure nnzero counts still refer to cooc freqs rather than G2 scores
Goods <- subset(Goods, nnzero >= 3) # only include nouns that occur with all 3 verbs (implies row marginal >= 6 for this data set)

res <- prcomp(Goods$S) # determine centrality of data poins by whitening coordinates
S.white <- scaleMargins(res$x, cols=1 / res$sdev)
fringeness <- rowNorms(S.white, method="manhattan") # good indicator whether data point is near the fringe of the distribution
fringeness <- rank(fringeness) / length(fringeness) # transform to quantiles in range 0 .. 1

if (FALSE) {
  ## visualize data set during development
  plot3d(Goods$S, col=1+(fringeness >= .8))
  idx <- fringeness >= .8
  text3d(Goods$S[idx,], texts=Goods$rows$term[idx], adj=c(0.5,1))
}

DSM_GoodsMatrix <- cbind(Goods$S, fringe=fringeness) # first 3 columns are coordinates, 4th column shows fringeness
save(DSM_GoodsMatrix, file="../pkg/wordspace/data/DSM_GoodsMatrix.rda", compress="gzip")


