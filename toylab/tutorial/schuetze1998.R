##
## Bonus practice session:
## Words Sense Induction as proposed by Schuetze (1998)
##

library(wordspace)
library(cluster) # for better clustering algorithm than kmeans()

## We need large-coverage word embeddings for the context words, so download
## the pre-compiled word2vec model (200k word forms as target terms, 300 dimensions)
load("GoogleNews300_wf200k.rda", verbose=TRUE)
M <- GoogleNews300_wf200k

## The wordspace package includes some WSD examples derived from the SemCor corpus
?SemCorWSD
View(SemCorWSD) # can also view SemCorWSD in RStudio editor
with(SemCorWSD, table(sense, target)) # overview of target words and their senses

## A simple viewer function for the SemCorWSD test items
display <- function (x, n=20, width=80, random=FALSE, show=c("sentence", "hw", "lemma")) {
  show <- match.arg(show)
  sents <- x[[show]]
  sents <- ifelse(nchar(sents) > width, paste0(substr(sents, 1, width-3), "..."), sents)
  formatted <- sprintf("%-14s  %s", x$sense, sents)
  if (length(formatted) > n) {
    idx <- if (random) sort(sample(seq_along(formatted), n)) else 1:n
    formatted <- formatted[idx]
  }
  cat(formatted, sep="\n")
}

display(SemCorWSD, 10) # first ten items
display(subset(SemCorWSD, target=="plant"), 10, random=TRUE) # 10 random items for "plant"
display(subset(SemCorWSD, target=="plant"), 10, random=TRUE, show="lemma")

## Let us illustrate the Schuetze procedure for the noun "vessel"
Vessel <- subset(SemCorWSD, target == "vessel")
display(Vessel)

## There are two senses in the data set, each with 6 example sentences
with(Vessel, table(paste(sense, gloss, sep=": ")))

## The target word forms in M are lowercased, so we need to convert the example sentences
Vessel$sentence
sent.lc <- tolower(Vessel$sentence)

## Compute centroid vectors representing the sentences
##  - context.vectors() splits on whitespace and ignores unknown words
##  - we don't exclude "vessel" from the contexts here,
##    but that shouldn't make much of a difference because it adds the same vector to each centroid
centroids <- context.vectors(M, sent.lc, row.names=Vessel$id) # returns matrix with 12 rows
round(centroids[, 1:5], 3) # first 5 dimensions of centroids

## Schuetze's approach is to cluster the sentence vectors (i.e. the centroids)
## based on their distance matrix
distances <- dist.matrix(centroids) # of course, you can try other distance measures as well

## A wordspace distance matrix can easily be visualized
plot(distances, labels=Vessel$sense) # senses seem to separate quite well
senses <- factor(Vessel$sense) # this will allow us to map senses to colours in the plot
plot(distances, labels=NULL, col=senses, cex=2) # easier to understand than the first plot
legend("topleft", inset=.02, pch=20, col=1:2, pt.cex=2, legend=levels(senses))

## Another possible visualization of the distance matrix is a heatmap plot,
## which rearranges sentences according to their mutual distances
heatmap(distances, symm=TRUE, labRow=Vessel$sense) # suggests two clearly separated clusters

## The dendrograms in the margins show a hierarchical clustering of the sentences
plot(hclust(as.dist(distances))) # hclust() expects 'dist' object rather than 'dist.matrix'

## PAM = partitioning around medoids is a good algorithm for flat clustering
res <- pam(distances, 2, diss=TRUE, keep.diss=TRUE)
plot(res, which.plots=2, main="silhouette plot") # long bars = points that fit well into their clusters
clusplot(res) # visualization of clustering (NB: symbols are clusters, not senses)

clusplot(res, col.p=factor(Vessel$sense)) # colours indicate word sense
table(res$clustering, Vessel$sense) # one sentence has been misassigned

## We will use purity as measure of clustering quality
purity <- function (x, y) {
  d <- table(x, y)  # row of table = sense distribution of cluster
  100 * sum(apply(d, 1, max)) / sum(d)
}
purity(res$clustering, Vessel$sense) # 92% purity


## Your task:
##  - try the same clustering procedure for the other target words
##  - can you write a function or loop to simplify/automate this procedure?
##  - can you visualize all sentences from SemCorWSD in a single semantic map?
##  - try using one of the pre-compiled WP500 models
##    (NB: will need $lemma instead of $sentence, and no lowercasing)

