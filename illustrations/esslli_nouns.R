##
##  Data and plots illustrations based on the ESSLLI noun clustering task
##

source("utilities.R")

load("data/bnc_vobj_basic.rda") # verb-object DSM for basic English nouns (BNC corpus)
task <- read.delim("data/concrete_nouns_categories.tbl") # concrete nouns (from ESSLLI 2008 task)

bnc_vobj_basic <- dsm.score(bnc_vobj_basic, score="t-score", sparse=TRUE, normalize=TRUE) # score & normalize row vectors

M <- bnc_vobj_basic$M # raw co-occurrence frequencies
S <- bnc_vobj_basic$S # scored and normalized row vectors

## -- clustering and semantic map examples --
S1 <- S[as.character(task$noun), ] # scored/normalized vectors for concrete nouns from ESSLLI task
S1 <- scale(S1, scale=FALSE) # center columns (for PCA-style dimensionality reduction)

Proj1 <- dsm.projection(S1, method="svd", n=2, with.basis=TRUE)

cats <- levels(task$class1)
cols <- rainbow(length(cats), v=.8)

## TODO: plots need complete revisions because SVD projection has changed (even with centering)
##  - should use opportunity to switch to a sensible DSM instead (as standard model to be include in package)
##  - compute semantic map with sammon() or so rather than SVD dimensions?

dev.new(width=8, height=6, bg="white")
par(cex=1.2, mar=c(2,2,2,1)+.1, xaxs="r", yaxs="r")

plot(Proj1, col=cols[task$class1], pch=20, cex=1.6, xlab="SVD dim 1", ylab="SVD dim 2", main="Semantic map (V-Obj from BNC)", xlim=c(-.5, .7), ylim=c(-.4,.6))
legend("topright", inset=.02, col=cols, pch=20, pt.cex=1.6, legend=cats, bg="white")
.pos <-c(
	4,3,1,3,4,1,3,2,3,1, #  1-10
        3,3,3,3,1,2,4,4,3,3, # 11-20
	1,3,3,3,3,3,1,3,3,1, # 21-30
	3,4,4,3,3,3,2,1,3,3, # 31-40
	1,3,4,4)             # 41-44
text(Proj1, col="black", labels=task$noun, pos=.pos, cex=1.0, font=2)
dev.copy2pdf(file="img/hieroglyph_semantic_map.pdf", bg="white", onefile=FALSE)


M3.dist <- dist(SVD$M[,1:2], method="euclidean") # standard agglomerative clustering with Euclidean distance
clusters <- hclust(M3.dist, method="complete")

par(cex=1.1, mar=c(1,4,2,1)+.1, xaxs="r")
plot(clusters, labels=labels(M3.dist), hang=-1, font=2, cex=1, main="Word space clustering of concrete nouns (V-Obj from BNC)", xlab="", ylab="Cluster size")
points(1:44, rep(0,44), pch=20, cex=2, col=cols[task$class1[clusters$order]])
dev.copy2pdf(file="img/hieroglyph_clustering.pdf", bg="white", onefile=FALSE)

dev.off()

## -- nearest neighbours on full matrix, with neighbourhood graph --
Mx <- norm.rows(am.score("t", M, row.freqs, col.freqs, N, sparse=TRUE, log=FALSE)) # gave nicest semantic map
Mx <- scale(Mx, scale=FALSE)
SVDx <- svd.decomp(Mx, n=100)

dev.new(width=6, height=6, bg="white")
par(cex=1.1, mar=rep(0,4), xaxt="n", yaxt="n", xaxs="i", yaxs="i")

neighbours(Mx, "dog", 40)
neighbours(SVDx$M, "dog", 30, plot="mds")
dev.copy2pdf(file="img/neighbourhood_dog.pdf", bg="white")

neighbours(SVDx$M, "coffee", 40)
neighbours(SVDx$M, "coffee", 20, plot="sammon", edges=TRUE)
dev.copy2pdf(file="img/neighbourhood_coffee.pdf", bg="white")

neighbours(SVDx$M, "hand", 30, plot="sammon", edges=TRUE)
dev.copy2pdf(file="img/neighbourhood_hand.pdf", bg="white")

neighbours(SVDx$M, "school", 30, plot="sammon", edges=TRUE)
dev.copy2pdf(file="img/neighbourhood_school.pdf", bg="white")

neighbours(SVDx$M, "trousers", 30, plot="sammon", edges=TRUE)
dev.copy2pdf(file="img/neighbourhood_trousers.pdf", bg="white")

neighbours(SVDx$M, "school", 30, plot="sammon", edges=TRUE, iterate=TRUE)
neighbours(Mx, "coffee", 20, plot="sammon", iterate=TRUE, edges=TRUE, aspect=1)
neighbours(SVDx$M, "trousers", 30, plot="sammon", edges=TRUE, iterate=TRUE)

dev.off()
