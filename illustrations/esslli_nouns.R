##
##  Data and plots illustrations based on the ESSLLI noun clustering task
##

source("utilities.R") # includes wordspace package

## use verb-object co-occurrences from C&C-parsed BNC
## (can't use data set included in package because it doesn't cover all ESSLLI nouns)
VN <- read.dsm.triplet("gzip -cd ../datasets/orig/bnc_cc_vn_triples.txt.gz | awk -F\"\\t\" '$3 == \"obj\"' | cut -f 1,2,4 |", value.first=TRUE, freq=TRUE, encoding="ascii", verbose=TRUE)


## remove low-frequency items, score with log-transformed sparse log-likelihood
VN <- subset(VN, nnzero >= 10, nnzero >= 5, recursive=TRUE)
VN <- dsm.score(VN, score="simple-ll", transform="log", sparse=TRUE, normalize=TRUE) 
SVD <- dsm.projection(VN, n=100, method="svd") # SVD projection

## nearest neighbours with nice formatting
print.nn <- function (x, name=NA, digits=1) {
  if (is.list(x)) {
    for (i in seq_along(x)) print.nn(x[[i]], name=names(x)[i], digits=digits)
  } else {
    fmt <- sprintf("%s (%g)", names(x), round(x, digits))
    if (!is.na(name)) cat(sprintf("%s: ", name))
    cat(paste(fmt, collapse=", "), "\n\n", sep="")
  }  
}

nearest.neighbours(VN, c("coffee", "plant", "dog", "book", "school"), 20)
nearest.neighbours(SVD, c("coffee", "plant", "dog", "book", "school"), 20)
print.nn(nearest.neighbours(SVD, c("coffee", "plant", "dog", "book", "school"), 20))
print.nn(nearest.neighbours(SVD, c("bucket", "rose", "cat", "coffee"), 20), digits=0)

## semantic map and clustering for nouns from ESSLLI 2008 task
task <- ESSLLI08_Nouns
task <- transform(task, class=factor(class, levels=c("bird", "groundAnimal", "green", "fruitTree", "tool", "vehicle"))) # sensible ordering of categories
nouns <- convert.lemma(task$word, "HWLC")
cats <- levels(task$class)
cols <- rainbow(length(cats), v=.8)


## compute semantic map from SVD-reduced space (works better)
S1 <- SVD[nouns, ]
Proj1 <- dsm.projection(scale(S1, scale=FALSE), method="svd", n=2, with.basis=TRUE) # PCA as semantic map coordinates
coord <- isoMDS(dist.matrix(S1, method="euclidean"), k=2, y=Proj1)$points # use MDS to spread points better

dev.new(width=8, height=6, bg="white")
par(cex=1.2, mar=c(2,2,2,1)+.1, xaxs="r", yaxs="r")

plot(coord, col=cols[task$class], pch=20, cex=1.6, xlab="SVD dim 1", ylab="SVD dim 2", main="Semantic map (V-Obj from BNC)", xlim=c(-2, 2.5), ylim=c(-1.5,1.7), xaxs="i", yaxs="i")
legend("topright", inset=.02, col=cols, pch=20, pt.cex=1.6, legend=cats, bg="white")
.pos <-c(
	2,1,3,3,1,1,1,4,2,4, #  1-10
    3,4,3,4,1,4,2,3,2,3, # 11-20
    2,4,1,4,1,1,4,3,3,3, # 21-30
	3,4,1,1,3,2,4,1,1,3, # 31-40
	3,4,3,3)             # 41-44
text(coord, col="black", labels=nouns, pos=.pos, cex=1.0, font=2)
dev.copy2pdf(file="img/vobj_semantic_map.pdf", bg="white", onefile=FALSE)


## corresponding hierarchical clustering (we cheat a bit with angular distance, which produces nicer clusters)
M3.dist <- dist.matrix(S1, method="cosine", as.dist=TRUE) # standard agglomerative clustering with Euclidean distance
clusters <- hclust(M3.dist, method="complete")

par(cex=1.1, mar=c(1,4,2,1)+.1, xaxs="r")
plot(clusters, labels=labels(M3.dist), hang=-1, font=2, cex=1, main="Clustering of concrete nouns (V-Obj from BNC)", xlab="", ylab="Cluster size")
points(1:44, rep(0,44), pch=20, cex=2, col=cols[task$class[clusters$order]])
dev.copy2pdf(file="img/vobj_clustering.pdf", bg="white", onefile=FALSE)

par.save <- par(lwd=3)
rect.hclust(clusters, k=6, border="#DD0000") # show split into 6 clusters
par(par.save)
dev.copy2pdf(file="img/vobj_clustering_6.pdf", bg="white", onefile=FALSE)

dev.off()


## nearest neighbours for selected nouns (list and neighbourhood plots)
format.NN <- function (x) {
  nn <- names(x)
  cat(paste(sprintf("%s (%.1f)", nn, x), collapse=", "))
  cat("\n")
}

dev.new(width=6, height=6, bg="white")
par(cex=1.1, mar=c(0,0,0,0), xaxt="n", yaxt="n", bty="ns")

format.NN(nearest.neighbours(SVD, "dog", n=30))
plot(nearest.neighbours(SVD, "dog", n=30, dist.matrix=TRUE))
dev.copy2pdf(file="img/neighbourhood_dog.pdf", bg="white")

format.NN(nearest.neighbours(SVD, "coffee", n=30))
plot(nearest.neighbours(SVD, "coffee", n=25, dist.matrix=TRUE), expand=.1)
dev.copy2pdf(file="img/neighbourhood_coffee.pdf", bg="white")

format.NN(nearest.neighbours(SVD, "hand", n=30))
plot(nearest.neighbours(SVD, "hand", n=25, dist.matrix=TRUE))
dev.copy2pdf(file="img/neighbourhood_hand.pdf", bg="white")

format.NN(nearest.neighbours(SVD, "school", n=30))
plot(nearest.neighbours(SVD, "school", n=25, dist.matrix=TRUE), expand=.1)
dev.copy2pdf(file="img/neighbourhood_school.pdf", bg="white")

format.NN(nearest.neighbours(SVD, "trousers", n=30))
plot(nearest.neighbours(SVD, "trousers", n=25, dist.matrix=TRUE))
dev.copy2pdf(file="img/neighbourhood_trousers.pdf", bg="white")

format.NN(nearest.neighbours(SVD, "rage", n=30))
plot(nearest.neighbours(SVD, "rage", n=25, dist.matrix=TRUE), expand=.1)
dev.copy2pdf(file="img/neighbourhood_rage.pdf", bg="white")

## need some more cheating (unreduced space) to see NN for two senses of plant
format.NN(nearest.neighbours(VN, "plant", n=30))
plot(nearest.neighbours(VN, "plant", n=25, dist.matrix=TRUE), expand=.1)
dev.copy2pdf(file="img/neighbourhood_plant.pdf", bg="white")

dev.off()
