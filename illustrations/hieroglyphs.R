##
##  Data and plots for the "hieroglpyhs" illustration
##

source("utilities.R")

load("data/bnc_vobj_basic.rda") # verb-object DSM for basic English nouns (BNC corpus)
task <- read.delim("data/concrete_nouns_categories.tbl") # concrete nouns (from ESSLLI 2008 task)
task <- transform(task, class1=factor(class1, levels=unique(as.character(class1))))

M <- bnc.vobj.basic$M # extract matrix for selected nouns and verbs
selected.nouns <- c("knife","cat","dog","boat","cup","pig","banana")
selected.verbs <- c("get","see","use","hear","eat","kill")
M2 <- M[selected.nouns, selected.verbs]

row.freqs <- bnc.vobj.basic$f1 # extract marginal frequencies and sample size
col.freqs <- bnc.vobj.basic$f2
N <- bnc.vobj.basic$N

## -- the cooccurrence matrix ---
print(M2)

## -- similarity scores for rows -- 
round(cosine(log(M2 + 1)), 3) # cosine similarities for sparse log-transformed frequencies

## -- 2D geometric illustration --
words2x <- c("cat","dog","knife","boat")
M2x <- M2[words2x, c("get","use")]
len2x <- p.norm(M2x)
print(round(dist(M2x, method="euclidean"), 2)) # distance matrix for illustration below
print(round(dist(norm.rows(M2x), method="euclidean"), 4)) # after normalisation
print(round(cosine(M2x), 4)) # cosine similarities and corresponding angles
print(round(cosine(M2x, angles=TRUE), 2))

# font scaling and other size adjustments for Ubuntu Linux (Mac OS X may differ!)
dev.new(width=6, height=6, bg="white", pointsize=10, type="Xlib") # only Xlib produces sensible PDF copies

par(cex=1.1, mar=c(4,4,2,2)+.1, xaxs="i", yaxs="i") 
plot(M2x, xlim=c(0,125), ylim=c(0,125), pch=21, cex=1.8, bg=1:4, main="Two dimensions of English V-Obj DSM")
text(M2x, labels=words2x, pos=3, cex=1.4)
dev.copy2pdf(file="img/hieroglyph_2d_1.pdf", bg="white", onefile=FALSE) # -- overlay 1

draw.arrow(M2x["dog",1], M2x["dog",2], M2x["cat",1], M2x[c("cat"),2], cut1=5, cut2=5, head1=TRUE, length=.15, label="d = 63.3", label.pos=.6, label.offset=-5, cex=1.2, lwd=2, col="#444444")
draw.arrow(M2x["dog",1], M2x["dog",2], M2x["boat",1], M2x[c("boat"),2], cut1=5, cut2=5, head1=TRUE, length=.15, label="d = 57.5", label.pos=.5, label.offset=-5, cex=1.2, lwd=2, col="#000088")
dev.copy2pdf(file="img/hieroglyph_2d_2.pdf", bg="white", onefile=FALSE) # -- overlay 2 (then restart)

plot(M2x, xlim=c(0,125), ylim=c(0,125), pch=21, cex=1.8, bg=1:4, main="Two dimensions of English V-Obj DSM")
text(M2x, labels=words2x, pos=3, cex=1.4)
.reduce <- 1 - 3 / len2x # reduce length of vector arrows by 3 units
arrows(0, 0, .reduce * M2x[,1], .reduce * M2x[,2], col=1:4, lwd=3, length=.2)
dev.copy2pdf(file="img/hieroglyph_2d_3.pdf", bg="white", onefile=FALSE) # -- overlay 3

curve(sqrt(80^2 - x^2), from=0, to=80, add=TRUE)
points(diag(80/len2x) %*% M2x, pch=1, cex=1.6, lwd=3, col=1:4)
dev.copy2pdf(file="img/hieroglyph_2d_4.pdf", bg="white", onefile=FALSE) # -- overlay 4

.start <- atan2(M2x["dog",2], M2x["dog",1]) # draw arc segment from "dog" to "knife"
.end <- atan2(M2x["knife",2], M2x["knife",1])
.phi <- seq(.start, .end, length.out=100)
.x <- 90*cos(.phi)
.y <- 90*sin(.phi)
lines(.x, .y, lwd=2, lty="22")
arrows(.x[99],.y[99], .x[100],.y[100], lwd=2, lty="22", length=.2)
text(80,56, expression(alpha == 54.3 * degree), cex=1.5, pos=4)
dev.copy2pdf(file="img/hieroglyph_2d_5.pdf", bg="white", onefile=FALSE) # -- overlay 5

dev.off()

## -- clustering and semantic map examples --
M3 <- M[as.character(task$noun), ] # extract full vectors for concrete nouns from ESSLLI task
M3 <- norm.rows(am.score("t", M3, row.freqs, col.freqs, N, sparse=TRUE, log=FALSE)) # gives nice plot
M3 <- scale(M3, scale=FALSE) # center columns (for PCA-style dimensionality reduction)

SVD <- svd.decomp(M3) # transform into SVD space

cats <- levels(task$class1)
cols <- rainbow(length(cats), v=.8)

# font scaling and other size adjustments for Ubuntu Linux (Mac OS X may differ!)
dev.new(width=8, height=6, bg="white", pointsize=10, type="Xlib") # only Xlib produces sensible PDF copies
par(cex=1.1, mar=c(2,2,2,1)+.1, xaxs="i", yaxs="i")

plot(svd.projection(SVD, 2), col=cols[task$class1], pch=20, cex=1.5, xlab="SVD dim 1", ylab="SVD dim 2", main="Semantic map (V-Obj from BNC)", xlim=c(-.5, .8), ylim=c(-.5,.6))
legend("topright", inset=.03, col=cols, pch=20, pt.cex=1.5, legend=cats, bg="white")
.pos <-c(
	4,3,1,3,4,1,3,2,3,1, #  1-10
        3,3,3,3,1,2,4,4,3,3, # 11-20
	1,3,3,3,3,3,3,3,3,1, # 21-30
	3,4,4,3,3,3,2,1,3,3, # 31-40
	1,3,4,4)             # 41-44
text(svd.projection(SVD, 2), col="black", labels=task$noun, pos=.pos, cex=1.0, font=2)
dev.copy2pdf(file="img/hieroglyph_semantic_map.pdf", bg="white", onefile=FALSE)


M3.dist <- dist(SVD$M[,1:2], method="euclidean") # standard agglomerative clustering with Euclidean distance
clusters <- hclust(M3.dist, method="complete")

par(cex=1.2, mar=c(1,4,2,1)+.1, xaxs="r")
plot(clusters, labels=labels(M3.dist), hang=-1, font=2, cex=1.2, main="Word space clustering of concrete nouns (V-Obj from BNC)", xlab="", ylab="Cluster size")
points(1:44, rep(0,44), pch=20, cex=2.5, col=cols[task$class1[clusters$order]])
dev.copy2pdf(file="img/hieroglyph_clustering.pdf", bg="white", onefile=FALSE)

dev.off()

## -- nearest neighbours on full matrix, with neighbourhood graph --
Mx <- norm.rows(am.score("t", M, row.freqs, col.freqs, N, sparse=TRUE, log=FALSE)) # gave nicest semantic map
Mx <- scale(Mx, scale=FALSE)
SVDx <- svd.decomp(Mx, n=100)

dev.new(width=6, height=6, bg="white", pointsize=10, type="Xlib") # only Xlib produces sensible PDF copies
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
