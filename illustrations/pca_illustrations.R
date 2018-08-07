##
##  Illustrations of PCA dimensionality reduction
##

source("utilities.R")


######################################################################
## buy-sell illustration of PCA-based dimensionality reduction

## use full verb-object data from BNC, on which DSM_VerbNounTriples_BNC is based
VN_BNC <- read.delim("../datasets/orig/bnc_cc_vn_triples.txt.gz", quote="",
                     header=FALSE, col.names=qw("f noun rel verb mode"))

VObj <- with(subset(VN_BNC, rel == "obj"),
             dsm(target=noun, feature=verb, score=f, raw.freq=TRUE))

## select buy & sell dimensions and nouns that occur with both verbs
BS <- subset(VObj, select=(term %in% qw("buy sell")), update.nnzero=TRUE)
table(BS$rows$nnzero)
BS <- subset(BS, nnzero >= 2 & rowSums(BS$M) >= 25)

## score with SketchEngine logDice measure (slightly shifted)
BS <- dsm.score(BS, score=function(O, R1, C1, ...) 14 + log2(2 * O / (R1 + C1)))

## now extract coordinates and target nouns
x <- BS$S[, 1]
y <- BS$S[, 2]
nouns <- BS$rows$term

## stats and examples for slides
dim(BS)
BSD <- data.frame(noun=nouns, buy=x, sell=y, row.names=NULL)
print(BSD[round(seq(5, 165, 16)), ], digits=3)

## suitable square graphics device
quartz(width=7, height=7)
par(cex=1.3,mar=c(4,4,0,0)+.2)

## helper function: plot x/y-axes as large arrows
plot.axes <- function (col="black", lwd=2, start=-2, end=4) {
  arrows(start,0, end,0, lwd=lwd, col=col)
  arrows(0,start, 0,end, lwd=lwd, col=col)
}

## plot the data as text labels to motivate latent dimensions
plot(0,0, type="n", xlim=c(0, 10), ylim=c(0, 10), xlab="buy", ylab="sell")
plot.axes(col="#0000AA", start=-.2, end=8)
text(x, y, labels=nouns, cex=0.7, font=2, col=rgb(0,0,0,.6)) 
dev.copy2pdf(file="img/buy_sell_labels_only.pdf", onefile=FALSE, out.type="cairo")

arrows(0,0, 9, 9, lwd=4, col="red") # the latent dimension
dev.copy2pdf(file="img/buy_sell_labels_latent.pdf", onefile=FALSE, out.type="cairo")

## helper functions for visualization
show.data <- function (x, y, xlim=c(-5.5, 4.5), ylim=c(-5.5, 4.5), axis.end=4, axis.start=-axis.end) {
  plot(0,0, type="n", xlim=xlim, ylim=ylim, xlab="buy", ylab="sell")
  plot.axes(col="#0000AA", start=axis.start, end=axis.end)
  points(x, y, pch=20)
}

## add label showing (total or projected) variance
print.variance <- function (variance, digits=2, x=-5.5, y=-5, pos=4, cex=1.3, col="black") {
  text(x, y, pos=pos, cex=cex, col=col, labels=paste("variance =", round(variance, digits)))
}

## p-norm of a vector (for scaling to unit vectors)
norm <- function (x, p=2) { sum(abs(x) ^ p) ^ (1/p) }

## visualise projection onto line (one-dimensional subspace) with arbitrary slope
show.projection <- function (x, y, v, basename=NULL) {
  v <- v / norm(v) # unit vector (direction of line)
  show.data(x, y)
  abline(0, v[2]/v[1], lwd=3)
  if (!is.null(basename)) dev.copy2pdf(file=paste(basename, "axis.pdf", sep="_"), onefile=FALSE, out.type="cairo")
  offsets <- cbind(x, y) %*% v
  projected <- offsets %*% t(v)
  segments(x, y, projected[,1], projected[,2], lwd=1)
  points(projected, pch=20, col="red")
  print.variance(var(offsets))  # show residual variance on line
  if (!is.null(basename)) dev.copy2pdf(file=paste(basename, "projection.pdf", sep="_"), onefile=FALSE, out.type="cairo")
}

## original uncentered data
show.data(x, y, xlim=c(-5, 10), ylim=c(-5, 10), axis.start=-4, axis.end=8)
dev.copy2pdf(file="img/buy_sell_0_uncentered.pdf", onefile=FALSE, out.type="cairo")

## after centring
cx <- x - mean(x)
cy <- y - mean(y)
show.data(cx, cy, xlim=c(-5, 10), ylim=c(-5, 10), axis.start=-4, axis.end=8)
dev.copy2pdf(file="img/buy_sell_0_centered.pdf", onefile=FALSE, out.type="cairo")

## with variance = distance information
show.data(cx, cy)
print.variance( sum(diag(var(cbind(cx, cy)))) ) # 2-d variance = trace of covariance matrix
dev.copy2pdf(file="img/buy_sell_1_variance.pdf", onefile=FALSE, out.type="cairo")

## show projections into several one-dimensional subspaces
show.projection(cx, cy, c(3,-2), basename="img/buy_sell_2") # small variance
show.projection(cx, cy, c(3, 1), basename="img/buy_sell_3")  # medium variance
show.projection(cx, cy, c(2, 3), basename="img/buy_sell_4")  # large variance

## show the optimal PCA dimensions 
pca <- prcomp(cbind(cx, cy))
pca$rotation <- pca$rotation * sign(pca$rotation[1, 1]) # so PC1 points into first quadrant
b1 <- pca$rotation[, 1]
sd1 <- pca$sdev[1]
b2 <- pca$rotation[, 2]
sd2 <- pca$sdev[2]

show.projection(cx, cy, b1, basename="img/buy_sell_5")

## now with PCA axes and labels
show.data(cx, cy)
abline(0, b1[2]/b1[1])
abline(0, b2[2]/b2[1])
arrows(0, 0, sd1*b1[1], sd1*b1[2], lwd=4, col="red")
arrows(0, 0, sd2*b2[1], sd2*b2[2], lwd=4, col="red")
dev.copy2pdf(file="img/buy_sell_6_pca.pdf", onefile=FALSE, out.type="cairo")

label.nouns <- qw("soul story liquor idea information people year man way coat hat supply present goods product share house ticket car land copy asset business stock hat dress pair bottle food book home")
label.pos <-   qw("3    3     3      3    2           2      3    2   1   1    2   4      1       3     2       4     3     4      4   3    2    2     4        2     2   1     4    4      4    1    4")

idx <- match(label.nouns, nouns)
label.pos <- as.integer(label.pos)
                  
points(cx[idx], cy[idx], pch=20, col="red")
text(cx[idx], cy[idx], labels=label.nouns, pos=label.pos, cex=.9, offset=.3, font=2)
dev.copy2pdf(file="img/buy_sell_6_pca_labels.pdf", onefile=FALSE, out.type="cairo")

dev.off()
