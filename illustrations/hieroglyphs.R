##
##  Data and plots for the "hieroglpyhs" illustration
##

source("utilities.R")

load("data/bnc_vobj_basic.rda") # verb-object DSM for basic English nouns (BNC corpus)

M <- bnc.vobj.basic$M # extract matrix for selected nouns and verbs
selected.nouns <- c("knife","cat","dog","boat","cup","pig","banana")
selected.verbs <- c("get","see","use","hear","eat","kill")
M2 <- M[selected.nouns, selected.verbs]

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
