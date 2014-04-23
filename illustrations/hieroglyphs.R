##
##  Data and plots for the "hieroglpyhs" illustration
##

source("utilities.R")

M <- DSM_HieroglyphsMatrix

## -- the cooccurrence matrix ---
print(M)

## -- similarity scores for rows -- 
round(dist.matrix(log(M + 1), method="cosine", convert=FALSE), 3) # cosine similarities for sparse log-transformed frequencies

## -- 2D geometric illustration --
words2x <- c("cat", "dog", "knife", "boat")
Mx <- M[words2x, c("get", "use")]
len2x <- rowNorms(Mx)
print(round(dist(Mx, method="euclidean"), 2)) # distance matrix for illustration below
print(round(dist(normalize.rows(Mx), method="euclidean"), 4)) # after normalisation
print(round(dist.matrix(Mx, method="cosine", convert=FALSE), 4)) # cosine similarities and corresponding angles
print(round(dist.matrix(Mx, method="cosine", convert=TRUE), 2))

dev.new(width=6, height=6, bg="white")
par(cex=1.2, mar=c(4,4,2,2)+.1, xaxs="i", yaxs="i") 

plot(Mx, xlim=c(0,125), ylim=c(0,125), pch=21, cex=1.6, bg=1:4, main="Two dimensions of English V-Obj DSM")
text(Mx, labels=words2x, pos=3, cex=1.2)
dev.copy2pdf(file="img/hieroglyph_2d_1.pdf", bg="white", onefile=FALSE) # -- overlay 1

draw.arrow(Mx["dog",1], Mx["dog",2], Mx["cat",1], Mx[c("cat"),2], cut1=5, cut2=5, head1=TRUE, length=.15, label="d = 63.3", label.pos=.6, label.offset=-5, cex=1.2, lwd=2, col="#444444")
draw.arrow(Mx["dog",1], Mx["dog",2], Mx["boat",1], Mx[c("boat"),2], cut1=5, cut2=5, head1=TRUE, length=.15, label="d = 57.5", label.pos=.5, label.offset=-5, cex=1.2, lwd=2, col="#000088")
dev.copy2pdf(file="img/hieroglyph_2d_2.pdf", bg="white", onefile=FALSE) # -- overlay 2 (then restart)

plot(Mx, xlim=c(0,125), ylim=c(0,125), pch=21, cex=1.6, bg=1:4, main="Two dimensions of English V-Obj DSM")
text(Mx, labels=words2x, pos=3, cex=1.2)
.reduce <- 1 - 3 / len2x # reduce length of vector arrows by 3 units
arrows(0, 0, .reduce * Mx[,1], .reduce * Mx[,2], col=1:4, lwd=3, length=.2)
dev.copy2pdf(file="img/hieroglyph_2d_3.pdf", bg="white", onefile=FALSE) # -- overlay 3

curve(sqrt(80^2 - x^2), from=0, to=80, add=TRUE)
points(diag(80/len2x) %*% Mx, pch=1, cex=1.6, lwd=3, col=1:4)
dev.copy2pdf(file="img/hieroglyph_2d_4.pdf", bg="white", onefile=FALSE) # -- overlay 4

.start <- atan2(Mx["dog",2], Mx["dog",1]) # draw arc segment from "dog" to "knife"
.end <- atan2(Mx["knife",2], Mx["knife",1])
.phi <- seq(.start, .end, length.out=100)
.x <- 90*cos(.phi)
.y <- 90*sin(.phi)
lines(.x, .y, lwd=2, lty="22")
arrows(.x[99],.y[99], .x[100],.y[100], lwd=2, lty="22", length=.2)
text(80,56, expression(alpha == 54.3 * degree), cex=1.4, pos=4)
dev.copy2pdf(file="img/hieroglyph_2d_5.pdf", bg="white", onefile=FALSE) # -- overlay 5

dev.off()
