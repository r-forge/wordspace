##
##  Miscellaneous graphs (e.g. transformation functions)
##

source("utilities.R")

## visualize transformation functions (from ?dsm.score)
dev.new(width=10, height=4)
par(mar=c(4,4,1,1))

seaborn.pal <- c("#4C72B0","#55A868","#C44E52","#8172B2","#CCB974","#64B5CD","black")
x <- seq(-3, 7, .025)

plot(0, 0, type="n", xaxs="i", yaxs="i", xlab="x", ylab="f(x)", xlim=c(-3,7), ylim=c(-1.5,2.5))
abline(h=-2:3, lwd=.5, col="gray60"); abline(v=-3:7, lwd=.5, col="gray60")
abline(h=0, lwd=1); abline(v=0, lwd=1)
legend("topleft", inset=.05, bg="white", lwd=5, col=seaborn.pal[c(3,1,2)],
       legend=c("log", "sigmoid", "sparse"))
dev.copy2pdf(file="img/transformation_functions_0.pdf", onefile=FALSE)

lines(x, sign(x) * log(abs(x) + 1), lwd=5, col=seaborn.pal[3])
dev.copy2pdf(file="img/transformation_functions_1.pdf", onefile=FALSE)

lines(x, tanh(x), lwd=5, col=seaborn.pal[1])
dev.copy2pdf(file="img/transformation_functions_2.pdf", onefile=FALSE)

lines(x, pmax(x, 0), lwd=5, col=seaborn.pal[2])
dev.copy2pdf(file="img/transformation_functions.pdf", onefile=FALSE)

dev.off()

## unit circle according to various p-norms
dev.new(width=6, height=6)
par(cex=1.2, mar=c(2,2,2,2))

unit.circle <- function (p, steps=45, lwd=4, ...) {
  phi <- seq(0, 2, length.out=(8*steps + 1)) # make sure we hit the corners as points
  v <- cbind(cospi(phi), sinpi(phi)) # unit circle
  if (p == Inf) 
  	r <- apply(abs(v), 1, max)
  else
	r <- (rowSums(abs(v) ^ p)) ^ (1/p) # p-norms (even for p < 1)
  v <- scaleMargins(v, rows=(1 / r)) # scale vectors to unit length
  lines(v, lwd=lwd, ...)
}

plot(0, 0, type="n", xlim=c(-1, 1), ylim=c(-1, 1), main="Unit circles for different p-norms")
abline(h=0); abline(v=0)
unit.circle(Inf, col=seaborn.pal[6])
unit.circle(5, col=seaborn.pal[3])
unit.circle(2, col=seaborn.pal[7])
unit.circle(1, col=seaborn.pal[1])
unit.circle(0.5, col=seaborn.pal[2], lty="22")
legend("topleft", inset=.07, bg="white", legend=expression(p == infinity, p == 5, p == 2, p == 1, p == 1/2), lwd=4, col=seaborn.pal[c(6,3,7,1,2)], lty=c(rep("solid", 4), "22"))

dev.copy2pdf(file="img/unit_circle_pnorms.pdf", onefile=FALSE)
dev.off()
