##
##  Miscellaneous graphs (e.g. transformation functions)
##

source("utilities.R")


######################################################################
## visualize transformation functions (from ?dsm.score)
dev.new(width=10, height=4)
par(mar=c(4,4,1,1))

seaborn.pal <- c("#4C72B0","#55A868","#C44E52","#8172B2","#CCB974","#64B5CD","black")
x <- seq(-3, 7, .025)

show.plot <- function (step) {
  idx <- seq_len(step)
  plot(0, 0, type="n", xaxs="i", yaxs="i", xlab="x", ylab="f(x)", xlim=c(-3,7), ylim=c(-1.5,2.5))
  abline(h=-2:3, lwd=.5, col="gray60"); abline(v=-3:7, lwd=.5, col="gray60")
  abline(h=0, lwd=1); abline(v=0, lwd=1)
  if (step >= 1) legend("topleft", inset=.05, bg="white", lwd=5,
                        col=seaborn.pal[c(3,1,2,2)][idx], lty=c("solid","solid","solid","22")[idx],
                        legend=c("log", "sigmoid", "sparse", "shifted")[idx])
  if (step >= 1) lines(x, sign(x) * log(abs(x) + 1), lwd=5, col=seaborn.pal[3])
  if (step >= 2) lines(x, tanh(x), lwd=5, col=seaborn.pal[1])
  if (step >= 3) lines(x, pmax(x, 0), lwd=5, col=seaborn.pal[2])
  if (step >= 4) lines(x, pmax(x-3, 0), lwd=5, col=seaborn.pal[2], lty="22")
}

show.plot(0)
dev.copy2pdf(file="img/transformation_functions_0.pdf", onefile=FALSE, out.type="cairo")

show.plot(1)
dev.copy2pdf(file="img/transformation_functions_1.pdf", onefile=FALSE, out.type="cairo")

show.plot(2)
dev.copy2pdf(file="img/transformation_functions_2.pdf", onefile=FALSE, out.type="cairo")

show.plot(3)
dev.copy2pdf(file="img/transformation_functions_3.pdf", onefile=FALSE, out.type="cairo")

show.plot(4)
dev.copy2pdf(file="img/transformation_functions.pdf", onefile=FALSE, out.type="cairo")

dev.off()


######################################################################
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


######################################################################
## scaling of latent dimensions
sigma <- DSM_SingularValues^2 # typical singular values of term-term matrix

sigma.plot <- function (alpha=1, range=NULL, n=1000, main=sprintf("alpha=%.1f", alpha), xlim=c(1, n), ylim=c(0, 1)) {
  s <- sigma[1:n]
  s <- s^alpha # relative weights of dimensions with power scaling
  s <- 100 * s / sum(s)
  k <- if (is.null(range)) seq_along(s) else seq(range[1], range[2])
  plot(0, 0, type="n", xlim=xlim, ylim=ylim, log="", xaxs="i", yaxs="i", bty="l", main=main,
       xlab="latent SVD dimensions", ylab="amount of distance information (%)")
  polygon(c(k, max(k), min(k)), c(s[k], 0, 0), col="#C44E52AA", border="#C44E52AA")
  lines(k, s[k], lwd=4, col="#C44E52")
}

dev.new(width=10, height=4)
par(mar=c(4,4,2,1), cex=1.2)

sigma.plot(alpha=1, xlim=c(1, 400), main=expression(paste("typical singular values ", sigma)))
dev.copy2pdf(file="img/singular_values_1.pdf", out.type="cairo")

sigma.plot(alpha=1, range=c(50, 500), xlim=c(0, 400), main=expression("skip first 50 dimensions"))
dev.copy2pdf(file="img/singular_values_2.pdf", out.type="cairo")

sigma.plot(alpha=1/2, xlim=c(0, 400), main=expression(paste("power scaling ", P == 1/2)))
dev.copy2pdf(file="img/singular_values_3.pdf", out.type="cairo")

sigma.plot(alpha=0, xlim=c(0, 400), main=expression(paste("power scaling ", P == 0)))
dev.copy2pdf(file="img/singular_values_4.pdf", out.type="cairo")

dev.off()
