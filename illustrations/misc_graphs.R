##
##  Miscellaneous graphs (e.g. transformation functions)
##

source("utilities.R")

## visualize transformation functions (from ?dsm.score)
dev.new(width=10, height=4)
par(mar=c(4,4,1,1))

seaborn.pal <- c("#4C72B0","#55A868","#C44E52","#8172B2","#CCB974","#64B5CD")
x <- seq(-3, 7, .025)

plot(0, 0, type="n", lwd=3, col=seaborn.pal[2], xaxs="i", yaxs="i", xlab="x", ylab="f(x)", xlim=c(-3,7), ylim=c(-1.5,2.5))
abline(h=-2:3, lwd=.5, col="gray60"); abline(v=-3:7, lwd=.5, col="gray60")
abline(h=0, lwd=1); abline(v=0, lwd=1)
lines(x, tanh(x), lwd=3, col=seaborn.pal[1])
lines(x, pmax(x, 0), lwd=3, col=seaborn.pal[2])
lines(x, sign(x) * log(abs(x) + 1), lwd=3, col=seaborn.pal[3])

legend("topleft", inset=.05, bg="white", lwd=4, col=seaborn.pal[c(2,1,3)],
       legend=c("sparse", "sigmoid", "log"))
dev.copy2pdf(file="img/transformation_functions.pdf", onefile=FALSE)