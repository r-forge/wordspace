library(Matrix)
library(RSpectra)
library(svd)
library(irlba)
library(sparsesvd)
library(microbenchmark)
library(dplyr)
library(ggplot2)

n = 5000
p = 2000
rate = 0.1
nu = nv = k = 20

set.seed(42)
x = matrix(rnorm(n * p), n)
x[sample(n * p, floor(n * p * (1 - rate)))] = 0
xsp = as(x, "dgCMatrix")

## For svd package
f = function(v) as.numeric(xsp %*% v)
tf = function(v) as.numeric(crossprod(xsp, v))
extx = extmat(f, tf, nrow(xsp), ncol(xsp))

res = microbenchmark(
    "svd"             = svd(x, nu, nv),
    "svds"            = svds(x, k, nu, nv, opts = list(tol = 1e-8)),
    "propack"         = propack.svd(x, k, opts = list(tol = 1e-8)),
    "trlan"           = trlan.svd(x, k, opts = list(tol = 1e-8)),
    "irlba"           = irlba(x, nu, nv, tol = 1e-8),
    "svds[sparse]"    = svds(xsp, k, nu, nv, opts = list(tol = 1e-8)),
    "propack[sparse]" = propack.svd(extx, k, opts = list(tol = 1e-8)),
    "trlan[sparse]"   = trlan.svd(extx, k, opts = list(tol = 1e-8)),
    "irlba[sparse]"   = irlba(xsp, nu, nv, tol = 1e-8),
    "sparsesvd[sparse]" = sparsesvd(xsp, rank=k),
    times = 3
)

dat = as.data.frame(res)
dat$type = ifelse(grepl("sparse", dat$expr), "Matrix type: Sparse",
                                             "Matrix type: Dense")
dat$expr = factor(gsub("\\[sparse\\]", "", dat$expr),
                  levels = c("svd", "sparsesvd", "irlba", "propack", "trlan", "svds"))
dat = dat %>% group_by(expr, type) %>% summarize(medtime = median(time) / 1e6)
dat$Package = ifelse(grepl("svds", dat$expr), "RSpectra", "Other")

if (!interactive()) cairo_pdf(file="svd_solvers.pdf", width=12, height=6)
ggplot(dat, aes(x = expr, y = medtime)) +
    facet_wrap(~ type, scales = "free", ncol = 2) +
    geom_bar(aes(fill = Package), stat = "identity", width = 0.7) +
    scale_y_continuous("Median Elapsed time (ms)") +
    scale_x_discrete("Functions") +
    ggtitle(sprintf("Top-%d SVD on %dx%d Matrix with %g%% fill rate", k, n, p, 100 * rate)) +
    theme_bw(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5)) |>
    print()
if (!interactive()) dev.off()
