##
## Code examples from tutorial part 4
##

library(wordspace)
options(digits=3)

M <- DSM_TermTerm$M
M
M[2, ] # row vector m2 = "dog"
M[, 5] # column vector for "important"

r <- DSM_TermTerm$rows$f
c <- DSM_TermTerm$cols$f
N <- DSM_TermTerm$globals$N

t(r)    # r' = row vector
t(t(r)) # r'' = column vector
r       # plain vector != row/column vector in R

log(M + 1) # discounted log frequency
(M["cause", ] + M["effect", ]) / 2 # mean of cause and effect

outer(r, c) / N # matrix of expected frequencies

(E <- outer(r, c) / N)      # all 3 should be equivalent
(E <- (r %*% t(c)) / N)     # %*% is the matrix product in R
(E <- tcrossprod(r, c) / N) # crossprod(A, B) = t(A) %*% B, tcrossprod(A, B) = A %*% t(B)

log2(M / E)
S <- pmax(log2(M / E), 0)   # = PPMI; pmax = "parallel max"s
S

x <- S[2, ]
b <- sqrt(sum(x ^ 2))
x0 <- x / b
sqrt(sum(x0 ^ 2))  # shows that x0 is normalized

b <- sqrt(rowSums(S^2))
b <- rowNorms(S, method="euclidean")  # more efficient (and less RAM)

S0 <- diag(1 / b) %*% S
S0 <- scaleMargins(S, rows=(1 / b))   # much more efficient (and less RAM)

S0 <- normalize.rows(S, method="euclidean") # convenience function

S0 %*% t(S0)
tcrossprod(S0) # same as the above, but may be more efficient
dist.matrix(S0, convert=FALSE) # check that inner products = cosine similarities

sim <- tcrossprod(S0)
angles <- acos(pmin(sim, 1)) * (180 / pi) # pmin needed in case sim > 1 because of rounding errors

b <- t(t(c(1, 1, 1, 1, .5, 0, 0))) # column basis vector for "animal" subspace
b <- normalize.cols(b) # basis vectors must be normalized

(x <- M %*% b)   # projection of data points into subspace coordinates
x %*% t(b)       # projected points in original space
tcrossprod(x, b) # NB: outer() only works for plain vectors

P <- b %*% t(b)  # projection operator
P # note that P is symmetric
M %*% P          # projected points in original space

fact <- svd(S0)  # SVD decomposition of S0 
round(fact$u, 3) # left singular vectors (columns) = U
round(fact$v, 3) # right singular vectors (columns) = V
round(fact$d, 3) # singular values = diagonal of Sigma
## note that S0 has effective rank 6 because last sigma is practically zero
barplot(fact$d ^ 2) # R2 contributions

r <- 2  # truncated rank-2 SVD
(U.r <- fact$u[, 1:r])
(Sigma.r <- diag(fact$d[1:r], nrow=r))
(V.r <- fact$v[, 1:r])

(X.r <- S0 %*% V.r) # project into latent coordinates
U.r %*% Sigma.r     # same result
scaleMargins(U.r, cols=fact$d[1:r]) # the wordspace way (faster, less RAM)

rownames(X.r) <- rownames(S0) # NB: will want to keep row labels

S0r <- U.r %*% Sigma.r %*% t(V.r) # rank-2 matrix approximation
round(S0r, 3)       # now compare with S0: where are the differences?

round(X.r %*% t(V.r), 3) # same result

## not a PCA because we didn't center the matrix columns
res <- prcomp(S0, rank=2)
res$x              # latent coordinates
res$sdev          # singular values (compare with SVD above)
barplot(res$sdev^2)
S0pca <- res$x %*% t(res$rotation) # matrix approximation
round(S0pca + outer(rep(1, 7), res$center), 3) # undo centering to match original matrix
