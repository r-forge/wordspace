##
## how to do RI efficiently
##

system.time(Q1 <- matrix(rbinom(nc*nd, 1, .01) * (2*rbinom(nc*nd,1,.5)-1), nc, nd))
#   user  system elapsed 
#  5.980   0.850   6.835 

system.time(Q2 <- apply(matrix(1:nd, ncol=1), 1, function (d) {.idx <- sample.int(nc, 180); x <- numeric(nc); x[.idx] <- 1 - 2*rbinom(180, 1, .5); x} ))
#   user  system elapsed 
#  0.370   0.360   0.726 

system.time({Q2 <- matrix(0, nc, nd); for (d in 1:nd) {.idx <- sample.int(nc, 180); Q2[.idx,d] <- 1 - 2*rbinom(180, 1, .5)} }) # this is faster & more memory-efficient
#   user  system elapsed 
#  0.230   0.000   0.224 

system.time({Qx <- matrix(0, nd, nc); for (d in 1:nd) {.idx <- sample.int(nc, 180); Qx[d,.idx] <- 1 - 2*rbinom(180, 1, .5)} })
# for dense matrix, transposed version is equally fast

system.time(Q3 <- { i <- as.vector(apply(matrix(1:nd, ncol=1), 1, function (d) sort(sample.int(nc, 180)-1L))); p <- 180L*(0:nd); x <- 1 - 2*rbinom(180*nd, 1, .5); new("dgCMatrix", Dim=as.integer(c(nc, nd)), p=p, i=i, x=x)} )
#   user  system elapsed 
#  0.290   0.000   0.289 

system.time(Q4 <- { i <- integer(180*nd); for (d in 0:(nd-1)) { i[(180L*d)+(1:180)] <- sort(sample.int(nc, 180)-1L) }; p <- 180L*(0:nd); x <- 1 - 2*rbinom(180*nd, 1, .5); new("dgCMatrix", Dim=as.integer(c(nc, nd)), p=p, i=i, x=x)} )
#   user  system elapsed 
#  0.290   0.000   0.299 

system.time(PM1 <- M %*% Q1)
system.time(PM2 <- M %*% Q2)
#   user  system elapsed 
# 16.110   0.060  16.191 

system.time(PM3 <- M %*% Q3)
system.time(PM4 <- M %*% Q4)
#   user  system elapsed 
#  0.780   0.050   0.834 

ortho <- diag(180, 1000,1000) # check orthogonality of random projection
norm(ortho, "f") # ca. 5700
norm(crossprod(Q1)-ortho, "f")
norm(crossprod(Q2)-ortho, "f")
norm(crossprod(Q3)-ortho, "f")
norm(crossprod(Q4)-ortho, "f")
# Q1 = 1420, Q2 .. Q4 = 1330

# check whether cosines are preserved in the projection
system.time(C0 <- cosine(M))
#   user  system elapsed 
# 30.410   2.050  32.509 

system.time(C1 <- cosine(PM1))
system.time(C2 <- cosine(PM2))
#   user  system elapsed 
#  5.070   1.590   5.508 

system.time(C3 <- cosine(as.matrix(PM3))) # v. slow without as.matrix()
system.time(C4 <- cosine(as.matrix(PM4)))
#   user  system elapsed 
#   2.84    0.37    1.98 

norm(C0, "f") # ca. 255.5
norm(C1 - C0, "f")
norm(C2 - C0, "f")
norm(C3 - C0, "f")
norm(C4 - C0, "f")
# C1 = 97, C2 .. C4 = 95


##
## read.dsm.triplet() without readr (slower, but a little less wasteful with RAM)
##
read.dsm.triplet <- function (filename, freq=FALSE, value.first=FALSE, tokens=FALSE, sep="\t", quote="", nmax=-1, sort=FALSE, encoding=getOption("encoding"), verbose=FALSE) {
  if (verbose) cat(sprintf("Loading DSM triplets from '%s' ... ", filename))
  is.pipe <- grepl("\\|\\s*$", filename, perl=TRUE)
  if (is.pipe) {
    filename <- sub("\\s*\\|\\s*$", "", filename, perl=TRUE)
    fh <- pipe(filename, encoding=encoding) # don't open connection yet, so scan() always converts to UTF-8
  } else {
    fh <- file(filename, encoding=encoding)
  }
  
  if (tokens) {
    if (freq || value.first) warning("freq= and value.first= options are ignored with tokens=TRUE")
    triplets <- scan(fh, what=list(l1=character(), l2=character()), sep=sep, quote=quote, nmax=nmax, quiet=TRUE, multi.line=FALSE, na.strings="", comment.char="")
    triplets$val <- 1 # will automatically be replicated by dsm()
    freq <- TRUE
  } else {
    format <- if (value.first) list(val=double(), l1=character(), l2=character()) else list(l1=character(), l2=character(), val=double())
    triplets <- scan(fh, what=format, sep=sep, quote=quote, nmax=nmax, quiet=TRUE, multi.line=FALSE, na.strings="", comment.char="")
  }
  close(fh)
  
  if (verbose) cat(sprintf("%.2fM %s\n", length(triplets$l1) / 1e6, if (tokens) "tokens" else "items"))

  dsm(target=triplets$l1, feature=triplets$l2, score=triplets$val, raw.freq=freq, sort=sort, verbose=verbose)
}
