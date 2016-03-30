\name{sparsesvd}
\alias{sparsesvd}
\title{Singular Value Decomposition of a Sparse Matrix.}
\usage{
  sparsesvd(M, rank=0L, tol=1e-15, kappa=1e-6)
}
\arguments{
  \item{M}{a sparse real matrix in \bold{Matrix} package format.
    The preferred format is a \code{\link[=dgCMatrix-class]{dgCMatrix}} and other storage formats will automatically be converted if possible.
  }
  \item{rank}{an integer specifying the desired number of singular components, i.e. the rank of the truncated SVD.
    Specify 0 to return all singular values of magnitude larger than \code{tol} (default).
  }
  \item{tol}{exclude singular values whose magnitude is smaller than \code{tol}}
  \item{kappa}{accuracy parameter \eqn{\kappa} of the SVD algorithm (with SVDLIBC default)}
}
\value{
  The truncated SVD decomposition
  \deqn{
    M_r = U_r D V_r^T
  }
  where \eqn{M_r} is the optimal rank \eqn{r} approximation of \eqn{M}.
  Note that \eqn{r} may be smaller than the requested number \code{rank} of singular components.
  
  The returned value is a list with components
  \item{d}{
    a vector containing the first \eqn{r} singular values of \code{M}
  }
  \item{u}{
    a column matrix of the first \eqn{r} left singular values of \code{M}
  }
  \item{v}{
    a column matrix of the first \eqn{r} right singular values of \code{M}
  }
}
\description{
  Compute the (usually truncated) singular value decomposition (SVD) of a sparse real matrix.
  This function is a shallow wrapper around the SVDLIBC implementation of Berry's (1992) single Lanczos algorithm.
}
\references{
  \url{http://tedlab.mit.edu/~dr/SVDLIBC/}

  Berry, Michael~W. (1992). Large scale sparse singular value computations.
  \emph{International Journal of Supercomputer Applications}, \bold{6}, 13--49.
}
\seealso{
  \code{\link{svd}}, \code{\link{sparseMatrix}}
}
\examples{
M <- rbind(
  c(20, 10, 15,  0,  2),
  c(10,  5,  8,  1,  0),
  c( 0,  1,  2,  6,  3),
  c( 1,  0,  0, 10,  5))
M <- Matrix::Matrix(M, sparse=TRUE)
print(M)

res <- sparsesvd(M, rank=2L) # compute first 2 singular components
print(res, digits=3)

M2 <- res$u \%*\% diag(res$d) \%*\% t(res$v) # rank-2 approximation
print(M2, digits=1)

print(as.matrix(M) - M2, digits=2) # approximation error
}
