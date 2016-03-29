\name{sparsesvd}
\alias{sparsesvd}
\title{Singular Value Decomposition of a Sparse Matrix.}
\usage{
  sparsesvd(M)
}
\arguments{
  \item{M}{a sparse real matrix in \bold{Matrix} package format.
    The preferred format is a \code{\link[=dgCMatrix-class]{dgCMatrix}} and other storage formats will automatically be converted.
  }
}
\value{
  The truncated SVD decomposition
  \deqn{
    M_r = U_r D V_r^T
  }
  where \eqn{M_r} is the optimal rank \eqn{r} approximation of \eqn{M}.
  
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
  TODO
}
\details{
  TODO
}
\references{
  \url{http://tedlab.mit.edu/~dr/SVDLIBC/}

  Berry, Michael~W. (1992). Large scale sparse singular value computations.
  \emph{International Journal of Supercomputer Applications}, \bold{6}, 13--49.
}
\seealso{
  \code{\link{svd}}, \code{\link{sparseMatrix}}
}

