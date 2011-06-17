\name{dsm.projection}
\alias{dsm.projection}
\title{
  Reduce Dimensionality of DSM by Subspace Projection (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

  dsm.projection(model, method=c("svd", "rsvd", "asvd", "ri", "ri+svd"),
                 n=NA, oversampling=NA, q=2, rate=.01, verbose=FALSE, with.basis=FALSE)

}
\arguments{

  \item{model}{
    a numeric matrix, or  an object of class \code{dsm} that has already been scored/weighted with \code{dsm.score};
  }

  \item{method}{
    projection method to use for dimensionality reduction (see "DETAILS" below)
  }
  
  \item{n}{
    an integer specifying the number of target dimensions.  If unspecified, the maximal number of latent dimensions is used (i.e. the minimum of the number of rows and columns of the DSM matrix)
  }

  \item{oversampling}{
    oversampling factor for stochastic dimensionality reduction algorithms (\code{rsvd}, \code{asvd}, \code{ri+svd}).  If unspecified, the default value is 2 for \code{rsvd}, 10 for \code{asvd} and 10 for \code{ri+svd} (subject to change).
  }

  \item{q}{
    number of power iterations in the randomized SVD algorithm (Halko \emph{et al.} 2009 recommend \code{q=1} or \code{q=2})
  }

  \item{rate}{
    fill rate of random projection vectors.  Each random dimension has exactly \code{rate * ncol(model)} nonzero components in the original space
  }
    
  \item{verbose}{
    if \code{TRUE}, some methods display progress messages during execution
  }
  
  \item{with.basis}{
    if \code{TRUE}, also returns orthogonal basis of the subspace as attribute of the reduced matrix (not available for random indexing methods)
  }
  
}

\details{
  The following dimensionality reduction algorithms can be selected with the \code{method} argument:
  \describe{

    \item{svd}{singular value decomposition (SVD).  Note that special efficient algorithms for sparse matrices are not implemented yet.  If the DSM has been scored with \code{scale="center"}, this method is equivalent to principal component analysis (PCA).}

    \item{rsvd}{randomized SVD (Halko \emph{et al.} 2009, p. 9) based on a factorization of rank \code{oversampling * n} with \code{q} power iterations.}

    \item{asvd}{approximate SVD, which determines latent dimensions from a random sample of matrix rows including  \code{oversampling * n} data points.}

    \item{ri}{random indexing (RI) -- \bold{TODO}}

    \item{ri+svd}{RI to \code{oversampling * n} dimensions, followed by SVD to pick the final \code{n} dimensions -- \bold{TODO}}

  }
  
}

\value{

  A numeric matrix with \code{n} columns (latent dimensions) and the same number of rows as the original DSM.  If \code{with.projection=TRUE}, an orthogonal basis \eqn{B} of the latent subspace is returned as attribute \code{"basis"} of the return value.  \eqn{B} is column-orthogonal, hence \eqn{B^T} projects into latent dimensions and \eqn{B B^T} is an orthogonal projection in the original coordinate system.

}

\references{
  Halko, N., Martinsson, P. G., and Tropp, J. A. (2009).
  Finding structure with randomness: Stochastic algorithms for constructing approximate matrix decompositions. Technical Report 2009-05, ACM, California Institute of Technology.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
