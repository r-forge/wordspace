\name{as.distmat}
\alias{as.distmat}
\alias{as.distmat.matrix}
\alias{as.distmat.sparseMatrix}
\title{
  Mark an arbitrary matrix as a pre-computed \code{dist.matrix} object (wordspace)
}
\description{

  Mark an arbitrary dense or sparse matrix as a pre-computed \code{\link{dist.matrix}} object, so it can be used with \code{\link{nearest.neighbours}} and \code{\link{pair.distances}}.  Default methods are provided for regular dense \code{\link{matrix}} objects as well as any type of \code{\link{sparseMatrix}} from the \bold{Matrix} package.

}
\usage{

as.distmat(x, \dots)

\method{as.distmat}{matrix}(x, similarity=FALSE, symmetric=FALSE, \dots)
\method{as.distmat}{sparseMatrix}(x, similarity=FALSE, symmetric=FALSE, force.dense=FALSE, \dots)

}
\arguments{

  \item{x}{a matrix-like object of a suitable class (for which a method implementation is available)}
  
  \item{similarity}{whether the matrix contains similarity or distance values}
  
  \item{symmetric}{whether the distance or similarity is symmetric (i.e. it has the same rows and columns in the same order and \eqn{d(x, y) = d(y, x)}). Methods trust the specified value and do not check whether this is actually true.}
  
  \item{force.dense}{whether to convert a sparse distance matrix into a dense \code{matrix} object. Keep in mind that the resulting matrix may be extremely large.}
  
  \item{...}{additional arguments passed on to the method implementations (see respective manpages for details)}
}
\value{

  If \code{x} is a dense matrix or \code{force.dense=TRUE}, it is assigned to class \code{\link{dist.matrix}} so it can be used with \code{\link{nearest.neighbours}} and \code{\link{pair.distances}} as well as the \code{\link[=plot.dist.matrix]{plot}} and \code{\link[=head.dist.matrix]{head}} methods.
  
  If \code{x} is a sparse matrix, it is marked with an attribute \code{dist.matrix} recognized by \code{\link{nearest.neighbours}} and \code{\link{pair.distances}}; however, method implementations for \code{dist.matrix} objects will not apply.
  
  In either case, attributes \code{similarity} and \code{symmetric} are set as specified.

}
\details{

  This method is called \code{as.distmat} because the regular name \code{as.dist.matrix} would collide with the \code{\link{as.dist}} method for \code{\link{matrix}} objects.

}
\seealso{

   \code{\link[=plot.dist.matrix]{plot}} and \code{\link[=head.dist.matrix]{head}} methods for distances matrices; \code{\link{nearest.neighbours}} and \code{\link{pair.distances}} 
  
}
\examples{

  # interpret co-occurrence frequency as similarity measure
  M <- as.distmat(DSM_HieroglyphsMatrix, similarity=TRUE)
  nearest.neighbours(M, "cat")
  nearest.neighbours(M, "hear", byrow=FALSE)

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{ methods }
