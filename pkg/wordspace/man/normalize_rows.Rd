\name{normalize.rows}
\alias{normalize.rows}
\alias{normalize.cols}
\title{
  Normalize Rows or Columns of Matrix to Unit Length (wordspace)
}
\description{

  Efficiently normalize the row or column vectors of a dense or sparse matrix to unit length.

}
\usage{

normalize.rows(M, \dots)

normalize.cols(M, \dots)

}
\arguments{

  \item{M}{a dense or sparse numeric matrix}

  \item{\dots}{optional \code{method} and \code{p} arguments are passed to \code{\link{rowNorms}} (or \code{\link{colNorms}})}

}
\details{
  Note that the generalized Minkowski length measure for \eqn{p < 1} cannot be normalized with \code{normalize.rows} because it doesn't satisfy the homogeneity axiom.  The Hamming length obtained for \eqn{p = 0} cannot be normalized at all.  See \code{\link{rowNorms}} for more information.
}
\value{

  A row-normalized (or column-normalized) matrix with the same dimensions as \eqn{M}.

}
\seealso{
  See \code{\link{rowNorms}} for details on available norms and their parameters.
}
%% \examples{
%% }
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
