\name{normalize.rows}
\alias{normalize.rows}
\title{
  Normalize Rows of Matrix to Unit Length (wordspace)
}
\description{

  Efficiently normalize the row vectors of a dense or sparse matrix to unit length.

}
\usage{

normalize.rows(M, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{M}{
    a dense or sparse numeric matrix
  }

  \item{\dots}{
    optional \code{method} and \code{p} arguments are passed to \code{\link{rowNorms}}
  }

}
\details{
  Note that the generalized Minkowski length measure for \eqn{p < 1} cannot be normalized with \code{normalize.rows} because it doesn't satisfy the homogeneity axiom.  The Hamming length obtained for \eqn{p = 0} cannot be normalized at all.  See \code{\link{rowNorms}} for more information.
}
\value{

  A row-normalized matrix with the same dimensions as \eqn{M}.

}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
  See \code{\link{rowNorms}} for details on available norms and their parameters.
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
