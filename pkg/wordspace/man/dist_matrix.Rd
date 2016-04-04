\name{dist.matrix}
\alias{dist.matrix}
\title{
  Distances/Similarities between Row or Column Vectors (wordspace)
}
\description{

  Compute a symmetric matrix of distances (or similarities) between the rows or columns of a matrix;
  or compute cross-distances between the rows or columns of two different matrices.
  This implementation is faster than \code{\link{dist}} and can operate on sparse matrices (in canonical DSM format).

}
\usage{

dist.matrix(M, M2 = NULL,
            method = c("cosine", "euclidean", "maximum", "manhattan", "minkowski", "canberra"), p = 2, 
            normalized = FALSE, byrow = TRUE, convert = TRUE, as.dist = FALSE, 
            terms = NULL, terms2 = terms, skip.missing = FALSE)

}
\arguments{

  \item{M}{a dense or sparse matrix representing a scored DSM, or an object of class \code{dsm}}
  
  \item{M2}{an optional dense or sparse matrix representing a second scored DSM, or an object of class \code{dsm}.
    If present, cross-distances between the rows (or columns) of \code{M} and those of \code{M2} will be computed.
  }

  \item{method}{distance or similarity measure to be used (see \dQuote{Distance Measures} below for details)}
  
  \item{p}{exponent of the \code{minkowski} \eqn{L_p}-metric, a numeric value in the range \eqn{0 \le p < \infty}.
    The range \eqn{0 \le p < 1} represents a generalization of the Minkowski distance, which cannot be derived from a proper mathematical norm (see details below).
  }
  
  \item{normalized}{if \code{TRUE}, assume that the row (or column) vectors of \code{M} and \code{M2} have been appropriately normalised (depending on the selected distance measure) in order to speed up calculations.
    This option is often used with the \code{cosine} metric, for which vectors must be normalized wrt. the Euclidean norm.  It is currently ignored for other distance measures.
  }

  \item{byrow}{whether to calculate distances between row vectors (default) or between column vectors (\code{byrow=FALSE})}
  
  \item{convert}{if \code{TRUE}, similarity measures are automatically converted to distances in an appropriate way (see \dQuote{Distance Measures} below for details).
    Note that this is the default setting and \code{convert=FALSE} has to be specified explicitly in order to obtain a similarity matrix.
  }
  
  \item{as.dist}{convert the full symmetric distance matrix to a compact object of class \code{dist}.
    This option cannot be used if cross-distances are calculated (with argument \code{M2}) or if a similarity measure has been selected (with option \code{convert=FALSE}).
  }

  \item{terms}{a character vector specifying rows of \code{M} for which distance matrix is to be computed (or columns if \code{byrow=FALSE})}

  \item{terms2}{a character vector specifying rows of \code{M2} for which the cross-distance matrix is to be computed (or columns if \code{byrow=FALSE}).
    If only the argument \code{terms} is specified, the same set of rows (or columns) will be selected from both \code{M} and \code{M2}; you can explicitly specify \code{terms2=NULL} in order to compute cross-distances for all rows (or columns) of \code{M2}.
  }
  
  \item{skip.missing}{if \code{TRUE}, silently ignores terms not found in \code{M} (or in \code{M2}). By default (\code{skip.missing=FALSE}) an error is raised in this case.}

}
%% \details{
%% }
\value{

  By default, a matrix of class \code{dist.matrix}, specifying distances or similarities between term vectors.
  A similarity matrix is marked by an additional attribute \code{similarity} with value \code{TRUE}.
  If the distance or similarity matrix is symmetric (i.e. neither a cross-distance matrix nor based on an asymmetric distance measure), it is marked by an attribute \code{symmetric} with value \code{TRUE}.
  
  If \code{as.dist=TRUE}, the matrix is compacted to an object of class \code{dist}.
}
\section{Distance Measures}{

  Given two DSM vectors \eqn{x} and \eqn{y}, the following distance metrics can be computed:

  \describe{
    \item{\code{euclidean}}{The Euclidean distance given by \deqn{
        d_2(x, y) = \sqrt{ \sum_i (x_i - y_i)^2 }}{ d_2(x, y) = sqrt( SUM(i) (x_i - y_i)^2 ) }
    }
    \item{\code{manhattan}}{The Manhattan (or \dQuote{city block}) distance given by \deqn{
        d_1(x, y) = \sum_i |x_i - y_i|}{ d_1(x, y) = SUM(i) |x_i - y_i| }
    }
    \item{\code{maximum}}{The maximum distance given by \deqn{
        d_{\infty}(x, y) = \max_i |x_i - y_i|}{ d_Inf(x, y) = MAX(i) |x_i - y_i| }
    }
    \item{\code{minkowski}}{The Minkowski distance is a family of metrics determined by a parameter \eqn{0 \le p < \infty}, which encompasses the Euclidean, Manhattan and maximum distance as special cases.  Also known as \eqn{L_p}-metric, it is defined by \deqn{
        d_p(x, y) = \left( \sum_i |x_i - y_i|^p \right)^{1/p}}{ d_p(x, y) = [ SUM(i) |x_i - y_i|^p ]^(1/p)  }
      for \eqn{p \ge 1} and by \deqn{
        d_p(x, y) = \sum_i | x_i - y_i |^p}{ d_p(x, y) = SUM(i) |x_i - y_i|^p }
      for \eqn{0 \le p < 1}.  In the latter case, it is not homogeneous and cannot be derived from a corresponding mathematical norm (cf. \code{\link{rowNorms}}).

      Special cases include the Euclidean metric \eqn{d_2(x, y)} for \eqn{p = 2} and the Manhattan metric \eqn{d_1(x, y)} for \eqn{p = 1}, but the dedicated methods above provide more efficient implementations.  For \eqn{p \to \infty}{p -> Inf}, \eqn{d_p(x, y)} converges to the maximum distance \eqn{d_{\infty}(x, y)}, but \code{p=Inf} is not allowed.  For \eqn{p = 0}, \eqn{d_p(x, y)} corresponds to the Hamming distance, i.e. the number of differences \deqn{
        d_0(x, y) = \#\{ i | x_i \ne y_i \}}{ #{i | x_i != y_i} }
    }
    \item{\code{cosine}}{The angular distance between two vectors can be computed by specifying \code{method="cosine"} in combination with \code{convert=FALSE}.  Angles are given in degrees, ranging from 0 to 180.}
    \item{\code{canberra}}{The Canberra metric has been implemented for compatibility with the \code{\link{dist}} function, even though it is probably not very useful for DSM vectors.  It is given by \deqn{
        \sum_i \frac{|x_i - y_i|}{|x_i| + |y_i|}}{ SUM(i) |x_i - y_i| / |x_i| + |y_i| }
      (see \url{https://en.wikipedia.org/wiki/Canberra_distance}).  Terms with \eqn{x_i = y_i = 0} are silently dropped from the summation.

      Note that \code{\link{dist}} uses a different formula https://en.wikipedia.org/wiki/Canberra_distance \deqn{
        \sum_i \frac{|x_i - y_i|}{|x_i + y_i|}}{ SUM(i) |x_i - y_i| / |x_i + y_i| }
      which is highly problematic unless \eqn{x} and \eqn{y} are guaranteed to be non-negative.  Terms with \eqn{x_i = y_i = 0} are imputed, i.e. set to the average value of all nonzero terms.
    }
  }

  In addition, the following similarity measures can be computed by setting \code{convert=FALSE}:
  \describe{
    \item{\code{cosine}}{The cosine similarity given by the normalized inner product of the vectors \eqn{x} and \eqn{y}.  Cosine similarity is related to the angle \eqn{\alpha} between the vectors by \deqn{
        \cos \alpha = \frac{x^T y}{\|x\|\cdot \|y\|}}{ cos(alpha) = x' y / |x| |y|}
    }
  }
  Keep in mind that a similarity matrix cannot be converted into a distance structure with \code{as.dist=TRUE}.
  
}
%% \references{
%% }
\seealso{
  
  \code{\link{rowNorms}}

}
\examples{

M <- DSM_TermTermMatrix
dist.matrix(M, as.dist=TRUE)                     # angular distance
dist.matrix(M, method="euclidean", as.dist=TRUE) # Euclidean distance
dist.matrix(M, method="manhattan", as.dist=TRUE) # Manhattan distance
dist.matrix(M, method="minkowski", p=1, as.dist=TRUE)  # L_1 distance
dist.matrix(M, method="minkowski", p=99, as.dist=TRUE) # almost L_Inf
dist.matrix(M, method="maximum", as.dist=TRUE)         # L_Inf (maximum)
dist.matrix(M, method="minkowski", p=.5, as.dist=TRUE) # L_0.5 distance
dist.matrix(M, method="minkowski", p=0, as.dist=TRUE)  # Hamming distance

round(dist.matrix(M, method="cosine", convert=FALSE), 3) # cosine similarity

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
