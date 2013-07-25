\name{pair.distances}
\alias{pair.distances}
\title{

  Semantic Distances Between Word Pairs (wordspace)

}
\description{

  Compute semantic distances (or similarities) between pairs of words based on a scored DSM matrix \code{M},
  according to any of the distance measures supported by \code{\link{dist.matrix}}.
  If one of the words in a pair is not represented in the DSM, the distance is set to \code{Inf}
  (or to \code{-Inf} in the case of a similarity measure)

}
\usage{
  pair.distances(w1, w2, M, rank=FALSE, \dots, batchsize=10e6, verbose=FALSE)
}
\arguments{

  \item{w1}{a character vector specifying the first word of each pair}

  \item{w2}{a character vector of the same length as \code{w1}, specifying the second word of each pair}
  
  \item{M}{a scored DSM matrix, suitable for passing to \code{\link{dist.matrix}}}
  
  \item{rank}{if \code{TRUE}, return neighbour rank of second word rather than distance}
  
  \item{\dots}{further arguments are passed to \code{dist.matrix} and determine the distance or similarity measure to be used (see \code{\link{dist.matrix}} for details)}

  \item{batchsize}{maximum number of similarity values to compute per batch. This parameter has an essential influence on efficiency and memory use of the algorithm and has to be tuned carefully for optimal performance.}
  
  \item{verbose}{if \code{TRUE}, display some progress messages indicating how data are split into batches}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  
  If \code{rank=FALSE} (the default), a numeric vector of the same length as \code{w1} and \code{w2}
  specifying the distances between the word pairs, according to the selected distance measure.
  
  If \code{rank=TRUE}, an integer or numeric vector of the same length as \code{w1} and \code{w2}
  specifying the neighbour rank of the second word.
  
  In either case, a distance of \code{Inf} (or similarity of \code{-Inf}) is returned for any word pair not represented in the DSM.

}
\seealso{
  \code{\link{dist.matrix}}
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
