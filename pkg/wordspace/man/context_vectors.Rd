\name{context.vectors}
\alias{context.vectors}
\encoding{utf8}
\title{
  Compute Bag-of-Words Context Vectors (wordspace)
}
\description{

  Compute bag-of-words context vectors as proposed by Schütze (1998) for automatic word sense disambiguation and induction.  Each context vector is the centroid of the DSM vectors of all terms occurring in the context.

}
\usage{

context.vectors(M, contexts, split = "\\\\s+",
                drop.missing = TRUE, row.names=NULL)

}
\arguments{

  \item{M}{
    numeric matrix of row vectors for the terms specified by \code{rownames(M)}, or an object of class \code{dsm}
  }

  \item{contexts}{
    the contexts for which bag-of-words representations are to be computed.  Either a character vector, for which each item is split into a bag of terms, or a list of vectors that can be used as row indices into \code{M}.  Vector representations for all terms in a bag are then looked up in \code{M} and averaged.
  }

  \item{split}{
    Perl regular expression determining how contexts given as a character vector are split into terms.  The default behaviour is to split on whitespace.
  }

  \item{drop.missing}{
    if \code{TRUE} (default), contexts that do not contain any known terms are silently dropped; otherwise the corresponding context vectors will be all zeroes.
  }

  \item{row.names}{
    a character vector of the same length as \code{contexts}, specifying row names for the resulting matrix of centroid vectors
  }
}
\details{

  Bag-of-words context vectors are computed by taking the centroid of the term vectors of all known terms in each context.  Neither word order nor any other structural properties of the contexts are taken into account.

}
\value{

  A numeric matrix with the same number of columns as \code{M} and one row for each context (excluding contexts without known terms if \code{drop.missing=TRUE}).  If the vector \code{contexts} has names or \code{row.names} is specified, the matrix rows will be labelled accordingly.  Otherwise the row labels correspond to the indices of the respective entries in \code{contexts}, so matrix rows can always be identified unambiguously if \code{drop.missing=TRUE}.
  
}
\references{

  Schütze, Hinrich (1998). Automatic word sense discrimination. \emph{Computational Linguistics}, \bold{24}(1), 97--123.  

}
\seealso{

  \code{\link{SemCorWSD}}
  
}
\examples{

# illustration of WSD algorithm: 6 sentences each for two senses of "vessel"
VesselWSD <- subset(SemCorWSD, target == "vessel")
with(VesselWSD, cat(paste0(sense, ": ", sentence, "\n")))

# provide sense labels in case some contexts are dropped b/c of too many missing words
Centroids <- with(VesselWSD, context.vectors(DSM_Vectors, lemma, row.names=sense))
Centroids[, 1:5]

(res <- kmeans(Centroids, 2)$cluster) # flat clustering with k-means
table(rownames(Centroids), res)       # ... works perfectly

\dontrun{
plot(hclust(dist.matrix(Centroids, as.dist=TRUE)))
}
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
