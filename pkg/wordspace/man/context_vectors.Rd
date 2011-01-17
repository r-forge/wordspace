\name{context.vectors}
\alias{context.vectors}
\title{
  Compute Bag-of-Words Context Vectors (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

  context.vectors(M, contexts, split="\\s+", drop.missing=TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{M}{
    numeric matrix of term vectors for the terms specified by \code{rownames(M)}
  }

  \item{contexts}{
    character vector containing the contexts for which bag-of-words representations are to be computed.  Each context is split into a bag of terms, whose representations are then looked up in \code{M} and averaged.
  }

  \item{split}{
    Perl regular expression determining how contexts are split into terms.  The default behaviour is to split on whitespace.
  }

  \item{drop.missing}{
    if \code{TRUE} (default), contexts that do not contain any known terms are silently dropped; otherwise the corresponding context vectors will be all zeroes.
  }

}
\details{

  Bag-of-words context vectors are computed by taking the centroid of the term vectors of all known terms in each context.  Neither word order nor any other structural properties of the contexts are taken into account.

}
\value{

  A numeric matrix with the same number of columns as \code{M} and one row for each context (excluding contexts without known terms if \code{drop.missing=TRUE}).  If the vector \code{contexts} has names, the matrix rows will be labelled accordingly.  Otherwise the row labels correspond to the indices of the corresponding contexts in \code{contexts}, so matrix rows can be identified unambiguously if \code{drop.missing=TRUE}.
  
}
\references{
%% ~put references to the literature/web site here ~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
