\name{dsm.score}
\alias{dsm.score}
\title{
  Weighting, Scaling and Normalisation of Co-occurrence Matrix (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

dsm.score(model,
          score = c("frequency", "simple-ll", "t-score", "z-score", "Dice", "MI", "tf.idf", "reweight"),
          sparse = TRUE, negative.ok = NA,
          transform = c("none", "log", "root", "sigmoid"),
          scale = c("none", "standardize", "center", "scale"),
          normalize = FALSE, method = "euclidean", p = 2,
          matrix.only = FALSE, update.nnzero = FALSE)

}
\arguments{
  
  \item{model}{
    a DSM model, i.e. an object of class \code{dsm}
  }

  \item{score}{
    the association measure to be used for feature weighting (see \dQuote{Details} below)
  }

  \item{sparse}{
    if \code{TRUE} (the default), compute sparse non-negative association scores (see \dQuote{Details} below).  Non-sparse association scores are only allowed if \code{negative.ok=TRUE}.
  }

  \item{negative.ok}{
    whether operations that introduce negative values into the score matrix (non-sparse association scores, standardization of columns, etc.) are allowed.  The default (\code{negative.ok}) is \code{TRUE} if the co-occurrence matrix \eqn{M} is dense, and \code{FALSE} if it is sparse.
  }
    
  \item{transform}{
    scale transformation to be applied to association scores (see \dQuote{Details} below)
  }

  \item{scale}{
    if not \code{"none"}, standardize columns of scored matrix by z-transformation (\code{"standardize"}), center without rescaling (\code{"center"}, or scale to unit RMS without centering (\code{"scale"})
  }

  \item{normalize}{
    if \code{TRUE} normalize row vectors of scored matrix to unit length, according to the norm indicated by \code{method} and \code{p}
  }

  \item{method, p}{
    norm to be used with \code{normalize=TRUE}.  See \code{\link{rowNorms}} for admissible values and details on the corresponding norms
  }
  
  \item{matrix.only}{
    whether to return updated DSM model (default) or only the matrix of scores (\code{matrix.only=TRUE})
  }
  
  \item{update.nnzero}{
    if \code{TRUE} and a full DSM model is returned, update the counts of nonzero entries in rows and columns according to the matrix of scores (there may be fewer nonzero entries with sparse association scores, or more from dense association scores and/or column scaling)
  }
}
\details{
  \subsection{Association measures}{
    \bold{*** TODO ***}
    
    For the \code{tf.idf} measure, relative document frequencies for feature terms should be provided as a variable \code{df}
    in the model's column information.  Otherwise, nonzero counts (variable \code{nnzero}) will be used and divided by the
    total number of rows of the cooccurrence matrix.
    
    IDF values are calculated as \eqn{\log \frac{N + 1}{\mathit{df} + 1}}, where the discounting avoids division-by-zero errors.
    Term frequencies are always unscaled; use an additional \code{log}-transformation for an effect similar to logarithmic tf.
    
  }
  \subsection{Scale transformations}{
    \bold{*** TODO ***}
  }
}
\value{
  Either an updated DSM model of class \code{dsm} (default) or the matrix of (scaled and normalised) association scores (\code{matrix.only=TRUE}).
  
  Note that updating DSM models may require a substantial amount of temporary memory (because of the way memory management is implemented in \R).  This can be problematic when running a 32-bit build of \R and when dealing with very large DSM models, so it may be better to return only the scored matrix in such cases.
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
