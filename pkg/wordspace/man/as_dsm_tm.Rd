\name{as.dsm.tm}
\alias{as.dsm.TermDocumentMatrix}
\alias{as.dsm.DocumentTermMatrix}
\title{
  Create DSM Object From \code{\link{tm}} Package (wordspace)
}
\description{

  Convert a \code{\link{tm}} term-document or document-term matrix into a \code{wordspace} DSM object.

}
\usage{

\method{as.dsm}{TermDocumentMatrix}(obj, \dots, verbose=FALSE)
\method{as.dsm}{DocumentTermMatrix}(obj, \dots, verbose=FALSE)

}
\arguments{

  \item{obj}{
    an term-document or document-term matrix from the \code{\link{tm}} package, i.e.
    an object of a class \code{\link{TermDocumentMatrix}} or \code{\link{DocumentTermMatrix}}.
  }
  
  \item{...}{
    additional arguments are ignored
  }
  
  \item{verbose}{
    if \code{TRUE}, a few progress and information messages are shown
  }
}
\value{

  An object of class \code{\link{dsm}}.

}
\seealso{

  TODO
  
}
\examples{
## TODO
}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
