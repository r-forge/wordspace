\name{print.dsm}
\alias{print.dsm}
\title{
  Print Information About DSM Object (wordspace)
}
\description{

  Prints a short summary describing the \code{dsm} object, including the number of rows and columns.

}
\usage{

\method{print}{dsm}(x, \dots)

}
\arguments{

  \item{x}{
    an object of class \code{dsm}
  }

  \item{\dots}{
    all other arguments are silently ignored
  }

}
\details{

  The main purpose of this method is to keep users from accidentally trying to print out the internal data structures of a large DSM object.

  For compatibility with the generic method (and the documentation of \code{\link{print}}), the DSM object is returned invisibly.
  
}
\seealso{

  \code{\link{print}} for the generic method.
  
}
\examples{

## TODO (but needs to compile DSM object first)

}
\author{Stefan Evert (\url{http://purl.org/stefan.evert})}
\keyword{TODO}
