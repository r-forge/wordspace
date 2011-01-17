\name{dsm.score}
\alias{dsm.score}
\title{
  Scoring, Weighting and Normalisation of Cooccurrence Matrix (wordspace)
}
\description{

  %%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

  dsm.score(model,
            score=c("frequency", "simple-ll", "t-score", "z-score", "Dice", "MI", "reweight"),
            sparse=FALSE,
            transform=c("none", "log", "root", "sigmoid"),
            scale=c("none", "standardize", "center"),
            normalize=FALSE, method="euclidean", p=2,
            matrix.only=FALSE)

}
\arguments{
  
  \item{model}{
  }

  \item{score}{
  }

  \item{sparse}{
  }

  \item{transform}{
  }

  \item{scale}{
    if not \code{"none"}, scale columns of scored matrix by z-transformation (\code{"standardize"}) or center values without scaling (\code{"center"})
  }

  \item{normalize}{
    whether to normalize rows vectors of scored matrix to unit length, according to the norm indicated by \code{method} and \code{p}
  }

  \item{method, p}{
    specific norm used if \code{normalize=TRUE}.  See \code{\link{rowNorms}} for admissible values and details on the corresponding norms
  }
  
  \item{matrix.only}{
    whether to return updated DSM object (default) or only the matrix of scores (\code{TRUE})
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  either an updated DSM object (default) or the matrix of (transformed, scaled, normalised) association scores (if \code{matrix.only=TRUE})
  
  note that updating large DSM objects requires considerably amounts of temporary memory (at least in R 2.12.x), which can be problematic on 32-bit platforms
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
