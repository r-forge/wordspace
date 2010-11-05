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
            score=c("frequency", "t-score", "z-score", "Dice", "MI", "reweight"), sparse=FALSE,
            transform=c("none", "log", "root", "sigmoid"),
            scale=c("none", "standardize"),
            normalize=FALSE, method="euclidean", p=2)

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
    if not \code{"none"}, scale columns of scored matrix by z-transformation (\code{"standardize"})
  }

  \item{normalize}{
    whether to normalize rows vectors of scored matrix to unit length, according to the norm indicated by \code{method} and \code{p}
  }

  \item{method, p}{
    specific norm used if \code{normalize=TRUE}.  See \code{\link{rowNorms}} for admissible values and details on the corresponding norms
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
