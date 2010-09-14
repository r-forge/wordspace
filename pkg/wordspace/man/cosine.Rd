\name{cosine}
\alias{cosine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cosine(M, M2 = M, angles = FALSE, normalised = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{M2}{
%%     ~~Describe \code{M2} here~~
}
  \item{angles}{
%%     ~~Describe \code{angles} here~~
}
  \item{normalised}{
%%     ~~Describe \code{normalised} here~~
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
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (M, M2=M, angles=FALSE, normalised=FALSE) {
  sim <- M \%*\% t(M2) # should calculate dot products between rows of M and rows of M2
  if (!normalised) {
    norms.M <- p.norm(M, 2) # norms of row vectors (if not normalised yet)
    norms.M2 <- p.norm(M2, 2)
    sim <- sim / outer(norms.M, norms.M2)
  }
  if (angles) {
    stopifnot(all(sim >= -(1+1e-6) & sim <= 1+1e-6)) # check that cosines are in appropriate range
    sim <- pmin(pmax(sim, -1), 1) # clamp to range [-1, 1]
    sim <- acos(sim) / pi * 180 # angles are returned in degrees
  }
  rownames(sim) <- rownames(M)
  colnames(sim) <- rownames(M2)
  return(sim)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
