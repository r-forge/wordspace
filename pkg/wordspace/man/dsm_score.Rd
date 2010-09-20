\name{dsm.score}
\alias{dsm.score}
\title{
TODO
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dsm.score(model, score=c("frequency", "t-score","z-score","Dice"), sparse=FALSE, transform=c("none", "log", "root", "sigmoid"))
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
function (AM=c("frequency", "t-score","z-score","Dice"), M, row.freqs, col.freqs, N, sparse=FALSE, log=FALSE) {
  AM <- match.arg(AM)
  f1 <- row.freqs[rownames(M)]    # first marginal frequencies (rows)
  stopifnot(all(!is.na(f1)))
  f2 <- col.freqs[colnames(M)]    # second marginal frequencies (columns)
  stopifnot(all(!is.na(f2)))

  O <- M # observed frequencies
  if (sparse) {
    scores <- matrix(0, nrow(M), ncol(M))
    nonzero <- M > 0
  }

  if (AM == "t-score") {
    E <- outer(f1, f2) / N # expected frequencies
    if (sparse) {
      idx <- nonzero | O > E
      scores[idx] <- (O[idx] - E[idx]) / sqrt(O[idx])
    }
    else {
      scores <- (O - E) / sqrt(O + 1) # "discounted" t-score
    }
  } else if (AM == "z-score") {
    E <- outer(f1, f2) / N # expected frequencies
    if (sparse) {
      idx <- nonzero | O > E
      scores[idx] <- (O[idx] - E[idx]) / sqrt(E[idx])
    }
    else {
      scores <- (O - E) / sqrt(E)
    }
  } else if (AM == "Dice") {
    scores <- 2 * O / outer(f1, f2, "+")
  } else {
    scores <- O # "frequency" measure = observed frequency
  }

  if (log) scores <- sign(scores) * log(abs(scores) + 1) # "signed log" transformation

  rownames(scores) <- rownames(M) # make sure that row and column names are preserved
  colnames(scores) <- colnames(M)
  return(scores)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
