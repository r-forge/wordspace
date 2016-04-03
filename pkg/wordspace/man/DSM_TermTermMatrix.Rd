\name{DSM_TermTermMatrix}
\alias{DSM_TermTermMatrix}
\docType{data}
\title{
  Example of a Term-Term Co-occurrence Matrix (wordspace)
}
\description{

  This matrix is a typical example of a term-term DSM co-occurrence matrix,
  derived from the English Wikipedia.

}
\usage{

DSM_TermTermMatrix

}
\format{

  A numeric matrix with 7 rows and 7 columns.

  Rows represent the target nouns \emph{cat}, \emph{dog}, \emph{animal}, \emph{time}, \emph{reason}, \emph{cause}, \emph{effect}.

  Columns specify co-occurrence frequencies of these nouns with the words \emph{breed}, \emph{tail}, \emph{feed}, \emph{kill}, \emph{important}, \emph{explain} and \emph{likely} in articles of the English Wikipedia.  Co-occurring words must appear within a distance of at most two word tokens of each other.

}
%% \details{
%% }
%% \source{
%% }
\seealso{

  This matrix covers the same target terms as the term-context matrix \code{\link{DSM_TermContextMatrix}}.

}
\examples{

\dontrun{
plot(hclust(dist.matrix(DSM_TermTermMatrix, as.dist=TRUE)))
}

}
\keyword{datasets}
