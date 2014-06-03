    \name{DSM_VerbNounTriples_BNC}
\alias{DSM_VerbNounTriples_BNC}
\docType{data}
\title{
  Verb-Noun Co-occurrence Frequencies from British National Corpus (wordspace)
}
\description{

  A table of co-occurrence frequency counts for verb-subject and verb-object pairs in the British National Corpus (BNC).
  Subject and object are represented by the respective head noun.  Both verb and noun entries are lemmatized.
  Separate frequency counts are provided for the written and the spoken part of the BNC.

}
\usage{
  
DSM_VerbNounTriples_BNC

}
\format{
  A data frame with 775283 rows and the following columns:
  \describe{
    \item{\code{noun}:}{noun lemma}
    \item{\code{rel}:}{syntactic relation (\code{subj} or \code{obj})}
    \item{\code{verb}:}{verb lemma}
    \item{\code{f}:}{co-occurrence frequency of noun-rel-verb triple in subcorpus}
    \item{\code{mode}:}{subcorpus (\code{written} for the writte part of the BNC, \code{spoken} for the spoken part of the BNC)}
  }
  
}
\details{
  In order to save disk space, triples that occur just once in the respective subcorpus (\emph{hapax legomena}) have been omitted from the table.
  This should not lead to a substantial distortion of the frequency counts for mid- to high-frequency nouns and verbs.
}
\source{
  TODO
}
\references{
  TODO
}
\examples{

## TODO: how to compile different DSMs from the table

}
\keyword{datasets}
