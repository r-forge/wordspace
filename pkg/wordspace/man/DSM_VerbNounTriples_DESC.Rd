\name{DSM_VerbNounTriples_DESC}
\alias{DSM_VerbNounTriples_DESC}
\docType{data}
\title{
  Verb-Noun Co-occurrence Frequencies from the DESC Corpus of English Movie Subtitles (wordspace)
}
\description{

  A table of co-occurrence frequency counts for verb-subject and verb-object pairs in English movie subtitles,
  obtained from the Darmstadt English Subtitles Corpus v2 (DESC).
  Subject and object are represented by the respective head noun.  Both verb and noun entries are lemmatized.

}
\usage{
  
head(DSM_VerbNounTriples_DESC, 20)

}
\format{
  A data frame with 389121 rows and the following columns:
  \describe{
    \item{\code{noun}:}{noun lemma}
    \item{\code{rel}:}{syntactic relation (\code{subj} or \code{obj})}
    \item{\code{verb}:}{verb lemma}
    \item{\code{f}:}{co-occurrence frequency of noun-rel-verb triple}
  }
  
}
\details{
  In order to save disk space, triples that occur just once in the corpus (\emph{hapax legomena}) have been omitted from the table.
  This should not lead to a substantial distortion of the frequency counts for mid- to high-frequency nouns and verbs.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
## TODO: how to compile different DSMs from the table
}
\keyword{datasets}
