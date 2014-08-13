\name{SemCorWSD}
\alias{SemCorWSD}
\docType{data}
\title{
  SemCor Word Sense Disambiguation Task (wordspace)
}
\description{

  A collection of sentences containing ambiguous words manually labelled with WordNet senses.
  The data were obtained from the SemCor corpus version 3.0.

}
\usage{

SemCorWSD

}
\format{

  A data frame with 253 rows and the following 8 columns (all of type character):

  \describe{
    \item{\code{id}}{Unique item ID}
    \item{\code{target}}{The target word (lemmatized)}
    \item{\code{pos}}{Word class of the target word (\code{n}, \code{v} or \code{a})}
    \item{\code{sense}}{Sense of the target word in this sentence (given as a WordNet lemma)}
    \item{\code{gloss}}{WordNet definition of this sense}
    \item{\code{sentence}}{The sentence containing the ambiguous word}
    \item{\code{hw}}{Lemmatized form of the sentence (\dQuote{headwords}); punctuation marks are excluded and all remaining words are case-folded}
    \item{\code{lemma}}{Lemmatized and POS-disambiguated form in CWB/Penn format, e.g. \code{move_N} for the headword \emph{move} used as a noun}
  }

}
\details{

  TODO

}
\source{

  TODO (SemCor reference, NLTK extraction)

}
%% \references{
%%  TODO
%% }
\examples{

with(SemCorWSD, table(target, sense))

with(SemCorWSD, sort(unique(paste0(target, " ", sense, ": ", gloss))))

}
\keyword{datasets}
