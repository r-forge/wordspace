## Practice session for part 2:
## How to read your own co-occurrence data into 'wordspace'

library(wordspace)

## The most convenient input format is a sparse co-occurrence matrix in triplet representation,
## as illustrated by the table DSM_VerbNounTriples_BNC. Such a table may containt additional 
## information, e.g. subcorpus (here: spoken vs. written BNC) or the relation between target
## and feature term (here: noun is the obj[ect] or subj[ect] of the verb).
subset(DSM_VerbNounTriples_BNC, noun == "man" & verb == "see")

## Once you have read such a table into an R data frame, various co-occurrence matrices can be
## computed by filtering the table and/or aggregating frequency counts. For example, to 
## describe verbs by the nouns they take as objects based on data from the spoken BNC:
tbl <- subset(DSM_VerbNounTriples_BNC, rel == "obj" & mode == "spoken")
VObj <- dsm(target=tbl$verb, feature=tbl$noun, score=tbl$f, raw.freq=TRUE, verbose=TRUE)
## Don't forget to specify raw.freq=TRUE so the marginal frequencies and sample size are
## computed automatically from the co-occurrence matrix. By default, dsm() assumes that the
## entries of your sparse matrix are pre-computed scores.

## The triplet representation may contain multiple entries for the same target-feature pair.
## In this case, frequency counts are automatically added up. For example, to describe nouns
## by verbs that take them as objects or subjects in the entire BNC:
NV <- with(DSM_VerbNounTriples_BNC,
           dsm(target=noun, feature=verb, score=f, raw.freq=TRUE, verbose=TRUE))
## The dsm() constructor sums over up to 4 entries for each noun-verb combination.

## Exercise:
##  - How many different DSMs can you compile from DSM_VerbNounTriples_BNC?
##  - Do you remember how to weight DSM features, normalize vectors and find nearest neighbours?


## If you want to load your own data, it is easiest to read the triplet representation from
## a TAB-delimited text file, which may be compressed (.gz, .bz2, ...). As has been pointed out
## in the slides, it is sufficient to provide co-occurrence frequency counts in the case of a 
## syntactic term-term matrix or a term-document matrix; marginal frequencies and sample size
## can be computed automatically by summing over the matrix.

## Loading TAB-delimited files efficiently into a data frame for further filtering is a little
## tricky, so it is often best to prepare a pre-filtered file with just the relevant 3 columns
## externally and load it with read.dsm.triplet(), which does its best to load data efficiently.
## Note that the input file should (i) not contain a header row and (ii) list targets and features
## in the first two TAB-delimited columns, followed by the frequency count (or pre-computed score).
## See ?read.dsm.triplet() for accepted file formats.

## Download the file "verb_dep.txt.gz", which is a typical example of such an input table. Note
## that frequency counts are listed in the first column here, which may be easier to produce with 
## standard tools (e.g. cwb-scan-corpus). We therefore need to specify value.first=TRUE below.
## Also note that we have freq=TRUE instead of raw.freq=TRUE for this function.
VDep <- read.dsm.triplet("verb_dep.txt.gz", freq=TRUE, value.first=TRUE, verbose=TRUE,
                         encoding="UTF-8")
## When working with language data, make sure to always specify the character encoding of the input
## file(s). Note that the option is simply called encoding=, but corresponds to the option known as
## fileEncoding= in read.delim() and similar functions.

## NB: you probably want to discard some of the 285k dimensions before you do anything with this DSM
## (cells with f=1 have already been discarded to make the download faster).

## read.dsm.triplet() expects targets in the first column and features in the second. If you want
## to build the "opposite" model (e.g. noun-verb instead of verb-noun), you can simply transpose
## the co-occurrence matrix:
DepV <- t(VDep)
NdepV <- subset(DepV, grepl("_N$", term)) # just nouns as targets


## For relatively small data sets, it is also possible to read in a list of co-occurrence tokens
## and compute the frequency counts in R.  In this case, the input file has only two columns
## specifying the target and feature term for each co-occurrence token.  You can download the
## file "adj_noun_tokens.txt.gz" from the course homepage as an example of this format:
AdjN <- read.dsm.triplet("adj_noun_tokens.txt.gz", tokens=TRUE, sort=TRUE, verbose=TRUE,
                         encoding="UTF-8")
## Note that freq=TRUE is implied by tokens=TRUE; sort=TRUE ensures that target and feature terms 
## are sorted alphabetically in the matrix representation.

## The same input format can be used for a term-document matrix: each document is tokenized in
## one-word-per-line format, then the tokens are annotated with a document ID in the second column.
## Download the file "delta_de_termdoc.txt.gz" for such a term-document representation of 75 German
## 19th-century novel that have been used for experiments in authorship attribution.
Delta <- read.dsm.triplet("delta_de_termdoc.txt.gz", tokens=TRUE, sort=TRUE, verbose=TRUE,
                          encoding="UTF-8")

## Exercise:
##  - One of the most successful text similarity measures in authorship attribution is Burrows's
##    Delta, which computes Manhattan distances between bag-of-words vectors containing 
##    z-transformed (standardized) relative frequencies of the 500-5000 most frequent words.
##  - Can you use 'wordspace' functions to create the required document-term matrix, determine
##    the most frequent words (MFW) as features, compute standardized relative frequencies, 
##    measure Manhattan distances between texts, and generate a hierarchical clustering?
##  - If not, what are the crucial limitations? Can you find a trick to circumvent them?
##  - As a starting point, lets us annotate the columns with author names:
Delta$cols <- transform(Delta$cols, author=sapply(strsplit(term, ":\\s+", perl=TRUE), `[`, 1))


## Term-document matrices can also be created efficiently with 'tm', the text mining package for R.
library(tm)
data(crude) # news messages on crude oil from Reuters corpus
cat(as.character(crude[[1]]), "\n") # a text example

corpus <- tm_map(crude, stripWhitespace) # some pre-processing
corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeWords, stopwords("english"))
cat(as.character(corpus[[1]]), "\n") # pre-processed text

dtm <- DocumentTermMatrix(corpus) # document-term matrix
inspect(dtm[1:5, 90:99])   # rows = documents

wordspace_dtm <- as.dsm(dtm, verbose=TRUE) # convert to DSM
print(wordspace_dtm$S[1:5, 90:99]) # same part of dtm as above

wordspace_tdm <- t(wordspace_dtm) # convert to term-document matrix
print(wordspace_tdm)


## In the case of surface co-occurrence, marginal frequencies should be provided in separate table files;
## the sample size is usually specified manually when loading the data set.
## As an example, download the file "potter_l2r2.txt.gz", which contains co-occurrence within a span of
## L2/R2 words extracted from a corpus of Harry Potter fan fiction. You will also need a separate file
## "potter_lemmas.txt.gz" with the overall corpus frequencies of all relevant lemmas, which we will use
## as marginals. The sample size (determined during extraction) is N = 63515435 tokens.
## Since the file "potter_lemmas.txt.gz" does not start with a suitable header row, we will have to 
## specify its column labels with the [row,col]info.header option. Note that we use the lemma frequencies
## both as row and as column marginals for our symmetric matrix and set span.size=4 (for L2/R2 span) to 
## obtain correct expected frequencies.
Potter <- read.dsm.triplet("potter_l2r2.txt.gz", freq=TRUE, sort=TRUE, verbose=TRUE, 
                           rowinfo="potter_lemmas.txt.gz", rowinfo.header=c("term", "f"),
                           colinfo="potter_lemmas.txt.gz", colinfo.header=c("term", "f"),
                           N=63515435, encoding="UTF-8")


## If you do not have suitable tools for corpus processing and the extraction of surface co-occurrence data,
## you might be able to use the 'quanteda' package to do everything within R. Unfortunately, its companion
## package 'readtext', which allows you to load corpus texts in many different format, does not seem to be
## available on CRAN yet and has to be installed from https://github.com/kbenoit/readtext
library(quanteda)

## For this example, we will use a small data sample included in the 'quanteda' package itself:
cat(substr(data_char_mobydick, 0, 384))

## A character vector of text samples can easily be converted into a corpus.
Moby <- corpus(data_char_mobydick)
summary(Moby)
kwic(Moby, "necessit*") # matches entire words

## Clean up tokens and apply stemmer (see ?tokens for options)
Moby.tok <- tokens(Moby, removePunct=TRUE, removeNumbers=TRUE)
Moby.tok <- tokens_tolower(Moby.tok)
Moby.tok <- removeFeatures(Moby.tok, stopwords("english"))
Moby.tok <- tokens_wordstem(Moby.tok, language="english")
head(Moby.tok[[1]], 20)

## Surface co-occurrence counts are obtained with the fcm() function, which stands for "feature co-occurrence
## matrix". Be sure to specify tri=FALSE so the full symmetric matrix is returned.
Moby.M <- fcm(Moby.tok, context="window", window=10, tri=FALSE)

## The 'fcm' object is an extension of a sparse matrix; after conversion to a suitable sparseMatrix format,
## it can be passed to the dsm() constructor. Note that this will compute marginal frequencies by summing
## over the rows and columns of the contingency table, leading to inflated sample size but correct E.
MobyDSM <- dsm(as(Moby.M, "sparseMatrix"), verbose=TRUE)
head(MobyDSM, 8)


## Exercise:
##  - Import your own co-occurrence data set (or term-document matrix) into R.
##  - Pick your own corpus, then extract syntactic co-occurrence tokens (or frequencies),
##    a term-document representation, or surface / textual co-occurrence with marginals.
##  - Read the co-occurrence data into R in a suitable way and explore the resulting DSM.
##  - Remember to always specify encoding of the input files, usually encoding="UTF-8".
