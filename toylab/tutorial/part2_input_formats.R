##
## Practice session for part 2:
## How to read your own co-occurrence data into 'wordspace'
##

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


################################################################################
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
## See ?read.dsm.triplet for accepted file formats.

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
## (cells with f=1 have already been dropped to make the download faster).

## read.dsm.triplet() expects targets in the first column and features in the second. If you want
## to build the "opposite" model (e.g. noun-verb instead of verb-noun), you can simply transpose
## the co-occurrence matrix:
DepV <- t(VDep)
NdepV <- subset(DepV, grepl("_N$", term)) # just nouns as targets
## -> this is a good moment to read up on ?subset.dsm

## For relatively small data sets, it is also possible to read in a list of co-occurrence tokens
## and compute the frequency counts in R.  In this case, the input file has only two columns
## specifying the target and feature term for each co-occurrence token.  You can download the
## file "adj_noun_tokens.txt.gz" from the course homepage as an example of this format:
AdjN <- read.dsm.triplet("adj_noun_tokens.txt.gz", tokens=TRUE, sort=TRUE, verbose=TRUE,
                         encoding="UTF-8")
## Note that freq=TRUE is implied by tokens=TRUE; sort=TRUE ensures that target and feature terms 
## are sorted alphabetically in the matrix representation.


################################################################################
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

## Don't peek now, but a solution to the exercise can be found at the end of the script file.


################################################################################
## Term-document matrices can also be created efficiently with 'tm', the text mining package for R.
library(tm)
data(crude) # news messages on crude oil from Reuters corpus
cat(as.character(crude[[1]]), "\n") # a text example

crude.corp <- tm_map(crude, stripWhitespace) # some pre-processing
crude.corp <- tm_map(crude.corp, content_transformer(tolower))
crude.corp <- tm_map(crude.corp, removePunctuation)
crude.corp <- tm_map(crude.corp, removeWords, stopwords("english"))
cat(as.character(crude.corp[[1]]), "\n") # pre-processed text

dtm <- DocumentTermMatrix(crude.corp) # document-term matrix
inspect(dtm[1:5, 90:99])              # rows = documents

wordspace_dtm <- as.dsm(dtm, verbose=TRUE) # convert to DSM
print(wordspace_dtm$S[1:5, 90:99]) # same part of dtm as above

wordspace_tdm <- t(wordspace_dtm) # convert to term-document matrix
print(wordspace_tdm)


################################################################################
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
                           span.size=4, N=63515435, encoding="UTF-8")


################################################################################
## If you do not have suitable tools for corpus processing and the extraction of surface co-occurrence data,
## you might be able to use the 'quanteda' package to do everything within R. Its companion package
## 'readtext' makes it easy to load corpus texts in many different formats.
library(quanteda)

## For this example, we will use a small data sample included in the 'quanteda' package itself:
crude.eda <- corpus(crude.corp)      # convert TM corpus from above into quanteda format
texts(crude.eda)[[1]]                # the first text in the corpus
kwic(crude.eda, "economy", window=3) # concordance search ("keyword in context")
## -> see quanteda documentation and tutorials for more information: http://quanteda.io/ 

## Split text on whitespace and apply stemmer
crude.tok <- tokens(crude.eda, "fasterword")
crude.tok <- tokens_wordstem(crude.tok, language="english")
crude.tok[[1]]

## Surface co-occurrence counts are obtained with the fcm() function, which stands for "feature co-occurrence
## matrix". Be sure to specify tri=FALSE so the full symmetric matrix is returned.
crude.M <- fcm(crude.tok, context="window", window=10, tri=FALSE)

## fcm() also computes marginal frequencies, but they may be incompatible with the co-oc matrix
## (unless the corpus has been converted to lowercase), so we need to obtain the marginals separately
crude.tf <- colSums(dfm(crude.tok, tolower=FALSE))
crude.cols <- data.frame(term=names(crude.tf), f=crude.tf) # marginal frequency table
crude.rows <- transform(crude.cols, f = 20 * f) # adjust row marginals for span size
crude.N <- sum(crude.tf) # sample size

## The 'fcm' object is an extension of a sparse matrix; after conversion to a suitable sparseMatrix format,
## it can be passed to the dsm() constructor together with the marginal frequencies. 
crudeDSM <- dsm(as(crude.M, "sparseMatrix"), verbose=TRUE,
                rowinfo=crude.rows, colinfo=crude.cols, N=crude.N)
crudeDSM
head(crudeDSM, 15)

## Corpus preprocessing can also be carried out with "quanteda" rather than "tm"
Immi <- data_char_ukimmig2010 # texts from UK party manifestos
cat(substr(Immi[[2]], 1, 710))

Immi.tok <- tokens(Immi, remove_punct=TRUE)
Immi.tok <- tokens_tolower(Immi.tok)
Immi.tok <- tokens_remove(Immi.tok, stopwords("english"))
Immi.tok <- tokens_wordstem(Immi.tok, language="english")
head(Immi.tok[[2]], 50)

## Practice: build a wordspace DSM from Immi.tok (as shown above), using a L5/R5 surface span


## Further packages you might want to take into consideration:
##  - corpustools (for working with tokenized corpus data)
##  - polmineR    (analysis of CWB-indexed corpora, see cwb.sf.net)
## 
## NLP annotation pipelines (parameter files for various languages can be downloaded):
##  - spacyr
##  - udpipe
##  - openNLP
##  - coreNLP


################################################################################
## udpipe is relatively easy to install and use, so let's use it to procsss a 
## (very) small example corpus available from the tutorial homepage. Make sure that
## you have downloaded the file "VSS.txt" and placed it in your working directory.
library(udpipe)

## First we need to obtain a suitable udpipe model. Pre-compiled models are available
## for 61 different languages, often with several models for the same language, see
## https://github.com/jwijffels/udpipe.models.ud.2.3/blob/master/inst/udpipe-ud-2.3-181115/README
?udpipe_download_model # help page also lists known models

## Download model to your working directory. You can skip this step if you already have
## done so before, but overwrite=FALSE also ensures that the existing file is used.
ud_en_file <- udpipe_download_model("english-ewt", overwrite=FALSE)$file_model

## Load model file into memory, returning an object that can be passed to annotation functions
ud_en <- udpipe_load_model(ud_en_file)

## Our data set is a very small corpus of six short stories, provided as a TAB-delimited table
## with text ID in the first column and unformatted raw text in the second. This is a format 
## used by several NLP packages in R to represent a document collection.
## For real-life data sets, it is better to use more efficient table loaders from the "readr"
## or "iotools" packages, or the convenient fread() from "data.table".
vss_txt <- read.delim("VSS.txt", header=FALSE, quote="", col.names=c("doc_id", "text"), 
                      fileEncoding="UTF-8", stringsAsFactors=FALSE)

## udpipe processes such data frames, or you can pass the texts simply as a character vector
vss <- udpipe(vss_txt, object=ud_en, trace=TRUE) # note the model is called "object" here
vss <- vss[, !(colnames(vss) == "sentence")]     # remove full sentences from table for readability

## udpipe returns annotated text in "vertical" format, i.e. as a table with rows = tokens
head(vss, 16)

## NB: You can also annotate text with an external software package that generates CoNLL-U
## format and then load the corpu with udpipe_read_conllu().


## In order to obtain proper frequency counts (especially for marginals), we cannot rely on
## functions from the "udpipe" package but have to carry out frequency aggregations and
## database-like joins ourselves. The "data.table" package (which is also used internally by
## "udpipe") provides an efficient and powerful implementation of such operations. It is one
## of the most useful things you can learn about data processing in R.
## First, we convert the data frame back into a data.table object:
library(data.table)
setDT(vss)              # convert data.frame to data.table


## Creating a term-document matrix is easy from the original tokens (setting score=1 for all target-feature pairs).
## You could also compute by-text frequency counts with data.table and pass them on to the dsm() constructor.
## Try both approaches on a large corpus to see which is more efficient.
vss.tdm <- dsm(target=vss$lemma, feature=vss$doc_id, score=1, raw.freq=TRUE, verbose=TRUE)

## If the texts had paragraph breaks, we could also create a term-paragraph matrix by constructing a unique 
## identifier for each paragraph from vss$doc_id and vss$paragraph_id.


## The cooccurrence() function helps us to compute textual and surface co-occurrence counts.  However, correct
## application can be a bit tricky and often requires pre- or post-processing.

## Co-occurrence within sentence units: The cooccurrence() function counts all possible combinations
## between multiple occurrences in the same sentence, potentially leading to substantially inflated counts.
## We therefore remove duplicates within sentences, using an aggregation (GROUP BY) in data.table.
tmp <- vss[, .(f=1), by=.(doc_id, sentence_id, lemma)]

## From this table, we can compute both within-sentence cooccurrence and the correct marginal frequencies.
tmp.margin <- tmp[, .(f=.N), by=.(lemma)] # marginal frequency = number of sentences containing lemma
setnames(tmp.margin, "lemma", "term")     # format of row/column information expected by dsm()
head(tmp.margin)

tmp.cooc <- cooccurrence(tmp, group=c("doc_id", "sentence_id"), term="lemma") # cooccurrence counts within each "group" of tokens
setDT(tmp.cooc)         # make sure the data.table indexing below works
tmp.cooc

## The sample size N is the number of textual units, i.e. unique (doc_id, sentence_id) combinations
tmp.N <- uniqueN(tmp, by=c("doc_id", "sentence_id"))  # uniqueN = number of unique entries

## tmp.cooc counts all possible combinations, but always returns alphabetically ordered pairs. We therefore
## have to add the mirrored pairs (with same frequency count) to obtain the full symmetric cooccurrence matrix.
tmp.cooc <- rbind(tmp.cooc, 
                  tmp.cooc[, .(term1=term2, term2=term1, cooc)])

## Now we can construct the term-term matrix from the co-occurrence data in tmp.cooc and the marginal
## frequencies in tmp.margin; this step should look familiar.
vss.sent <- dsm(target=tmp.cooc$term1, feature=tmp.cooc$term2, score=tmp.cooc$cooc,
                rowinfo=tmp.margin, colinfo=tmp.margin, N=tmp.N,
                raw.freq=TRUE, sort=TRUE, verbose=TRUE)
head(vss.sent, 12) # NB: matrix won't be symmetric unless sort=TRUE!


## Surface co-occurrence can also be obtained from cooccurrence(), but the function only searches
## a one-sided span to the right of the target. Again, we have to mirror the co-occurrence data 
## for the full symmetric matrix. Can you work out how to compute co-occurrence frequencies for
## an asymmetric span, e.g. L2/R5 (2 tokens to the left of the target, 5 tokens to the right)?

## Note that we have to apply cooccurrence() to a single column rather than the entire dataframe
## in this case -- it's a weird way to distinguish textual and surface co-occurrence, isn't it?
##  - skpigram=4 selects a L0/R5 span (i.e. up to 4 tokens can be "skipped" to find a cooccurrence)
##  - relevant=(condition) allows us to filter tokens, e.g. to only accept content words
lexical.pos <- c("NOUN", "VERB", "ADJ", "ADV")
tmp.cooc <- cooccurrence(vss$lemma, skipgram=4,
                         relevant=(vss$upos %in% lexical.pos))
setDT(tmp.cooc)
tmp.cooc

tmp.cooc <- rbind(tmp.cooc,
                  tmp.cooc[, .(term1=term2, term2=term1, cooc)])

## Marginal frequencies are simply the individual corpus frequencies of the lemmas, but we need 
## to multiply the row marginals by the span size adjustment, which is 2 x 5 = 10 in this case
## (dsm() doesn't have a span.size= argument like read.dsm.triplet()).
tmp.rows <- tmp[, .(f=10 * .N), by=.(lemma)]
setnames(tmp.rows, "lemma", "term") # don't forget to rename "lemma" to "term"
tmp.cols <- tmp[, .(f=.N), by=.(lemma)]
setnames(tmp.cols, "lemma", "term")
tmp.N <- nrow(tmp) # sample size = corpus size in tokens

## We're ready to construct the distributional model now.
vss.win5 <- dsm(target=tmp.cooc$term1, feature=tmp.cooc$term2, score=tmp.cooc$cooc,
                rowinfo=tmp.rows, colinfo=tmp.cols, N=tmp.N,
                raw.freq=TRUE, verbose=TRUE)
head(vss.win5, 15) # allow non-symmetric matrix so the top left corner looks more interesting


## Syntactic co-occurrence: udpipe also produces simple dependency parses and it is easy to 
## create a (filtered or structured) distributional model based on direct dependency links.
## However, due to the basic dependency scheme, many interesting relations require multiple
## steps through the dependency graph, which is quite tricky to implement. We will therefore
## focus on direct dependencies in this example.

## Fortunately, udpipe has a convenience function to annotate the lemma, pos, ... of dependency
## parents together with the child token. This makes it very easy to filter by part-of-speech
## and extract lemmatized co-occurrence tokens (i.e. instances of a dependency relation).
tmp <- cbind_dependencies(vss)
tmp <- subset(tmp, upos %in% lexical.pos & upos_parent %in% lexical.pos)

## We can now directly obtain the co-occurrence tokens and let dsm() compute marginal 
## frequencies by summing over rows and columns.  Keep in mind that we have to swap 
## parent and child roles if we want both child->parent and parent<-child dependencies.
vss.dep <- dsm(target=c(tmp$lemma, tmp$lemma_parent), feature=c(tmp$lemma_parent, tmp$lemma),
               score=1, raw.freq=TRUE, sort=FALSE, verbose=TRUE)
head(vss.dep, 15) # again not symmetric so top-left corner shows something interesting


################################################################################
## Using coreNLP or openNLP requires the rJava interface, which can be a bit painful to get to work.
## If you manage to install everything, you can use coreNLP to obtain enhanced dependencies, which
## capture most relations of interest as direct links.
library(coreNLP)

## First, you need to download suitable CoreNLP models. Here we only take the base package for English
## (which is always required). Several other languages are also available, see ?downloadCoreNLP.
## The following line needs to be run only once -- the model will be stored inside the coreNLP package.
downloadCoreNLP(type="base")  # almost 400 MiB, will take quite long!!

## At the start of each CoreNLP session, you have to initialize the desired model.  On my computer,
## the default English model (and even more so "english_all") fails to process even our small
## 8000-word sample -- it hangs for more than an hour with high CPU load.  The minimal model
## "english_fast" works, but does not include syntactic analysis.  You will need to edit its 
## parameter file, adding "parse" (and possibly "ner") to the list of annotation modules.
edit(file=file.path(system.file("extdata", package="coreNLP"), "StanfordCoreNLP-english-fast.properties"))
## In RStudio, this should open a file edit dialog window. Click "Save" after you have made 
## the required changes.  The first line of the parameter file should now look like this
##   annotators = tokenize, ssplit, pos, lemma, ner, parse

## Initialize a "fast" English model that only runs the required components:
initCoreNLP(type="english_fast")

## Now you can simply use annotateString() to process plain text. Note that document structure
## is not preserved by the function. If you need the information, you have to process each text
## separately in a loop.
vss.corenlp <- annotateString(vss_txt$text) # about 1 minute on a fast laptop

## The main advantage of CoreNLP are the enhanced dependency relations, which are
## extracted by default with the getDependency() function.
vss.deprel <- getDependency(vss.corenlp)
vss.deprel[26:30, ]

## vss.deprel only lists the word forms of head and dependent. However, the govIndex
## and depIndex columns refer to rows in the token annotation table, from which lemma
## and other linguistic attributes can be extracted.
vss.tokens <- getToken(vss.corenlp)
vss.tokens[25:30, ]  # compare index columns above with this table

## We will extract unfiltered dependency relations here but restrict the head and dependent
## to nouns, verbs and adjectives. Before we look up the extra information, remove items
## with invalid govIndex or depIndex (marking the root node of a sentence).
vss.deprel <- subset(vss.deprel, !is.na(govIndex) & !is.na(depIndex))

## Add lemma and POS for head and dependent based on the govIndex and depIndex offsets into vss.tokens.
vss.deprel <- transform(vss.deprel,
                        head_lemma = vss.tokens$lemma[govIndex],
                        head_pos = vss.tokens$POS[govIndex],
                        dep_lemma = vss.tokens$lemma[depIndex],
                        dep_pos = vss.tokens$POS[depIndex],
                        stringsAsFactors=FALSE)

## Now reduce the table to the desired parts of speech (nouns, verbs, adjectives).
vss.deprel <- subset(vss.deprel, grepl("^[NVJ]", head_pos) & grepl("^[NVJ]", dep_pos))

## We are now ready to construct the DSM object, keeping in mind that we need to add
## the mirror pairs to include both "upward" and "downward" dependencies.
vss.coredep <- dsm(target=c(vss.deprel$head_lemma, vss.deprel$dep_lemma),
                   feature=c(vss.deprel$dep_lemma, vss.deprel$head_lemma),
                   score=1, raw.freq=TRUE, sort=TRUE, verbose=TRUE)

subset(vss.coredep, f >= 30, f >= 45, matrix.only=TRUE) # most frequent targets / features

## For large data sets, you will want to perform the operations more efficiently using only
## the required columns from the data frames, and with direct vector indexing rather than
## subset(); alternatively, perform the merge and frequency counts efficiently with data.table.


################################################################################
## Main exercise:
##  - Import your own co-occurrence data set (or term-document matrix) into R.
##  - Pick your own corpus, then extract syntactic co-occurrence tokens (or frequencies),
##    a term-document representation, or surface / textual co-occurrence with marginals.
##  - Read the co-occurrence data into R in a suitable way and explore the resulting DSM.
##  - Remember to always specify encoding of the input files, usually encoding="UTF-8".


## Further exercises:
##  - modify the udpipe example to compute asymmetric spans (e.g. L2/R5)
##  - modify the udpipe example to obtain _structured_ surface and syntactic context,
##    (i.e. distinguish between left/right co-occurrence and different dependency relations)
##  - modify the coreNLP example to focus on specific dependency relations: 
##    what is a sensible set of relations, depending on the part of speech of targets?


## Solution to the "Delta" exercise from above:
##  - the object Delta contains a term-document matrix, which needs to be transposed
##    to obtain the document-term matrix (where targets are the individual novels)
Delta <- t(Delta) 
##  - we want to scale the term frequency counts to relative frequencies, i.e. 
##    divide them by the length of the novel in tokens (equation: O / R1)
##  - this can be achieved with a user-defined association measure as follows
Delta <- dsm.score(Delta, score=function(O, R1, ...) O / R1)
##  - alternatively, notice that the relative frequencies of all terms in a document
##    must add up to 1.0 (because they are fractions of the tokens in the document)
##  - normalizing vectors to have a Manhattan length of 1.0 ensures this is the case;
##    it has exactly the same effect as the user-defined score above
Delta <- dsm.score(Delta, score="frequency", normalize=TRUE, method="manhattan")
##  - we will apply Delta here based on the 200 most frequent words (MFW); since we
##    do not know the corresponding frequency threshold, we filter the columns of 
##    Delta by ranks (note that rank(-f) ranks the highest frequencies first)
Delta <- subset(Delta, select=(rank(-f) <= 200))
##  - following Burrows (2002), the relative frequencies are converted to z-scores,
##    i.e. the distribution of values in each matrix _column_ is standardized
##  - we want to apply this column scaling to the realtive frequencies in Delta$S,
##    so we specify score="reweight"
##  - negative.ok=TRUE is required to enfoce creating a dense matrix 
Delta <- dsm.score(Delta, score="reweight", scale="standardize", negative.ok=TRUE)
##  - now compute the distance matrix between all 75 novels, again following 
##    Burrows (2002) to use the Manhattan metric (rather than better "cosine")
Delta.dm <- dist.matrix(Delta, method="manhattan")
##  - we can simply plot the distance matrix for a 2-d visualisation
plot(Delta.dm, labels=Delta$rows$author)
##  - a hierarchical clustering shows that Delta groups texts from the same author
##    fairly well; can you spot its mistakes?
plot(hclust(as.dist(Delta.dm)))
