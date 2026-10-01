##
## Hands-on for day 5: The FAST dataset and tasks
##

library(wordspace)
library(wordspaceEval)

## Take a first look at the FAST data set
head(FAST)
tail(FAST)
?FAST # and always read the documentation!

## Number of task items (stimulus, first, hapax, random) in EAT and USF subsets
table(FAST$norm, ifelse(FAST$in_test, "test", "training"))

## We could evaluate FAST as a single data set, but there are substantial differences
## between the two association norms (cultural differences between BrE and AmE, sampled
## in different time periods). In fact, many stimuli occur both in the EAT and USF
## subsets, but with different FIRST associates.
## Here, let us evaluate the EAT subset. You can practice by carrying out the same
## evaluation for USF.
EAT <- subset(FAST, norm == "EAT")


######################################################################
## Multiple choice with DSM

## We will focus on using DSM similarity for the multiple-choice task, using one of
## the pre-compiled models with SVD dimensionality reduction. You have probably
## already noticed in the hands-on sessions that nearest.neighbours() and similar functions
## are much more efficient in a dense SVD space than for the unreduced models.

## Feel free to pick another one of the pre-compiled models, of course, e.g. neural embeddings.
load("models/WP500_Win5_Lemma_svd500.rda", verbose=TRUE)
DSM <- WP500_Win5_Lemma_svd500
dim(DSM) # note that the reduced model is a raw matrix, not a 'dsm' object

## We don't need training data for this task, so let us extract the test set.
Test <- subset(EAT, in_test)

## We can now simply use the eval.multiple.choice() function, we just need to specify the
## appropriate column names for target, correct choice and distractors. Pick the lemma
## versions because our DSM uses POS-disambiguated lemmas as targets.
## NB: If you loaded a neural embedding, you will need to adapt the command below so that
##     it selects the word forms of stimulus and responses.
result <- eval.multiple.choice(Test, DSM,
            target.name="stimulus.lemma",               # target = stimulus
            correct.name="FIRST.lemma",                 # correct choice = FIRST associate
            distractor.name="^(HAPAX|RANDOM)\\.lemma")  # distractors: HAPAX and RANDOM (matched with regular expression)
result

## There are many test items for which our DSM doesn't provide enough information (if either
## the stimulus or the FIRST response are not among DSM targets, it's impossible to get the
## correct answer), shown as 'missing' in the result table. These are counted as false positives,
## leading to poor accuracy for models with low vocabulary coverage.
## While this is a realistic measure of performance, we may still want to know how our model
## performed on the items covered by its vocabulary. We can use transform() to compute this
## adjusted accuracy:
transform(result, acc.adj=100 * TP / (TP + FP - missing)) # this looks a lot better!

## You can now experiment with different distance measures or with using neighbour rank instead
## of raw distance values. The evaluation code uses pair.distances() internally, and you can
## pass any of the parameters accepted by pair.distances().
## E.g. evaluation on average neighbour rank (this will take much longer to execute!)
result <- eval.multiple.choice(Test, DSM,
                               rank="avg", avg.method="geometric",
                               target.name="stimulus.lemma",
                               correct.name="FIRST.lemma",
                               distractor.name="^(HAPAX|RANDOM)\\.lemma")
transform(result, acc.adj=100 * TP / (TP + FP - missing))


######################################################################
## Multiple choice with first-order (FO) collocations

## Where can we obtain collocations, i.e. first-order co-occurrence data?
## Recall that a DSM has at its heart a co-occurrence matrix, in our case for a L5/R5 surface span.
## dsm.score() can be used to compute various association measures that are commonly used in 
## collocation research. In other words: a scored DSM contains collocational profiles, i.e. 
## first-order data. It only becomes second-order when we compute distances between rows.

## Let us now load the co-occurrrence matrix for a L5/R5 surface span.
load("models/WP500_Win5_Lemma.rda", verbose=TRUE)
FO <- WP500_Win5_Lemma # FO = first order
FO # this is a 'dsm' object, of course

## Apply dsm.score() to compute association scores as a measure of collocational strength.
## The evaluation showed optimal results for MI2, but dsm.score() doesn't directly support
## this measure. Let us go with simple log-likelihood instead (NB: local-MI isn't provided
## because it's practically equivalent to simple-ll, and the latter has a more solid
## mathematical foundation).
FO <- dsm.score(FO, score="simple-ll") # always compute sparse association scores!

## If you really want to use MI2, you'll have to specify it as a user-defined equation:
## FO <- dsm.score(FO, score=(function (O, E, ...) O^2 / E), transform="log")

## The last challenge is how to trick pair.distances() into using the association scores
## as a similarity measure rather than computing distances between rows of the matrix.
## We do so by marking it as a pre-computed similarity matrix.
FO.mat <- as.distmat(FO, similarity=TRUE)

## Now evaluate the FO data on the EAT test set. Specifying a distance measure doesn't make
## sense, of course, but you can try neighbour rank.
## NB: If you pass FO instead of FO.mat, the evaluation code will compute distances between
##     rows in the unreduced DSM ... this is not what you want!
result <- eval.multiple.choice(Test, FO.mat,
                               target.name="stimulus.lemma",
                               correct.name="FIRST.lemma",
                               distractor.name="^(HAPAX|RANDOM)\\.lemma")
transform(result, acc.adj=100 * TP / (TP + FP - missing))


######################################################################
## Combining second-order (DSM) and first-order (FO) information

## Since distance metrics and assocation scores are on entirely different scales
## (and the same holds between different metrics and between different AMs),
## the only sensible way of combining them is at the level of neighbour rank.
## As long as neighbour ranks are based on the same vocabulary (which we ensure by
## using the same underlying co-occurrence matrix in both cases), they are directly
## comparable and it makes sense to average the DSM and FO neighbour rank.

## The standard evaluation functions and pair.distances() don't support such complex
## combinations, so you will have to carry out the entire evaluation yourself now.
## (You can actually smuggle the combination into eval.multiple.choice(), but this
## involves more advanced R programming and we will not attempt it here.)

## Let us compute all relevant neighbour ranks. pair.distances() returns rank Inf (infinity)
## if a pair does not occur in the DSM or FO matrix.
rank.first.DSM <- pair.distances(Test$stimulus.lemma, Test$FIRST.lemma,  DSM, rank="fwd")
rank.hapax.DSM <- pair.distances(Test$stimulus.lemma, Test$HAPAX.lemma,  DSM, rank="fwd")
rank.randm.DSM <- pair.distances(Test$stimulus.lemma, Test$RANDOM.lemma, DSM, rank="fwd")

rank.first.FO <- pair.distances(Test$stimulus.lemma, Test$FIRST.lemma,  FO, rank="fwd")
rank.hapax.FO <- pair.distances(Test$stimulus.lemma, Test$HAPAX.lemma,  FO, rank="fwd")
rank.randm.FO <- pair.distances(Test$stimulus.lemma, Test$RANDOM.lemma, FO, rank="fwd")

## Compute the average ranks. We will take the geometric mean to favour low ranks in the
## combination, but you could also try the arithmetic or harmonic mean - or even the
## minimum of the two ranks (what would be the intuition behind this?).
rank.first <- sqrt(rank.first.DSM * rank.first.FO)
rank.hapax <- sqrt(rank.hapax.DSM * rank.hapax.FO)
rank.randm <- sqrt(rank.randm.DSM * rank.randm.FO)

head(rank.first, 10) # take a quick look at the first 10 items
head(rank.hapax, 10)
head(rank.randm, 10)

## The combined model is correct if the FIRST response has a lower neighbour rank than
## either of the distractors. Note the use of & (instead of &&) for a vectorized AND.
is.correct <- (rank.first < rank.hapax) & (rank.first < rank.randm)
table(is.correct)      # number of TPs and FPs
100 * mean(is.correct) # and the accuracy of the model

## Your turn now: Identify the "missing" items (hint: look for combined rank = Inf) to
## compute adjusted accuracy. Then try different models and rank combinations.


######################################################################
## The open-vocabulary lexical access task

## Implementation of this task is left as an exercise for the enthusiastic reader. :o)

