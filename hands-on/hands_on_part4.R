######################################################################
## Example code for part 4: DS beyond NLP: Linguistic Theory
##

library(wordspace)
library(wordspaceEval)

######################################################################
## PART 1: sentence vectors
TT <- DSM_TermTermMatrix
TT
# Sentence 1: Cats and dogs need their time (manually lowercase and lemmatized)
s1 <- "cat and dog need their time" 
# Sentence 2: Time is the cause not the effect (manually lowercase and lemmatized)
s2 <- "time is the cause not the effect"

context.vectors(TT, s1)
context.vectors(TT, s2)

# Let's check we are really taking the average: 
(TT['cat', ] + TT['dog', ] + TT['time', ]) / 3
# We are :) 

# Context.vectors() can also take a list as an input. See ?context.vectors for more details
contexts <- context.vectors(TT, c(s1, s2))
contexts
# Let's make rownames more meaningful
rownames(contexts) <- c("s1", "s2")
# And then update the TT by appending our new sentence vectors...
TT <- rbind(TT, contexts)
# And have a look: 
TT

nearest.neighbours(TT, c("s1", "s2"), n=6)

plot(dist.matrix(TT))

######################################################################
# PART 2: compositionality in DSMs

# Let us try to use vector addition and multiplication for adjective/noun composition
# We select DSM_Vectors because lemmas are conveniently annotated by the POS (so, "book_V", "book_N")
load("models/WP500_Win5_Lemma.rda")
DSM <- as.matrix(WP500_Win5_Lemma)

# Vecchi et al. had "parliamentary tomato" as an example of deviant AN

# First we extract the vectors 
parliamentary <- DSM['parliamentary_J', ]
tomato <- DSM['tomato_N', ]

# Then we create the composed representations 
parltom_add <- parliamentary + tomato
parltom_mult <- parliamentary * tomato # * performs element-wise multiplication (Hadamard product) in R

# To compare them to the vectors in the main space, we have two options. 

# Option 1: Either we update DSM_Vectors with the two new vectors, creating a copy of it to leave the main model untouched
DSM_up <- rbind(DSM, parltom_add)
DSM_up <- rbind(DSM_up, parltom_mult)
# ... and calculate nearest neighbours then
nearest.neighbours(DSM_up, c("tomato_N", "parliamentary_J", "parltom_add", "parltom_mult"), n=10)
# ... let's visualize
plot(nearest.neighbours(DSM_up, "parltom_add", n=10, dist.matrix=TRUE))
plot(nearest.neighbours(DSM_up, "parltom_mult", n=10, dist.matrix=TRUE))

# Option 2: creating an extra matrix with the composed vectors
composed_dsm <- dsm(rbind(parltom_add, parltom_mult))
# ... and calculate nearest neighbours using the M2 option of nearest.neighbours()
nearest.neighbours(composed_dsm, c("parltom_add", "parltom_mult"), n=10, M2=DSM)
# Watch out: if you renormalize the composed vectors matrix, you should do the same to DSM_Vectors you use in M2

# Try by yourself with different input models, and different examples of AN (e.g., red {blood, wine, square})
# If you wanted to make the head of the compound count more in the output, which simple trick could you use? 


######################################################################
## BONUS example: disambiguating "bank", based on the two example sentences from the lecture slides
bank_1 <- "The broker went to the bank to secure his cash"
bank_2 <- "The river bank was steep and dangerous"

# And let us load a "real sized" model, word2vec, which is based on inflected words
# and has a large vocabulary. Which make it particularly suitable if we want to do 
# minimal preprocessing
load("models/GoogleNews300_wf200k.rda")
# As usual we give a more "human" name to the model
dsm <- GoogleNews300_wf200k

# Little detour: we want to check the coverage of our sentences in the model 

## Sentence 1
# Let's turn the string into a list:
bank_1_list <- unlist(strsplit(bank_1, " ")) # this perfoms a whitespace tokenization of bank_1
bank_1_list
bank_1_list %in% rownames(dsm)
# how many?
sum(bank_1_list %in% rownames(dsm)) 
# how many, in proportion to length?
sum(bank_1_list %in% rownames(dsm)) / length(bank_1_list)
# which ones?
bank_1_list[bank_1_list %in% rownames(dsm)]

## Sentence 2 (code variation)
bank_2_list <- unlist(strsplit(bank_2, " "))
sum(bank_2_list %in% rownames(dsm)) / length(bank_2_list)
# Negate the condition with ! to see directly which words are missing
bank_2_list[!(bank_2_list %in% rownames(dsm))]

# And now we create a matrix with the context vectors of our two sentences, using the 
# word2vec representations as our input
bank_candidates <- context.vectors(dsm, c(bank_1, bank_2))
rownames(bank_candidates) <- c("bank1", "bank2")

# We can use the M2 option of the nearest neighbor function ?nearest.neighbours
# to find the nearest neighbours of the selected row in the full embeddings matrix 'dsm'
nearest.neighbours(bank_candidates, "bank1", n=10, M2=dsm)
nearest.neighbours(bank_candidates, "bank2", n=10, M2=dsm)

# We can get the "usual" semantic map for the neighbors of "bank" in word2vec
plot(nearest.neighbours(dsm, "bank", n= 10, dist.matrix=TRUE))
# ... and also the semantic maps for our sentences, "bank1" and "bank2", 
# extracting their nearest neighbours from word2vec
plot(nearest.neighbours(bank_candidates, "bank1", n=10, M2=dsm, dist.matrix=TRUE))
plot(nearest.neighbours(bank_candidates, "bank2", n=10, M2=dsm, dist.matrix=TRUE))


######################################################################
## ONE DATASET, MANY TASKS: LAZARIDOU 2013

# We talked about the Lazaridou2013 dataset in the lecture. Conveniently, it is part of wordspaceEval
# Have a look at it:
View(Lazaridou2013)
# Read the documentation:
?Lazaridou2013
# And finally, get an idea of its size/scope with the summary function
summary(Lazaridou2013)

# This dataset hasn't been developed for standard NLP tasks. Yet, we can still construct such 
# evaluation tasks from it and solve them using the 3 functions we learnt today.

## STEP 1: multiple choice ##

# To carry out multiple choice on Lazaridou2013, we need to create the distractors 
# (remember? TOEFL has one word, one candidate synonym, and 3 distractor words)
# Here we have pairs of stems and derived words 
# (stem: drive, derived: driver; stem: happy, derived: unhappy)
# Base can be our input word, and the task for the model will be to find the derived word. 
# How can we generate distractors? Idea: given a pair of base/derived words, 
# we can take 3 random other derived words as distractors. Then the assumption will be:
# If my DSM is "good", it will be able to spot the relatedness between drive/driver, and assign
# it a higher similarity (or lower distance) than the other derived words. 
# Of course it is not a perfect strategy (can you think about potential issues with it?)
# But let's get started nevetherless :) 

# We set a random seed to be able to replicate the results (this is customary
# when doing random shuffling
set.seed(001)
candidates1 <- sample(Lazaridou2013$derived)
candidates1
candidates2 <- sample(Lazaridou2013$derived)
candidates2
candidates3 <- sample(Lazaridou2013$derived)
candidates3

# Have a look at the ?eval.multiple.choice function: it actually allows 
# you to specify the name of the target (for us it is "stem") 
# of correct choice (for us it is "derived") the column for targets and the star
# of the string for the column names specifying distractors. The default is distract*, 
# which is good enough for us. We will have add 3 columns to Lazaridou2013:
# distractor1, distractor2, distractor3

# Since we're going to manipulate the data set, let us create a copy. (Footnote: This isn't
# strictly necessary. You can add to or otherwise modify a read-only data set loaded from 
# a package. R will automatically create a copy in your workspace under the same name, which 
# "shadows" the version from the wordspaceEval package. But this can become quite confusing,
# so it's better to make our own copy to begin with.)
MyLazaridou2013 <- Lazaridou2013
MyLazaridou2013$distractor1 <- candidates1
MyLazaridou2013$distractor2 <- candidates2
MyLazaridou2013$distractor3 <- candidates3

# Let us experiment with word2vec, as it is the largest we have and state-of-the-art.
# We will lose a bit in performance because w2v contains inflected words (so, two vectors for
# dog and dogs, and we will be only using the singular (stem) in this case. 
# But we will gain in coverage so we are good for now. Have fun experimenting with other models;
# watch out for coverage though! 
load("models/GoogleNews300_wf200k.rda")
dsm <- GoogleNews300_wf200k

# And finally we can call multiple choice, making sure we indicate the right names
eval.multiple.choice(MyLazaridou2013, dsm, target.name="stem", correct.name="derived", distractor.name="distract")
# not bad... Further check: is the performance comparable in test vs. training?
# but can we use details=TRUE to investigate things a bit deeper?
Lazaridou2013_multiple <- eval.multiple.choice(MyLazaridou2013, dsm, target.name="stem", correct.name="derived", distractor.name="distract", details=TRUE)

# First of all, we can now have a look at the mistakes made by our model (some are missing items, marked with "Inf")
head(subset(Lazaridou2013_multiple, !correct)) # !correct matches items where the Boolean variable correct is FALSE

# But what is also interesting is the "correct.dist" field, which is actually telling us the distance 
# between base and derived word. Now, wouldn't it be nice to know whether some affixes have smaller/bigger 
# distances? We don't have this info here directly, but we can retrieve it from the main dataframe
# because details=TRUE obviously has kept the order! 
# Let's also keep track of the part of speech of the derived word, and of the stem

Lazaridou2013_multiple$affix <- Lazaridou2013$affix
Lazaridou2013_multiple$derivedPOS <- Lazaridou2013$derivedPOS
Lazaridou2013_multiple$stemPOS <- Lazaridou2013$stemPOS

head(Lazaridou2013_multiple)

# One more step: let us get rid of the Inf cases (remember that from 
# now on we cannot concatenate columns anymore, because we have altered the order)
nrow(Lazaridou2013_multiple)
Lazaridou2013_multiple <- subset(Lazaridou2013_multiple, correct.dist != Inf)
nrow(Lazaridou2013_multiple)

# Now we have a lot of information to play with, and we can go wild ;) 
# For example...
boxplot(correct.dist ~ derivedPOS, data=Lazaridou2013_multiple)
boxplot(correct.dist ~ paste(derivedPOS,stemPOS), data=Lazaridou2013_multiple)
boxplot(correct.dist ~ affix,data=Lazaridou2013_multiple)

## STEP 2: prediction of similarity ratings and of quality ratings. ##
# Remember: annotations are available only for the test set

# Base/derived relatedness score
eval.similarity.correlation(subset(Lazaridou2013, set=="test"), dsm, word1.name="stem", word2.name="derived", score.name="relatedness_score")
# Not so good, and with these vectors we have little room for parameter manipulation. How about using rank? (Careful: This is computationally very heavy!)
eval.similarity.correlation(subset(Lazaridou2013, set=="test"), dsm, word1.name="stem", word2.name="derived", score.name="relatedness_score", rank="fwd")
####                                            rho      p.value missing         r     r.lower    r.upper
#### subset(Lazaridou2013, set == "test") 0.2726184 8.405727e-17      65 0.0153527 -0.05004555 0.08061985
eval.similarity.correlation(subset(Lazaridou2013, set=="test"), dsm, word1.name="stem", word2.name="derived", score.name="relatedness_score", rank="avg")
####                                            rho      p.value missing          r     r.lower    r.upper
#### subset(Lazaridou2013, set == "test") 0.2821681 6.176956e-18      65 0.02236349 -0.04304745 0.08758353

# Quality score of the derived vectors (which makes little sense for us, given that it is the quality of the vectors the authors 
# used in their original study). 
eval.similarity.correlation(subset(Lazaridou2013, set=="test"), dsm, word1.name="stem", word2.name="derived", score.name="quality_score")

## STEP 3: categorization ##

# We can use the eval.clustering function to test how well affixes group together...
eval.clustering(Lazaridou2013, dsm, word.name="derived", class.name="affix")
# Or if we can categorize derived word in their part of speech... (less classes, gets easier)
eval.clustering(Lazaridou2013, dsm, word.name="derived", class.name="derivedPOS")
# Or if we can categorize derived word in their part of speech... (gets slightly easier, because stems tend to be more frequent and thus 
# have better representations)
eval.clustering(Lazaridou2013, dsm, word.name="stem", class.name="stemPOS")
