## Create a separate directory for the course, into which you will download
## all the materials (examples scripts, data sets and pre-compiled DSMs),
## and which you will also use for your own R scripts during the exercises.

## If you work in RStudio, it is best to create a new project based on this
## directory ("New Project" -> "Existing Directory"). RStudio will then make
## sure that the working directory of your R session is set correctly.
## You can easily switch between different projects from the drop-down menu
## in the top right corner of the RStudio window.

## If you do now want to use an RStudio project (e.g. because you prefer a
## different IDE such as Emacs with the excellent ESS interface), you will
## need to set your working directory at the start of every R session.
## You should include this in any R script you write (and will need to add
## it to the example scripts downloaded from the wiki), so that the scripts
## can run as standalone programs. In short: use RStudio projects!

## For example, Gabriella doesn't like RStudio projects, so she has to set
## the working directory manually. She has created a folder on her desktop
## called "loveDSMs" for the course materials. At the start of every session,
## she has to execute the command below to set her working directory -- you
## will need to insert the path to your course folder, of course.
setwd("~gabriella/Desktop/loveDSMs")

## We will organize data files into a few subdirectories in order to keep
## your course folder clean and structured. Please create subdirectories
##
##   models  ...  for pre-compiled DSMs
##   data    ...  for other data files
##
## Our example scripts will always assume that your data files are
## organized into these directories. Place all R scripts (examples
## downloaded from the course wiki as well as your own) in the main
## course directory.

## The hands-on examples and exercises in this course rely on "wordspace",
## a user-friendly R package for distributional semantics. We will also
## use some additional non-public data sets in the "wordspaceEval" package,
## which you can download from the course wiki and install locally, as 
## shown in the first lecture.

## R scripts in this course always start by loading the "wordspace" package,
## so its functions and data sets are available. Sometimes you will also
## need to load further packages (e.g. for visualization).
library(wordspace)


## Let us start with a little recap, and one of the toy models in wordspace
## So, what did we learn today?
# A DSM is a matrix
DSM_TermTermMatrix
## Which contains
# ... targets: rows
rownames(DSM_TermTermMatrix)
# ... contexts: columns
colnames(DSM_TermTermMatrix)
# ... co-occurrence counts (or fancier versions of co-occurrence) for target/context pairs: matrix cells
DSM_TermTermMatrix['cat', 'tail']
# The row corresponding to a target (vector) is the best approximation we have for it its meaning
DSM_TermTermMatrix['cause', ]
DSM_TermTermMatrix['effect', ]
# By the way, columns are vectors, too... 
DSM_TermTermMatrix[, 'important']
DSM_TermTermMatrix[, 'explain']
## (An aside: strings can either be enclosed in single quotes '...' or double
## quotes "..." in R, and there is no functional difference. Gabriella prefers
## single quotes, while Stephanie prefers double quotes, so you'll get used to both
## over the next few days. :-)


## While we can use our co-occurrence matrix as a collection of co-occurrences, 
## which is already an excellent basis to investigate some properties of the
## words we are interested in, what you really want to do with a DSM 
## is typically to COMPARE things

## Look at the vectors: you can already see that some dimensions will be the
## main responsibles for the similarity between 'cause' and 'effect'! 
## That's exactly what cosine similarity is sensitive to (you will learn more
## about it tomorrow): a measure of context overlap.

# ... and this is how you get similarities in wordspace
pair.distances('cause', 'effect', DSM_TermTermMatrix, convert=FALSE)
pair.distances('cause', 'time', DSM_TermTermMatrix, convert=FALSE)
## Don't worry about the extra output starting with attr(,"similarity").
## The "wordspace" package tags some results with additional information
## that other functions can use. In this case, the attribute distinguishes
## between distances (nearest neighbours = lowest values) and similarities
## (nearest neighbours = highest values).

## The angles between vectors are a distance measure that corresponds to cosine
## similarity, and by omitting convert=FALSE, we allow "wordspace" to translate
## the similarity values into angles.
pair.distances('cause', 'effect', DSM_TermTermMatrix)

## The "wordspace" package actually allows us compute the similarities or
## distances between all the target words in a DSM at once, and have them
## summarized in a matrix:
dist.matrix(DSM_TermTermMatrix)
dist.matrix(DSM_TermTermMatrix)['cat', 'dog']
pair.distances('cat', 'dog', DSM_TermTermMatrix)

## Q: Can you print the corresponding matrix of cosine similarities?
## (hint: try round(... , 3) to make the result more readable)

## By the way, remember than when we inspected the column vectors we had the
## impression that certain columns were more useful for comparing the targets
## cause and effect than others?
# Let's see what happens if we focus on them 
DSM_reduced <- DSM_TermTermMatrix[, c('important','likely','explain')]
# We gain a little bit of similarity for cause/effect
pair.distances('cause', 'effect', DSM_TermTermMatrix, convert=FALSE)
pair.distances('cause', 'effect', DSM_reduced, convert=FALSE)
# But we lose a lot of it for cat/dog
pair.distances('cat', 'dog', DSM_TermTermMatrix, convert=FALSE)
pair.distances('cat', 'dog', DSM_reduced, convert=FALSE)
## Nah, we gotta get smarter when selecting dimensions!
## Stephanie will tell you everything about this tomorrow.

## Let's get back at our "full" model, and inspect the nearest neighbors of some
## words we are interested in, i.e. its semantic environment.
nearest.neighbours(DSM_TermTermMatrix, 'cause')
# We can also visualize it
plot(nearest.neighbours(DSM_TermTermMatrix, 'cause', dist.matrix=TRUE))
# Let's also visualize it in our "badly reduced" matrix
plot(nearest.neighbours(DSM_reduced, 'cause', dist.matrix=TRUE))
# Yeah... it is pretty clear we broke it :) 


## Comparing distance and neighbour rank
# Distance (cosine/angle, unless specified otherwise)
pair.distances('cat', 'dog', DSM_TermTermMatrix)
pair.distances('dog', 'cat', DSM_TermTermMatrix)
# Neighbour rank (also based on the default distance unless specified)
#   fwd ... rank of second word among neighbors of the first one
#   bwd ... vice versa (same as switching the order in the pair.distances call)
#           pair.distances(a,b,M,rank="fwd") = pair.distances(b,a,M,rank="bwd")
pair.distances('cat', 'dog', DSM_TermTermMatrix, rank="fwd")
pair.distances('cat', 'dog', DSM_TermTermMatrix, rank="bwd")

# Cats and dogs are best buddies, so we cannot see a difference here
# Animal and cause? 
pair.distances('animal', 'cause', DSM_TermTermMatrix)
pair.distances('cause', 'animal', DSM_TermTermMatrix)
pair.distances('animal', 'cause', DSM_TermTermMatrix, rank="fwd")
pair.distances('animal', 'cause', DSM_TermTermMatrix, rank="bwd")
pair.distances('cause', 'animal', DSM_TermTermMatrix, rank="fwd")

# Let's understand this asymmetry better
plot(nearest.neighbours(DSM_TermTermMatrix, 'animal', dist.matrix=TRUE))
pair.distances('animal','cause', DSM_TermTermMatrix, rank="fwd")
nearest.neighbours(DSM_TermTermMatrix, 'animal')
plot(nearest.neighbours(DSM_TermTermMatrix, 'cause', dist.matrix=TRUE))
pair.distances('cause','animal', DSM_TermTermMatrix, rank="fwd")
nearest.neighbours(DSM_TermTermMatrix, 'cause')
## The very same distance means something very different if we take the
## perspective of "animal" or "cause" -- because "cause" lives in a very
## densely populated subspace.


## Now, enjoy your first steps with more realistically sized DSMs.
## For example, a small model DSM_Vectors is included in the "wordspace"
## package. Here is how you get more info about it:
?DSM_Vectors
dim(DSM_Vectors) # shows number of rows (= targets) and columns (= features)

## Try also one or more of the large-vocabulary models that you can
##   1. download from the course wiki
##   2. save in the folder "models" as discussed before
## For a first test, please download "WP500_DepFilter_Lemma.rda",
## which is the smallest model file at 31.1 MiB size.

## You can now load the model as follows:
load("models/WP500_DepFilter_Lemma.rda", verbose=TRUE)
## verbose=TRUE shows the name of the object that has been loaded frome
## the .rda file. You can now assign it to a shorter name for convenience.
## NB: R doesn't make an actual copy unless you modify one of the two 
##     variables, so don't hesitate to assign shorter names in this way.
model <- WP500_DepFilter_Lemma

## And now we know how to list and plot some neighbors
nearest.neighbours(DSM_Vectors, 'banana_N')
plot(nearest.neighbours(DSM_Vectors, 'banana_N', dist.matrix=TRUE))
nearest.neighbours(model, 'banana_N')
plot(nearest.neighbours(model, 'banana_N', dist.matrix=TRUE))

## You can play around with distance vs. rank to "hunt" for asymmetries if this
## phenomenon is interesting for you. You may notice that some operations become
## slower: this is because the models we are using are getting larger. If you
## want to get an idea of how much computation is happening in the background,
## you can actually time the operations yourself:
system.time(dist.matrix(DSM_TermTermMatrix))
system.time(dist.matrix(DSM_Vectors))
## You better had not try this with the bigger model ... :) 

## But have a look at the output of this: 
system.time(pair.distances("cherry_N", "banana_N", DSM_Vectors))
system.time(pair.distances("cherry_N", "banana_N", DSM_Vectors, rank="fwd"))
## And this:
system.time(pair.distances("cherry_N", "banana_N", model))
system.time(pair.distances("cherry_N", "banana_N", model, rank="fwd"))
## Can you guess why the difference between distance and rank in terms of
## computing time is much bigger in the second case?


## BONUS: Code for reproducing the hieroglyph example in the slides
M <- log2(DSM_HieroglyphsMatrix + 1)
round(M, 3)

pair.distances("dog", "cat", M, convert=FALSE)

nearest.neighbours(M, "dog", n=3)
plot(nearest.neighbours(M, "dog", n=5, dist.matrix=TRUE))
