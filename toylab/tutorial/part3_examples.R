library(wordspace)
options(width=60) # so output fits on slides

M <- DSM_Vectors
nearest.neighbours(M, "walk_V")

library(wordspaceEval)
head(TOEFL80)

eval.multiple.choice(TOEFL80, M) # manhattan is slightly better

RG65[seq(1,61,5), ]
head(WordSim353)

eval.similarity.correlation(RG65, M, convert=FALSE)
plot(eval.similarity.correlation(RG65, M, details=TRUE, convert=FALSE))

ESSLLI08_Nouns[seq(1,40,5), ]

eval.clustering(ESSLLI08_Nouns, M)
