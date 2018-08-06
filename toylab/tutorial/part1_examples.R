library(wordspace)

DSM_HieroglyphsMatrix
DSM_TermTermMatrix
DSM_TermContextMatrix

M <- log2(DSM_HieroglyphsMatrix + 1)
round(M, 3)

pair.distances("dog", "cat", M, convert=FALSE)

nearest.neighbours(M, "dog", n=3)
plot(nearest.neighbours(M, "dog", n=3, dist.matrix=TRUE))
