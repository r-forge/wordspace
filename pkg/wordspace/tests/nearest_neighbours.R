## Validate nearest.neighbours() implementation

library(wordspace)

## example data from tutorial: hieroglyphs matrix
M <- log(DSM_HieroglyphsMatrix + 1) # similarity measure in tutorial is cosine on log frequencies

## nearest neighbours of mystery word "dog" as shown in tutorial
nn <- nearest.neighbours(M, "dog", convert=FALSE)
stopifnot(all(names(nn) == c("cat", "pig", "cup", "boat", "banana", "knife")))
stopifnot(round(nn["cat"], 3) == 0.961) # verify similarity values in tutorial slides
stopifnot(round(nn["pig"], 3) == 0.939)
stopifnot(round(nn["knife"], 3) == 0.770)

## extract same nearest neighbours with vector as target
nn2 <- nearest.neighbours(M, M["dog", ], convert=FALSE)
stopifnot(all.equal(nn2[-1], nn)) # nn2 also finds "dog" as first NN

## return list of NN vectors for multiple targets
res <- nearest.neighbours(M, c("dog", "rat", "cat"), skip.missing=TRUE)
stopifnot(length(res) == 2 && all(names(res) == c("dog", "cat")))

## return full distance matrix
DM1 <- nearest.neighbours(M, "dog", n=3, dist.matrix=TRUE)
DM2 <- dist.matrix(M, terms=c("dog", "cat", "pig", "cup"))
stopifnot(all(DM1 == DM2)) # nearest.neighbours() calls dist.matrix(), so results should be exactly the same

## access matrix by columns
cM <- t(M)
nn <- nearest.neighbours(cM, "dog", byrow=FALSE, convert=FALSE)
stopifnot(all(names(nn) == c("cat", "pig", "cup", "boat", "banana", "knife")))
stopifnot(round(nn["cat"], 3) == 0.961)
cDM1 <- nearest.neighbours(cM, "dog", n=3, byrow=FALSE, dist.matrix=TRUE)
stopifnot(all(cDM1 == DM2))

## find nearest neighbours in pre-computed distance matrix
pDM <- dist.matrix(M, method="cosine")
nn3 <- nearest.neighbours(pDM, "dog")
stopifnot(all(names(nn3) == names(nn)))
stopifnot(all.equal(cos(nn3 * pi / 180), nn)) # angular distance vs. cosine similarity
nn4 <- nearest.neighbours(pDM, "dog", byrow=FALSE) # pDM is symmetric, so nn4 == nn3
stopifnot(all(names(nn4) == names(nn3)))

## search arbitrary similarity matrix for "nearest neighbours"
SMd <- as.distmat(DSM_TermTerm, similarity=TRUE)          # dense similarity matrix, extracted from DSM
SMs <- as.distmat(DSM_TermContextMatrix, similarity=TRUE) # sparse similarity matrix

res <- nearest.neighbours(SMd, "time", n=3)
stopifnot(length(res) == 3 && all(res == c(134, 100, 94)) && all(names(res) == c("kill", "likely", "important")))
res <- nearest.neighbours(SMd, "feed", n=3, byrow=FALSE)
stopifnot(length(res) == 3 && all(res == c(86, 32, 29)) && all(names(res) == c("animal", "dog", "time")))

res <- nearest.neighbours(SMs, "animal", n=5) # should only return 3 nonempty cells
stopifnot(length(res) == 3 && all(res == c(15, 10, 2)) && all(names(res) == c("Pet", "Feral", "Felidae")))
res <- nearest.neighbours(SMs, "Back_Pain", n=3, byrow=FALSE) # only 2 nonempty cells
stopifnot(length(res) == 2 && all(res == c(6, 1)) && all(names(res) == c("cause", "reason")))

