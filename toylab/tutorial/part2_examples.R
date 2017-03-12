library(wordspace)

TC <- DSM_TermContext
head(TC, Inf)

TT <- DSM_TermTerm
head(TT, Inf)

TT$rows
TT$cols
TT$globals$N

options(digits=3)
dsm.score(TT, score="MI", sparse=FALSE, matrix.only=TRUE)

dsm.score(TT, score="MI", matrix=TRUE)
dsm.score(TT, score="simple-ll", transform="log", matrix=TRUE)

TT <- dsm.score(TT, score="frequency", transform="log")

pair.distances(c("cat","cause"), c("animal","effect"), TT, method="euclidean")

dist.matrix(TT, method="euclidean")
dist.matrix(TT, method="minkowski", p=4)

rowNorms(TT$S, method="euclidean")
TT <- dsm.score(TT, score="freq", transform="log", normalize=TRUE, method="euclidean")
rowNorms(TT$S, method="euclidean")
dist.matrix(TT, method="euclidean")

TT2 <- dsm.projection(TT, n=2, method="svd")
TT2

x <- TT2[, 1]
y <- TT2[, 2]
plot(TT2, pch=20, col="red", xlim=extendrange(x), ylim=extendrange(y))
text(TT2, rownames(TT2), pos=3)

subset(DSM_VerbNounTriples_BNC, noun == "dog" & verb == "walk")

tri <- subset(DSM_VerbNounTriples_BNC, rel == "obj")
VObj <- dsm(target=tri$noun, feature=tri$verb, score=tri$f, raw.freq=TRUE)
VObj
head(VObj$rows, 20)

VObj <- dsm.score(VObj, score="MI", normalize=TRUE)
nearest.neighbours(VObj, "dog")

nearest.neighbours(VObj, "dog", method="manhattan")

VObj50 <- dsm.projection(VObj, n=50, method="svd")
nearest.neighbours(VObj50, "dog", method="euclidean")
