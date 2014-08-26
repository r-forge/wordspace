##
## Prepare R data set for SemCor WSD task
##

SemCorWSD <- read.delim("task/semcor_wsd.txt", stringsAsFactors=FALSE, quote="", fileEncoding="utf8")
if (FALSE) {
  str(SemCorWSD)
  with(SemCorWSD, table(paste(target, sense, sep=": "), target))
  unique(SemCorWSD$gloss)
}

## drop infrequent senses (f < 5 instances)
sense.freq <- table(SemCorWSD$sense) 
sense.ok <- names(sense.freq)[sense.freq >= 5]
SemCorWSD <- subset(SemCorWSD, sense %in% sense.ok)
if (FALSE) {
  with(SemCorWSD, table(paste(target, sense, sep=": "), target))
}

## drop uncommon senses (less than 10% of target occurrences)
sense.target <- with(SemCorWSD, prop.table(table(sense, target), margin=2))
sense.ok <- rownames(sense.target)[ apply(sense.target >= 0.1, 1, any) ]
SemCorWSD <- subset(SemCorWSD, sense %in% sense.ok)
if (FALSE) {
  with(SemCorWSD, table(paste(target, sense, sep=": "), target))
}

## drop targets with just a single remaining sense
sense.target <- with(SemCorWSD, table(sense, target))
target.ok <- colnames(sense.target)[ colSums(sense.target > 0) >= 2 ]
SemCorWSD <- subset(SemCorWSD, target %in% target.ok)
if (FALSE) {
  with(SemCorWSD, table(paste(target, sense, sep=": "), target))
}

## prepare final data frame
SemCorWSD <- SemCorWSD[, -2] # drop sentence ID
rownames(SemCorWSD) <- SemCorWSD$id

with(SemCorWSD, table(paste(target, sense, sep=": "), target))
with(SemCorWSD, sort(unique(paste0(target, " ", sense, ": ", gloss))))

str(SemCorWSD)
save(SemCorWSD, file="task/SemCorWSD.rda", compress="xz", compression_level=6)
