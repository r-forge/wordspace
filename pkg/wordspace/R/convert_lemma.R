convert.lemma <- 
function (lemma, format=c("CWB", "BNC", "DM", "HW", "HWLC")) {
  format <- match.arg(format)
  lemma <- as.character(lemma)

  tmp <- strsplit(lemma, "_", fixed=TRUE)
  errors <- sapply(tmp, function (x) !(length(x) == 2 && nchar(x[2]) == 1))
  if (any(errors)) stop("invalid lemma string(s): ", paste(lemma[errors], collapse=", "))
  hw <- sapply(tmp, function (x) x[1])
  pos <- sapply(tmp, function (x) x[2])
  
  if (format == "BNC") {
    penn <- c("N",     "V",    "J",   "R",   "D"  )
    bnc  <- c("SUBST", "VERB", "ADJ", "ADV", "ART", "UNC")
    pos.idx <- match(pos, penn, nomatch=length(bnc)) # map Penn POS code to BNC POS code, setting to UNC if not found
    bncpos <- bnc[pos.idx]
  }
  
  switch(format,
    CWB = lemma,
    BNC = paste(tolower(hw), bncpos, sep="_"),
    DM  = paste(hw, tolower(pos), sep="-"),
    HW  = hw,
    HWLC = tolower(hw)
  )
}