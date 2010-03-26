## J Toedling                                         1 February 2010

#### some code and functions to work with colourspace sequences

## Translation matrix:
NtoCSmatrix <- matrix(c(0:3, 1:0, 3:2, 2:3, 0:1, 3:0), ncol=4)
dimnames(NtoCSmatrix) <- list("Base1st"=c("A","C","G","T"),
                              "Base2nd"=c("A","C","G","T"))

## auxiliary recursive function for nucleotide to colorspace translation:
nucToCS.rec <- function(chars){
  # recursive function to get matrix elements
  stopifnot(all(chars %in% c("A","C","G","T")))
  lx <- length(chars) # for less typing length(chars)
  if (lx<2)
    return(chars[1]) # need first char for decoding
  else
    return(c(nucToCS.rec(chars[-lx]),
             NtoCSmatrix[chars[lx-1L], chars[lx]]))
}#nucToCS.pair

# test:
nucToCS.rec(c("A","C"))
nucToCS.rec(c("A","C","C","T"))

## actual function for nucleotide to colorspace translation:
translateNucToCS <- function(seqs){
  stopifnot(is.character(seqs))
  splitted <- strsplit(seqs, "")
  cs.seqs <- lapply(splitted, nucToCS.rec)
  return(sapply(cs.seqs, paste, collapse=""))
}#translateNucToCS  

## test:
translateNucToCS(c("ACCAT","GAGAT","TAAAGTGCTTCCATTTTGTGTGT"))


### translating back colour-space reads:
CStoNmatrix <- apply(NtoCSmatrix, 1,
                     function(x) c("A","C","G","T")[order(x)])
dimnames(CStoNmatrix) <- list("PrevBase"=c("A","C","G","T"), "Colour"=0:3)

csToNuc.rec <- function(baseChar, colourChars){
  # recursive function to get matrix elements
  # not: very initial character is not returned here
  stopifnot(length(baseChar)==1,
            baseChar %in% c("A","C","G","T"),
            (length(colourChars)==0) ||
            colourChars[1] %in% c("0","1","2","3"))
  if (length(colourChars)==0)
    return() # need first char for decoding
  else {
    thisNuc <- CStoNmatrix[baseChar, colourChars[1]]
    return(c(thisNuc, csToNuc.rec(thisNuc, colourChars[-1])))
  }
}#csToNuc.rec

csToNuc.rec("T", as.character(0:3))

translateCStoNuc <- function(cs.seqs, returnInitial=FALSE){
  stopifnot(is.character(cs.seqs), is.logical(returnInitial))
  ## what are the initial letters (last primer base,
  ##  needed for decoding ambiguos colour-space)
  initials <- substr(cs.seqs, 1, 1)
  stopifnot(all(initials %in% c("A","C","G","T")))
  splitted <- strsplit(cs.seqs, "")
  cs.seqs <- lapply(splitted, function(x) csToNuc.rec(x[1],x[-1]))
  cs.seqs <- sapply(cs.seqs, paste, collapse="")
  if (returnInitial)
    cs.seqs <- paste(initials, cs.seqs, sep="")
  return(cs.seqs)
}#translateNucToCS  

translateCStoNuc(c("T3002113202013000111111", "T0123"))
