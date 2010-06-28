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


## adapter trimming in colourspace:
trimCSAdapter <- function(fq, adapter, match.score=1,
                          mismatch.score=-1,score.threshold=2)
{
  ## argument 'fq' : BStringSet of input colourspace sequences:
  ### check arguments:
  stopifnot(inherits(fq,"BStringSet"),
            inherits(adapter, "BString")||
            (inherits(adapter, "character") && length(adapter)==1))
  if (is.character(adapter))
    adapter <- BString(adapter)

  ## is first letter of adapter a nucleotide for decoding?
  if (as.character(substr(adapter,1L,1L)) %in% c("A","C","G","T"))
    adapter <- substr(adapter, 3L, nchar(adapter))
    # cannot match those
  n.adjust <- 1L # we need to substract one more for the unknown start of the adapter
  read.length <- nchar(fq)
  #construct colourspace substitution matrix:
  mat <- matrix(mismatch.score, nrow=9, ncol=9)
  dimnames(mat) <- list(c(0:3,".","A","C","G","T"),
                        c(0:3,".","A","C","G","T"))
  diag(mat) <- match.score
  mat[".", 1:4] <- match.score
  mat[1:4,"."] <- match.score

  ## do pairwise alignment
  pa <- pairwiseAlignment(fq, adapter, type="overlap",
                          substitutionMatrix=mat,
                          gapOpening=0, gapExtension=-Inf)
  
  areCompleteOverlap <- (score(pa) > score.threshold) &
                        (start(pattern(pa))==1L) &
                        (end(pattern(pa))==read.length)
  kstarts <- integer(length(fq))+1L
  kends <- ifelse(score(pa)<score.threshold, read.length,
                  ifelse(end(pattern(pa))==read.length,
                         end(pattern(pa))-width(pattern(pa))-n.adjust,
                         read.length))
  if (sum(areCompleteOverlap)>0){
    kstarts[areCompleteOverlap] <- 1L
    kends[areCompleteOverlap] <- read.length
  }
  fq2 <- narrow(fq, start=kstarts, end=kends)
  stopifnot(length(fq)==length(fq2))
  return(fq2)
}#trimCSAdapter


## test:
seqs.cs <- BStringSet(c("shrek primer1 rev"=
                        "T233020103031311231200032032222031220201003000312",
                        "mmu-miR-294"="T3002113202002000111111",
                        "mmu-miR-294 plus ten nt primer with 1 error"=
                        "T30021132020020001111112330201032",
                        "mmu-miR-293 with 20 nt primer with 3 errors"=
                        "T321130331222100113211123302020303132223120",
                        "mmu-miR-293 with 20 nt primer with 10 errors"=
                        "T321130331222100113211111102211312132223122"))
trimCSAdapter(seqs.cs, adapter=seqs.cs[[1]])

