###################################################################
### auxiliary functions:
###################################################################

# 1. get a named integer vector of chromosome lengths, based on the
#  organism slot of AlignedGenomeIntervals objects
getChromLengths <- function(x){
  stopifnot(inherits(x, "AlignedGenomeIntervals"))
  if ("chrlengths" %in% slotNames(x) && length(x@chrlengths)>0) {
    chrlens <- x@chrlengths
  } else {
    if ( length(x@organism) != 1 ) {
      warning("Neither 'chrlengths' nor 'organism' are defined in ",
              "object '", deparse(substitute(x)), "'!")
      splitted <- split(x[,2], chromosome(x))
      chrlens <- sapply(splitted, max)
    } else {
      orgpackage <- paste("org",x@organism,"eg.db", sep=".")
      worked <- require(orgpackage, character.only=TRUE)
      if (!worked)
        stop("No package called '",orgpackage,"' found. ",
             "Install this package or check wether the organism annotation ",
             "of ",deparse(substitute(x))," is correct (e.g. 'Mm' or 'Hs')\n")
      chrlens <- get(paste(gsub("\\.db$", "", orgpackage),
                           "CHRLENGTHS", sep=""))
      names(chrlens) <- paste("chr",
                              gsub("^chr","", names(chrlens)),sep="")
      #names(chrlens) <- gsub("MT$", "M", names(chrlens))
      #maybe do something about the chrM/MT nuisance later
    }
  }
  unichrx <- unique(chromosome(x))
  areIn <- unichrx %in% names(chrlens)
  if (!all(areIn))
    warning("Chromosomes '", paste(unichrx[!areIn], collapse=","),
            "' mentioned in object but not found in names of ",
            "vector of chromosome lengths.")
  return(chrlens)
} #getChromLengths


#2. get the read counts for a supplied set of genomic features
# AI: AlignedGenomeIntervals object
# FG: genome.intervals objects of genomic features
# normalise: normalise counts by number of genomic copies per feature
getFeatureCounts <- function(AI, FG, nameColumn="Name",
                             fractionIncluded=1,
                             returnType="AlignedGenomeIntervals",
                             mem.friendly=FALSE)
{
  stopifnot(inherits(AI, "AlignedGenomeIntervals"),
            inherits(FG, "Genome_intervals"),
            nameColumn %in% names(annotation(FG)))
  returnType <- match.arg(returnType,
                          c("AlignedGenomeIntervals","integer"))
  ### get feature names and frequencies:
  feat <- annotation(FG)[[nameColumn]]
  featFreq <- table(feat)
  featMaxFreq <- max(featFreq)
  ## I. first drop intervals with more matches than maximal number
  ##     of feature copies:
  AI <- AI[AI@matches <= featMaxFreq]
  fo  <- fracOverlap(AI, FG, 0, mem.friendly=mem.friendly)
  fo  <- subset(fo, fraction1 >= fractionIncluded)
  # additional filtering for at most as many matches as feature copies:
  fo$featname <- feat[fo$Index2]
  featsPerAlignedInterval <- split(fo$featname, fo$Index1)
  candidateIntervals <- names(featsPerAlignedInterval)
  keep <- sapply(candidateIntervals, function(thisInt){
    overlappingFeats <- featsPerAlignedInterval[[thisInt]]
    ## we keep intervals if
    ## a. they overlap a single feature AND
    ## b. they aligned sequences has at most as many matches
    ##     as the feature has copies in the genome
    return(length(unique(overlappingFeats))==1 &&
           AI@matches[as.integer(thisInt)] <=
           featFreq[ overlappingFeats[1] ])
  })
  keptIntervals <- candidateIntervals[keep]
  fo <- subset(fo, Index1 %in% keptIntervals)
  ### split by mature featNA
  IperInd2 <- split(fo$Index1, fo$Index2)
  NperInd2 <- sapply(IperInd2, function(z)
                     sum(AI[z]@reads))
  ## option to return named integer vector of counts here:
  if (returnType=="integer"){
    names(NperInd2) <- feat[as.integer(names(NperInd2))]
    return(NperInd2)
  }
  ## prepare return object:
  Mo <- FG[as.integer(names(IperInd2))]
  MoFeatNames <- annotation(Mo)[[nameColumn]]
  xo <- AlignedGenomeIntervals(start=Mo[,1], end=Mo[,2],
                               chromosome=chromosome(Mo),
                               strand=strand(Mo),
                               sequence=MoFeatNames,
                               reads=as.integer(NperInd2),
                               matches=featFreq[MoFeatNames],
                               id=MoFeatNames
                               )
  return(xo)
}#getFeatureCounts

# auxiliary function to compute summed read number and GC content
#  over fixed-width windwos, faster but more restricted than "perWindow"
# before: used only + strand, but probably not necessary since
#  reverse complement has same GC content
windowCountAndGC <- function(G, chr, bspackage, wsize=1000L,
                             drop.zero=TRUE, verbose=TRUE)
{
  stopifnot(inherits(G, "AlignedGenomeIntervals"))
  G <- G[chromosome(G)==chr]
  if (nrow(G)==0){
    warning("No aligned genome intervals on chr",chr,"!")
    return(data.frame(start=integer(0), end=integer(0),
                      n.reads=integer(0), GC=numeric(0)))
  }
  G[,1] <- G[,2] <- rowMeans(G[,1:2])
  #G <- reduce(G, method="exact")
  if (verbose) cat("computing coverage...\n")
  covChr <- coverage(G)   # , byStrand=TRUE)
  stopifnot(length(covChr)==1L)
  covPlus <- covChr[[1]]  #grep("\\+",names(covChr))]]
  starts <- seq.int(1L, length(covPlus)-wsize+1L, by=wsize)
  ends   <- seq.int(wsize, length(covPlus), by=wsize)
  if (verbose) cat("summing read counts over running windows...\n")
  agPlus  <- aggregate(covPlus, FUN=sum,
                       start=starts, end=ends)
  W <- data.frame(start=starts, end=ends,n.reads=agPlus)
  if (drop.zero) W <- subset(W, n.reads>0)
  stopifnot(require(bspackage, character.only=TRUE))
  # get the genome sequence object from the package
  if (verbose) cat("computing GC content...\n")
  genseq <- get(ls(paste("package", bspackage, sep=":"))[1])
  ## normalise chromosome name:
  normchr <- gsub("MT$","M", chr)
  stopifnot(normchr %in% seqnames(genseq))
  Wseqs <- Views(unmasked(genseq[[normchr]]),
                 start=W$start, end=W$end)
  Waf <- alphabetFrequency(Wseqs, as.prob=TRUE, baseOnly=TRUE)
  W$GC <- round(rowSums(Waf[,c("C","G"), drop=FALSE]),digits=3)
  if (verbose) cat("Done.\n")
  return(W)
} #windowCountAndGC

## usage example:
#exW <-  windowCountAndGC(exAI, chr="chrX",
#                         bspackage="BSgenome.Mmusculus.UCSC.mm9")


## auxiliary function to compute consensus with number of aligned
##  reads as weights
weightedConsensusMatrix <- function(seqs, weights, shift=NULL,
                                    baseLetters=c("A", "C", "G", "T", "N"))
{
  stopifnot(is.character(seqs), is.integer(weights),
            length(seqs)==length(weights))
  ## any shift of sequences specified?
  if (is.null(shift))
    shift <- integer(length(seqs))
  else
    stopifnot(length(shift)==length(seqs))
  ## what is the final length of the consensus matrix?
  maxreadlen <- max(nchar(seqs)+shift)
  ## initialise result matrix:
  DW <- matrix(0L, nrow=length(baseLetters), ncol=maxreadlen)
  rownames(DW) <- baseLetters

  ## split sequences into individual letters
  splitted.seqs <- strsplit(seqs, split="")
  ## loop over individual sequences and fill in matrix:
  for (i in seq.int(length(seqs))){
    x.coords <- match(splitted.seqs[[i]], baseLetters)
    y.coords <- pmax.int(0L,
                         seq.int(length(splitted.seqs[[i]])) + shift[[i]])
    DW[cbind(x.coords, y.coords)] <-
      DW[cbind(x.coords, y.coords)] + weights[[i]]
  }
  return(DW)
}#weightedConsensusMatrix

### test new function:
### Align following sequences with weights:
###   ACATT       1
###    CGTTA     10
###     TTG       3
###  GACATT       4

#dweights <- c(1L, 10L, 3L, 4L)
#d <- c("ACATT","CGTTA", "TTG", "GACATT")
#dshifts <- c(0L, 1L, 2L, -1L)
#weightedConsensusMatrix(d, dweights, shift=dshifts)
#consensusString(.Last.value, ambiguityMap="N")
