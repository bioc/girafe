###################################################################
### auxiliary functions:
###################################################################

# 1. get a named integer vector of chromosome lengths, based on the
#  organism slot of AlignedGenomeIntervals objects
getChromLengths <- function(x){
  stopifnot(inherits(x, "AlignedGenomeIntervals"))
  if ( length(x@organism) != 1 )
    stop("Organism of '",deparse(substitute(x)),
         "' not defined! Set with organism <-")
  orgpackage <- paste("org",x@organism,"eg.db", sep=".")
  worked <- require(orgpackage, character.only=TRUE)
  if (!worked)
    stop("No package called '",orgpackage,"' found.",
         "Install this package or check wether the organism annotation",
         "of ",deparse(substitute(x))," is correct (e.g. 'Mm' or 'Hs')\n")
  chrlens <- get(paste(gsub("\\.db$", "", orgpackage),
                       "CHRLENGTHS", sep=""))
  names(chrlens) <- paste("chr",
                          gsub("^chr","", names(chrlens)),sep="")
  return(chrlens)
} #getChromLengths


#2. get the read counts for a supplied set of genomic features
# AI: AlignedGenomeIntervals object
# FG: genome.intervals objects of genomic features
# normalise: normalise counts by number of genomic copies per feature
getFeatureCounts <- function(AI, FG, nameColumn="Name",
                             fractionIncluded=1,
                             returnType="AlignedGenomeIntervals")
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
  fo  <- fracOverlap(AI, FG, 0)
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
