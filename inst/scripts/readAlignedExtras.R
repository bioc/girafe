## function to read in PASS alignment outputs in Gff3 format
##  into objects of class AlignedRead (defined in package ShortRead)

require("ShortRead")
require("genomeIntervals")

readAlignedGff3 <- function(file){
  stopifnot(file.exists(file))
  G <- readGff3(file, isRightOpen=FALSE)
  ids <- as.character(getGffAttribute(G, "Name"))
  widths <- G[,2]-G[,1]+1L
  baseseq <- "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
  A <- AlignedRead(sread=DNAStringSet(baseseq,
                     start=1, width=widths),
                   id=BStringSet(ids),
                   chromosome=factor(seqnames(G)),
                   position=as.integer(G[,1]),
                   strand=factor(G$strand, levels=c("-","+")))
  warning("Added read sequences that only consist of Ns.")
  return(A)
}#readAlignedGff3

### Note that here currently an artifical sequence only consisting of Ns for
##   is used for each read. This invalidates the matches slots of any
##   AlignedGenomeIntervals object the resulting AlignedRead object would
##   be converted to, amongst other problems.
###  Ideal would a subsequent addition of the acutal sequence using the
##    the reference sequences and aligned positions.


## BED file definition here:
##  http://genome.ucsc.edu/goldenPath/help/customTrack.html#BED
readAlignedBed <- function(file, ...){
  stopifnot(file.exists(file))
  G <- scan(file, flush=TRUE, strip.white=TRUE,
            what=list(chrom="",
                      chromStart=0L,
                      chromEnd=0L,
                      name="",
                      score=NULL, # not relevant
                      strand=""),
            ...)
  if (!all(G$strand %in% c("+","-")))
     stop("All strand entries (column 6) must be either '+' or '-'.")
  ids <- G$name
  widths <- G$chromEnd - G$chromStart # no +1! see BED description
  baseseq <- paste(rep("N", max(widths)),collapse="")
  A <- AlignedRead(sread=DNAStringSet(rep(baseseq, length(G$chrom)),
                     start=1, width=widths),
                   id=BStringSet(ids),
                   chromosome=factor(G$chrom),
                   position=G$chromStart + 1L,
                   strand=factor(G$strand, levels=c("-","+")))
  warning("Added read sequences that only consist of Ns.")
  return(A)
} # readAlignedBed 

## function using respective BSgenome package to add position sequences:
addSequences <- function(x, bspackage){
  stopifnot(inherits(x, "AlignedRead"))
  stopifnot(require(bspackage, character.only=TRUE))
  # get the genome sequence object from the package
  genseq <- get(ls(paste("package", bspackage, sep=":"))[1])

  seqs <- getSeq(genseq,
                 gsub("MT","M", as.character(chromosome(x))),
                 start=position(x),
                 end=position(x)+width(x)-1L,
                 strand=as.character(x@strand))
  x@sread <- DNAStringSet(seqs)
  return(x)
}#addSequences


## function to compute some alignment statistics from an
##  AlignedRead object created from a Bowtie mapping
##   (works with bwtmap objects created from Bowtie 0.12.0,
##    and likely some versions before and afterwards)
alignStat <- function(A, type="Bowtie",
                      bspackage="BSgenome.Mmusculus.UCSC.mm9",...)
{
  stopifnot(inherits(A, "AlignedRead"))
  type <- match.arg(type, c("Bowtie", "MAQ"))
  if (type=="Bowtie"){
    ## check Bowtie-specific column headers
    stopifnot(all(c("similar","mismatch") %in% # obj from Bowtie?
                  names(pData(alignData(A)))) )
    ## 1. number of alignments per read:
    nalign <- pData(alignData(A))$similar+1L
    ## 2. mismatch information:
    splitted <-  strsplit(pData(alignData(A))$mismatch,",")
    ## 2a. number of mismatches in each alignment:
    nmismatch <- listLen(splitted)
    ## 2b. position of mismatch in read:
    mmPosPerRead <- lapply(splitted, function(x)
                           as.integer(gsub("\\:.+$","", x))+1L)
  } # if (type=="Bowtie")
  ## new: also support for MAQ
  if (type=="MAQ"){
    stopifnot(all(c("nExactMatch24","nOneMismatch24") %in% # obj from MAQ
                  names(pData(alignData(A)))) )
    ## 1. number of alignments per read:
    ##NB: MAQ only always reports one hit as random
    nalign <- ifelse(pData(alignData(A))$nExactMatch24 > 0L,
                     pData(alignData(A))$nExactMatch24,
                     pData(alignData(A))$nOneMismatch24)
    ## 2. mismatch information:
    ## 2a. get read sequences:
    readSeqs <- ifelse(A@strand=="+", as.character(sread(A)),
                       as.character(reverseComplement(sread(A))))
    ##NB: MAQ reports reverse-complement in case of minus-strand matches
    ## 2b. genomic sequences of those positions
    getSequences <- function(x, bspackage){
      stopifnot(inherits(x, "AlignedRead"))
      stopifnot(require(bspackage, character.only=TRUE))
      # get the genome sequence object from the package
      genseq <- get(ls(paste("package", bspackage, sep=":"))[1])
      seqs <- getSeq(genseq,
                     gsub("MT","M", as.character(chromosome(x))),
                     start=position(x),
                     end=position(x)+width(x)-1L,
                     strand=as.character(x@strand))
      return(seqs)
    }## getSequences
    genomeSeqs <- getSequences(A, bspackage=bspackage)
    ## 2c: generate merged strings, in which '?' indicates mismatches
    merged <- compareStrings(readSeqs, genomeSeqs)
    ## mismatch positions
    mmPosPerRead <- lapply(strsplit(merged, ""),
                           function(z) grep("\\?", z))
    ## 2e: number of mismatches per read
    nmismatch <- listLen(mmPosPerRead)
  } #if (type=="MAQ")
  ## prepare and return result:
  res <- list(n.alignments=table(nalign),
              n.mismatch=table(nmismatch),
              mm.position=table(unlist(mmPosPerRead)))
  class(res) <- c("AlignmentStatistics", class(res))
  return(res)
}#alignStat
# testing:
# testAR <- readAligned("test.bwtmap", type="Bowtie")
# s <- alignStat(testAR)


## for plotting, you can use, e.g., 'lapply' with 'barplot'
setMethod("plot", signature=c(x="AlignmentStatistics",y="missing"),
  function(x, y, ...){
    oldpar <- par(mfrow=c(3,1), font.lab=2,
                  mar=c(4,4,4,0)+0.1, cex.main=1.3)
    on.exit(par(oldpar))
    barplot(x$n.alignments, ylab="Frequency",
            main="Number of alignments per read", ...)
    xpos <- barplot(x$n.mismatch, ylab="Frequency",
                    main="Number of mismatches per alignment", ...)
    percentages <- round(100*x$n.mismatch/sum(x$n.mismatch), digits=1)
    mtext(side=1, line=-1.3, at=xpos, font=2,
          text=paste(percentages,"%",sep=""))
    barplot(x$mm.position, ylab="Frequency",
            main="Positions of mismatches in reads", ...)
    invisible(NULL)
}) #plot

# plot(s, col="cornflowerblue") #test
