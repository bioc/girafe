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
  A <- AlignedRead(sread=DNAStringSet(rep(baseseq, nrow(G)),
                     start=1, width=widths),
                   id=BStringSet(ids),
                   chromosome=factor(seq_name(G)),
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

  seqs <- getSeq(genseq, as.character(chromosome(x)),
                 start=position(x),
                 end=position(x)+width(x)-1L,
                 strand=as.character(x@strand))
  x@sread <- DNAStringSet(seqs)
  return(x)
}#addSequences
