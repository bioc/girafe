## additional ways to import data into AlignedGenomeIntervals objects

agiFromOneChrBam <- function(from){
  stopifnot(is.list(from),
            all(c("rname","strand", "pos", "seq") %in% names(from)))
  from$seq  <- as.character(from$seq)
  fromWidth <- nchar(from$seq)
  readPos   <- paste(from$rname, from$strand,
                     from$pos, from$pos+fromWidth-1L,
                     from$seq, sep=".")
  tablePos  <- table(readPos)
  ## get index of unique positions
  idx       <- which(!duplicated(readPos))
  ### condense multiply mentioned alignments; element 'posfreq' preserves
  GI <- new("AlignedGenomeIntervals",
            .Data = cbind(from$pos[idx],
                          from$pos[idx]+fromWidth[idx]-1L),
            annotation=data.frame(
              "seq_name"      = from$rname[idx],
              "strand"        = from$strand[idx],
              "inter_base"    = vector("logical", length(idx))),
            sequence = from$seq[idx],
            reads    = as.integer(tablePos[match(readPos[idx],
              names(tablePos))]),
            matches  = rep.int(1L, length(idx)),
            id       = character(length(idx)))
  stopifnot(validObject(GI))
  return(GI)
}#agiFromOneChrBam 


agiFromBam <- function(bamfile, ...)
{
  require("Rsamtools")
  ## tests:
  stopifnot(file.exists(bamfile), length(bamfile)==1)
  ## read BAM header:
  H <- scanBamHeader(bamfile)
  allChr <- names(H[[1]]$targets)
  ## what to read from the BAM files:
  bamWhat <- c("flag", "rname", "strand", "pos","seq")
  param <- ScanBamParam(simpleCigar=TRUE,
                        reverseComplement=TRUE,
                        what=bamWhat)
  ## now read each chromosome seperately:
  res <- vector("list", length(allChr))
  names(res) <- allChr
  for (thisChr in allChr){
    thisRange <- RangesList(IRanges(1, H[[1]]$targets[thisChr]))
    names(thisRange) <- thisChr
    theseParams <- initialize(param, simpleCigar=TRUE,
                              flag=scanBamFlag(isUnmappedQuery=FALSE),
                              reverseComplement=TRUE,
                              what=bamWhat, which=thisRange)
    S <- scanBam(bamfile, param=theseParams, ...)[[1]]
    if (length(S[[1]])==0){
      res[[thisChr]] <- NULL; next}
    ## TO DO: build AlignedGenomeIntervals here
    res[[thisChr]] <- agiFromOneChrBam(S)
  }# for all chromosomes
  ## compute matches over all chromosomes
  res <- do.call("c", res)
  tabSeq <- table(res@sequence)
  res@matches <- as.integer(tabSeq[match(res@sequence, names(tabSeq))])
  return(res)
}#agiFromBam