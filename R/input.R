## additional ways to import data into AlignedGenomeIntervals objects

agiFromOneChrBam <- function(from){
  stopifnot(is.list(from),
            all(c("rname","strand", "pos", "seq") %in% names(from)))
  from$seq  <- as.character(from$seq)
  fromWidth <- nchar(from$seq)
  readPos   <- paste(from$rname, from$strand,
                     from$pos, from$pos+fromWidth-1L,
                     from$seq, sep=".")
  tablePos  <- table(sort(readPos))
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
  ## tests:
  stopifnot(file.exists(bamfile), length(bamfile)==1L)

  ## which function to use for each iteration:
  if ("package:parallel" %in% search())
    lFun <- mclapply
  else
    lFun <- lapply

  ## read BAM header:
  H <- scanBamHeader(bamfile)
  allChr <- names(H[[1]]$targets)
  ## what to read from the BAM files:
  selWhat <- c("flag", "rname", "strand", "pos","seq")
  param <- ScanBamParam(simpleCigar=TRUE,
                        reverseComplement=TRUE,
                        what=selWhat)
  ## now read each chromosome seperately:
  #res <- vector("list", length(allChr))
  #for (thisChr in allChr){
  res <- lFun(as.list(allChr), function(thisChr){
    thisRange <- IRangesList(IRanges(1L, H[[1]]$targets[thisChr]))
    names(thisRange) <- thisChr
    theseParams <- initialize(param, simpleCigar=TRUE,
                              flag=scanBamFlag(isUnmappedQuery=FALSE),
                              reverseComplement=TRUE,
                              what=selWhat, which=thisRange)
    S <- scanBam(bamfile, param=theseParams, ...)[[1]]
    if (length(S[[1]])==0){
      return(NULL)}
      #res[[thisChr]] <- NULL; next}
    #res[[thisChr]] <- agiFromOneChrBam(S)
    return(agiFromOneChrBam(S))
  })# for all chromosomes
  names(res) <- allChr
  res <- res[!sapply(res, is.null)]
  ## compute matches over all chromosomes
  res <- do.call("c", res)
  tabSeq <- table(sort(res@sequence))
  res@matches <- as.integer(tabSeq[match(res@sequence, names(tabSeq))])
  return(res)
}#agiFromBam
