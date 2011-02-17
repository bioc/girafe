
### function to remove 3' adapter remnants from reads
trimAdapter <- function(fq, adapter, match.score=1, mismatch.score=-1,
                        score.threshold=2)
{
  stopifnot(inherits(fq,"ShortReadQ"),
            inherits(adapter, "DNAString")|inherits(adapter, "character"))
  if (is.character(adapter))
    adapter <- DNAString(adapter)
  read.length <- unique(width(fq))
  if (length(read.length)>1)
    stop(paste("Expected all reads in object",
               deparse(substitute(fq)),
               "to be of the same lengths! Found lengths:",
               paste(read.length, collapse=", "),"\n"))
  mat <- nucleotideSubstitutionMatrix(match=match.score,
                                      mismatch=mismatch.score)
  pa <- pairwiseAlignment(sread(fq), adapter, type="overlap",
                          substitutionMatrix=mat,
                          gapOpening=0, gapExtension=-Inf)
  areCompleteOverlap <- (score(pa) >= score.threshold) &
                        (start(pattern(pa))==1) &
                        (end(pattern(pa))==read.length)
  kstarts <- integer(length(fq))+1L
  kends <- ifelse(score(pa)<score.threshold, read.length,
                  ifelse(end(pattern(pa))==read.length,
                         end(pattern(pa))-width(pattern(pa)),read.length))
  kwidths <- kends - kstarts + 1L
  if (sum(areCompleteOverlap)>0)
    kwidths[areCompleteOverlap] <- 0L
  fq2 <- narrow(fq, start=kstarts, width=kwidths)
  stopifnot(length(fq)==length(fq2))
  return(fq2)
}#trimAdapter
