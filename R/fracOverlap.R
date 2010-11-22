fracOverlap <- function(I1, I2, min.frac=0.0, both=TRUE, mem.friendly=FALSE)
{
  stopifnot(inherits(I1,"Genome_intervals"),
            inherits(I2,"Genome_intervals"))
  if ( inherits(I1, "AlignedGenomeIntervals") || inherits(I2, "AlignedGenomeIntervals") )
    ov <- interval_overlap(I1,I2, mem.friendly=mem.friendly)
  else
    ov <- interval_overlap(I1,I2)
  # get base pair overlap
  lens <- listLen(ov)
  overlap1 <- rep.int(1:length(ov), lens)
  overlap2 <- unlist(ov, use.names=FALSE)
  left <- pmax.int(I1[overlap1,1], I2[overlap2,1])
  right <- pmin.int(I1[overlap1,2], I2[overlap2,2])
  stopifnot(all(right >= left))
  bases <- right-left+1L
  len1 <- I1[overlap1,2]-I1[overlap1,1]+1L
  len2 <- I2[overlap2,2]-I2[overlap2,1]+1L
  frac1 <- round(bases/len1, digits=2)
  frac2 <- round(bases/len2, digits=2)
  #frac <- round(bases/pmin(len1, len2), digits=2)
  #min.len <- pmin.int(I1[overlap1,2]- I1[overlap1,1]+1,
  #                    I2[overlap2,2]- I2[overlap2,1]+1)
  #frac <- round(bases/min.len, digits=2)
  res <- data.frame("Index1"=overlap1, "Index2"=overlap2,
                    "n"=bases, "fraction1"=frac1, "fraction2"=frac2)
  if (min.frac > 0.0) # fraction specified?
    if (both) # shall both overlap partners meet that requirement?
      res <- subset(res, fraction1 >= min.frac & fraction2 >= min.frac)
    else # or is one sufficient?
      res <- subset(res, fraction1 >= min.frac | fraction2 >= min.frac)
  return(res)
}# fracOverlap
