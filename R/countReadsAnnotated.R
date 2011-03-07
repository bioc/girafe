# function to count reads per type of genomic feature (gene, miRNA, etc.)

# I. set type order, each read should only be counted once, e.g. if a read
#  is annotated for a miRNA inside a gene, it should only count for
#  miRNAs not for genes too); order from most specific to least specific
# II. normalise read counts by number of matches per read
genomeFeatureClassOrder <-
  list("miRNA"=c("miRNA","microRNA"),
       "7SK"=c("7SK","7SK-RNA","7SK_RNA", "RN7SK"),
       "Y_RNA"=c("yRNA","Y-RNA","Y_RNA"),
       "rRNA"=c("rRNA", "ribosomal-RNA"),
       "snoRNA"="snoRNA",
       "snRNA"="snRNA", "tRNA"=c("tRNA","transfer-RNA"),
       "ribozyme"="ribozyme",
       "vRNA"=c("Vault","vRNA","Vault-RNA","svRNA","VAULT"),
       "TERC"=c("Telomerase-vert","TERC"),
       "NRON"="NRON",
       "lincRNA"=c("lincRNA", "lncRNA", "long-ncRNA"),
       "ncRNA"=c("ncRNA", "NCRNA","non-coding"),
       "pseudogene"=c("pseudogene","Pseudogene"),
       "exon"="exon", "intron"="intron",
       "gene"=c("gene","Gene", "CDS"),
       "LINE"=c("LINE", "LINE_repeat"),
       "repeat"=c("repeat", "SINE", "LTR",
         "Satellite", "Simple_repeat", "other_repeat"),
       "other"="other"
       )

countReadsAnnotated <- function(GI, M, typeColumn="type", fractionGI=0.7,
                                mem.friendly=FALSE, showAllTypes=FALSE)
{
  stopifnot(inherits(GI, "Genome_intervals"), #"AlignedGenomeIntervals"),
            inherits(M, "Genome_intervals"))
  ## which function to use for each iteration:
  if ("package:multicore" %in% search())
    lFun <- mclapply
  else
    lFun <- lapply
  classOrder <- genomeFeatureClassOrder # separate to allow easy editing
    mClass <- factor(M@annotation[[typeColumn]])
  levels(mClass) <- classOrder
  if (any(is.na(mClass))){
    warning('Unrecognised classes in class factor replaced by "other".')
    mClass[is.na(mClass)] <- "other"
  }
  if (showAllTypes) print(table(mClass))
  fo <- fracOverlap(GI, M, mem.friendly=mem.friendly)
  fo <- subset(fo, fraction1 >= fractionGI)
  splitted <- split(fo$Index2, fo$Index1)
  classPerInt <- unlist(lFun(splitted, function(x)
                             return(sort(mClass[x])[1])), use.names=FALSE)
  if (inherits(GI, "AlignedGenomeIntervals")){
    nreadsPerClass <- sapply(split(as.integer(names(splitted)), classPerInt),
                             function(z) sum(reads(GI)[z]/matches(GI)[z]) )
    nreadsPerClass <- round(nreadsPerClass)
    nreadsPerClass["unannotated"] <-
      round(sum(reads(GI)/matches(GI)) - sum(nreadsPerClass))
  } else {
    nreadsPerClass <- sapply(split(as.integer(names(splitted)), classPerInt),
                             function(z) length(z) )
    nreadsPerClass["unannotated"] <- nrow(GI) - sum(nreadsPerClass)
  }
  nreadsPerClass <- nreadsPerClass[nreadsPerClass!=0L]
  return(nreadsPerClass)
}#countReadsAnnotated
