##########################################################################
### exporting methods for AlignedGenomeIntervals
###########################################################################

### auxiliary function for writing to file
writeExportData <- function(dat, attribs, format, filename,
                            writeHeader=TRUE, append=FALSE)
{
  trackDef <- paste('track type=', format, sep='')
  for (i in 1:length(attribs))
    trackDef <- paste(trackDef,' ',names(attribs)[i],'="',
                      attribs[[i]],'"', sep="")
  if (writeHeader)
    cat(trackDef, "\n", file=filename, sep="", append=append)
  suppressWarnings(write.table(dat, file=filename, sep=" ",
                               append=(append || writeHeader),
                               col.names=FALSE, row.names=FALSE,
                               quote=FALSE))
  message(paste("Result written to file:\n", filename,"\n"))
}#writeExportData


setMethod("export", signature("AlignedGenomeIntervals", "character", "character"), function(object, con, format, ..., writeHeader=TRUE, append=FALSE){
  format <- match.arg(format, c("wiggle_0", "bed", "bedGraph",
                                "bedStrand", "bedGraphStrand"))
  # maybe others later (esp. bigWig, bigBed may be of interest)
  fileSuffix <- switch(format,"wiggle_0"="wig","bed"="bed", "bedGraph"="beg",
                       "bedStrand"="bed", "bedGraphStrand"="beg", "txt")

  further.args <- lapply(as.list(match.call(expand.dots=FALSE)[["..."]]),eval)

  ## default track attributes
  attribs <- list(name="aligned data", description="aligned data",
                  color="200,100,00", altColor="0,100,200",
                  visibility="full", autoScale="on")
  if (length(further.args)>0)
    for (i in 1:length(further.args))
     attribs[names(further.args)[i]] <- further.args[[i]]

  ## chromosome information and chromsome replacements
  chroms <- as.character(seq_name(object))
  chroms <- gsub("chrMT","chrM", chroms)

  ## any scores for the intervals?
  ### for bedGraph export it may be preferable to use the number of reads
  ###  per interval anyway. In that case set the score of the object to
  ### NULL before exporting.
  if (is.null(score(object)) || is.na(score(object)))
    iscores <- reads(object)
  else
    iscores <- score(object)
 
  if (format=="wiggle_0"){
    # http://genome.ucsc.edu/goldenPath/help/wiggle.html
    ### need one file per chromosome?
    ### no, but same span throughout the entire data set
    # 1. compute coverage over all chromosomes
    covA <- coverage(object, byStrand=FALSE)
    ### start output:
    resultFile <- paste(gsub("\\..+$","", con),fileSuffix, sep=".")
    trackDef <- paste('track type=', format, sep='')
    ## include zero-line unless specified otherwise:
    if (!is.element("alwaysZero", names(attribs)))
      attribs[["alwaysZero"]] <- "on"
    ## write track header line
    for (i in 1:length(attribs))
      trackDef <- paste(trackDef,' ',names(attribs)[i],'="',
                        attribs[[i]],'"', sep="")
    cat(trackDef, "\n", file=resultFile, sep="")
    for (chr in names(covA)){
      # go through data for each chromosome
      startsA <- c(1, cumsum(runLength(covA[[chr]])))
      # last element actually isn't a start:
      startsA <- startsA[-length(startsA)]
      valA    <- runValue(covA[[chr]])
      lenA    <- runLength(covA[[chr]])
      toWrite <-  valA != 0
      if (any(toWrite)){
        cat("variableStep chrom=",chr,"\n",sep="",
            file=resultFile, append=TRUE)
        for (i in which(toWrite)){
          cat(paste(startsA[i]+(0:(lenA[i]-1))," ",
                    rep(valA[i], lenA[i]),"\n", sep=""),
              sep="", file=resultFile, append=TRUE)
        }
      }
    }#for (chr in names(covA))
    cat("Result written to file:\n", resultFile,"\n")
    dat <- covA
  }
  if (format=="bed") {
    ## five columns: 1.chrom 2.chromStart 3.chromEnd
    # 4. name (displayed next to element, use for number of reads
    # 5. score (not used) 6. strand (either '+' or '-')
    ## bed coordinates are zero-based, half-open
    dat <- data.frame(chr=chroms,
                      chromStart=sprintf("%.0f", object[,1]
                        - ifelse(object@closed[,1], 1L, 0L)),
                      chromEnd=sprintf("%.0f", object[,2]
                        - ifelse(object@closed[,2], 0L, 1L)),
                      name=sprintf("%.0f", object@reads),
                      score=iscores, # see above
                      strand=as.character(strand(object)),
                      stringsAsFactors=FALSE)
    resultFile <- paste(gsub("\\..+$","", con),fileSuffix, sep=".")
    writeExportData(dat, attribs, "bed", resultFile,
                    writeHeader=writeHeader, append=append)
  }
  if (format=="bedStrand") {
    ## see 'bed', but one file per strand
    dat <- data.frame(chr=chroms,
                      chromStart=sprintf("%.0f", object[,1]
                        - ifelse(object@closed[,1], 1L, 0L)),
                      chromEnd=sprintf("%.0f", object[,2]+
                        - ifelse(object@closed[,2], 0L, 1L)),
                      name=sprintf("%.0f", object@reads),
                      score=iscores, # see above
                      strand=as.character(strand(object)),
                      stringsAsFactors=FALSE)
    for (strand in c("-", "+")){
      resultFile <- paste(gsub("\\..+$","", con),"_",
                          ifelse(strand=="-", "minus", "plus"),
                          ".", fileSuffix, sep="")
      attribs2 <- attribs
      attribs2[["name"]] <-
        paste(attribs[["name"]], ifelse(strand=="-", "minus", "plus"),sep="_")
      attribs2[["description"]] <-
        paste(attribs[["description"]],", ", strand," strand",sep="")
      dat2 <- dat[as.character(strand(object))==strand, , drop=FALSE]
      writeExportData(dat2, attribs2, "bed", resultFile,
                      writeHeader=writeHeader, append=append)
    }#for (strand in c("-", "+"))}
  }
  if (format=="bedGraph") {
    ## what scores to plot in bedGraph, normally the number of reads
     
    # from definition of the format at:
    # http://genome.ucsc.edu/goldenPath/help/bedgraph.html
    # coordinates are zero-based, half-open.
    dat <- data.frame(chr=chroms,
                      start=sprintf("%.0f", object[,1]
                        - ifelse(object@closed[,1], 1L, 0L)),
                      end=sprintf("%.0f", object[,2]
                        - ifelse(object@closed[,2], 0L, 1L)),
                      score=formatC(iscores, format="fg"),
                      stringsAsFactors=FALSE)
    resultFile <- paste(gsub("\\..+$","", con),fileSuffix, sep=".")
    writeExportData(dat, attribs, "bedGraph", resultFile,
                    writeHeader=writeHeader, append=append)
  }
  if (format=="bedGraphStrand") {
    # like 'bedGraph' but one file per strand
    dat <- data.frame(chr=chroms,
                      start=sprintf("%.0f", object[,1]
                        - ifelse(object@closed[,1], 1L, 0L)),
                      end=sprintf("%.0f", object[,2]
                        - ifelse(object@closed[,2], 0L, 1L)),
                      score=formatC(iscores, format="fg"),
                      stringsAsFactors=FALSE)
    for (strand in c("-", "+")){
      resultFile <- paste(gsub("\\..+$","", con),"_",
                          ifelse(strand=="-", "minus", "plus"),
                          ".", fileSuffix, sep="")
      attribs2 <- attribs
      attribs2[["name"]] <-
        paste(attribs[["name"]],strand,"strand")
      attribs2[["description"]] <-
        paste(attribs[["description"]],strand,"strand")
      dat2 <- dat[as.character(strand(object))==strand, , drop=FALSE]
      writeExportData(dat2, attribs2, "bedGraph", resultFile,
                      writeHeader=writeHeader, append=append)
    }#for (strand in c("-", "+"))}
  }# bedGraphStrand
  invisible(dat)
})

# "..." not last because everything after "..." needs exact matching of argument name
setMethod("export", signature("Genome_intervals", "character", "ANY"),
          function(object, con, format="bed", ..., # ... not last
                   nameColumn=NULL, writeHeader=TRUE, append=FALSE)
          {
            if (!is.null(nameColumn))
              stopifnot(nameColumn %in% names(annotation(object)))
            
            ## additional arguments for UCSC track definition line?
            further.args <-
              lapply(as.list(match.call(expand.dots=FALSE)[["..."]]),eval)
            ## default track attributes
            attribs <- list(name="intervals data",
                            description="intervals data",
                            color="200,100,00", altColor="0,100,200",
                            visibility="full", autoScale="on")
            if (length(further.args)>0)
              for (i in 1:length(further.args))
                attribs[names(further.args)[i]] <- further.args[[i]]
 
            ## chromosome information and chromsome replacements
            chroms <- as.character(seq_name(object))
            chroms <- gsub("chrMT","chrM", chroms)

            ## BED has three required columns and up to nine optional
            ## five columns: 1.chrom 2.chromStart 3.chromEnd
            ## 4. OPTIONAL name (displayed next to element)
            ## 5. OPTIONAL score (not used)
            ## 6. OPTIONAL strand (either '+' or '-')
            
            ## bed coordinates are zero-based, half-open
            dat <- data.frame(chr=chroms,
                              chromStart=sprintf("%.0f", object[,1]
                                - ifelse(object@closed[,1], 1L, 0L)),
                              chromEnd=sprintf("%.0f", object[,2]
                                - ifelse(object@closed[,2], 0L, 1L)),
                              stringsAsFactors=FALSE)
            ## names for genomic intervals provided?
            if (!is.null(nameColumn))
              dat$name <- as.character(annotation(object)[[nameColumn]])
            ## strand information provided (Genome_intervals_stranded)?
            if (inherits(object, "Genome_intervals_stranded")){
              ## strand is only 6th column, fill in 'name' and 'score' before
              if (is.null(nameColumn))
                dat$name <- rep(".", nrow(object))
              dat$score  <- integer(nrow(object))
              dat$strand <- as.character(strand(object))
            }# if Genome_intervals_stranded
            resultFile <- paste(gsub("\\..+$","", con),"bed", sep=".")
            girafe:::writeExportData(dat, attribs, "bed", resultFile,
                            writeHeader=writeHeader, append=append)
            invisible(dat)
}) # setMethod("export", signature("Genome_intervals", "character"))
