### New S4 class for representing reads aligned to genome intervals
###  the idea is to use the functionality, especially the overlap function,
##   of the package genomeIntervals, and extending that class by
##   by additional features useful for aligned reads
##   class definition
setClass("AlignedGenomeIntervals",
         contains = "Genome_intervals",
         representation = representation(
           "sequence"     ="character",
           "reads"        ="integer",
           "matches"      ="integer",
           "organism"     ="character",
           "chrlengths"   ="integer", # to replace 'organism'
           "score"        ="numeric",
           "id"           ="character"),
         prototype = prototype(score=as.numeric(NA)),
         validity = function( object ) {
            fails <- character(0)
            if ( !('strand' %in% names(object@annotation) || !is.factor( object@annotation$strand) || nlevels(object@annotation$strand)!=2)){
              fails <- c(fails, "The 'annotation' slot should have column named 'strand' that is a factor with exactly two levels." )}
            if (length(object@reads)!=nrow(object) ||
                length(object@matches)!=nrow(object) ||
                length(object@sequence)!=nrow(object)){
              fails <- c(fails, "The vectors 'read_sequence','reads' and 'matches' must have exactly the same length as the number of intervals!\n") }
            if (length(object@score)!=1){
              if (length(object@score)!=nrow(object)){
                fails <- c(fails, "The vector 'score' must have exactly the same length as the number of intervals!\n") }}
            if (length(object@id)>0){
              if (length(object@id)!=nrow(object))
                fails <- c(fails, "The vector of interval identifiers 'id' must have exactly the same length as the number of intervals!\n") }
            if ("chrlengths" %in% slotNames(object) &&
                length(object@chrlengths)>0){
              unichrx <- unique(chromosome(object))
              areIn <- unichrx %in% names(object@chrlengths)
              if (!all(areIn))
                fails <- c(fails, paste("Chromosomes '", paste(unichrx[!areIn], collapse=","), "' mentioned in object but not found in names of", "vector of chromosome lengths.")) }
            if (length(fails) > 0) return(fails)
            return(TRUE)
          }
)#setClass("AlignedGenomeIntervals")

## NEW SLOTS OVER genomeIntervals class:
## sequence: sequence of the read aligned to the interval
##           (could think of using  a DNAString instead)
## reads: total number of reads that were aligned to this interval
## matches: the total number of genomic intervals that reads which were
##          aligned to this interval were aligned to. A value of '1' thus
##          means that this read sequence matches uniquely to this one
##          genome interval only
## id: optional identifier for each aligned genome interval
## chrlengths: optional vector of chromosome lengths; used by 'coverage'
##             and other methods in place of retrieving the chromosome
##             lengths from an organism annotation package

### alternative function to create objects of the class
AlignedGenomeIntervals <- function(start, end, chromosome, strand, reads, matches, sequence, id=character(), chrlengths=integer())
{
  stopifnot(length(start)==length(end),
            length(end)==length(chromosome),
            length(chromosome)==length(strand),
            length(strand)==length(sequence),
            is.character(sequence), is.character(id))
  if (missing(reads))
    reads <- rep.int(1L, length(start))
  if (missing(matches))
    matches <- rep.int(1L, length(start))
  strand <- factor(strand, levels=c("+","-"))
  if (any(is.na(strand)))
    stop("Argument 'strand' should consist of '+' and '-' entries only.\n")
  start <- as.integer(start)
  end <- as.integer(end)
  stopifnot(all(end >= start))
  if (length(id)==0)
    id <- character(length(start))
  else
    stopifnot(length(id)==length(start))
  GI <- new("AlignedGenomeIntervals",
            .Data = cbind(start, end),
            annotation=data.frame(
              "seq_name"      = factor(chromosome),
              "strand"        = strand,
              "inter_base"    = vector("logical", length(start))),
            sequence   = sequence,
            reads      = as.integer(reads),
            matches    = as.integer(matches),
            id         = as.character(id),
            chrlengths = as.integer(chrlengths))
  ## test new object:
  stopifnot(validObject(GI))
  return(GI)
}# alternative initiator function for AlignedGenomeIntervals

getReadPosDf <- function(from){
  stopifnot(inherits(from, "AlignedRead"))
  readSeq <- as.character(sread(from))
  readPos <- paste(chromosome(from), from@strand,
                   position(from), position(from)+width(from)-1,
                   readSeq, sep=".")
  tablePos <- table(readPos)
  readDat <- data.frame(chr=chromosome(from),
                        strand=from@strand,
                        start=position(from),
                        end=position(from)+width(from)-1L,
                        seq=readSeq,
                        posfreq=tablePos[match(readPos, names(tablePos))],
                        stringsAsFactors=FALSE)
  ### condense multiply mentioned alignments; element 'posfreq' preserves
  #### the information about number of matching reads
  readDat <- unique(readDat)
  tabSeq <- table(readDat$seq)
  readDat$matches <- tabSeq[match(readDat$seq, names(tabSeq))]
  readDat <- readDat[order(readDat$chr, readDat$start, readDat$end),]
  return(readDat)
}#getReadPosDf

### coercion from ShortRead's class AlignedRead
setAs("AlignedRead", "AlignedGenomeIntervals",
      function(from, to){
        # first drop unaligned reads from the AlignedRead object
        from <- from[!is.na(position(from))]
        readDat <- getReadPosDf(from=from)
        GI <- new("AlignedGenomeIntervals",
                  .Data = cbind(readDat$start, readDat$end),
                  annotation=data.frame(
                    "seq_name"      = factor(readDat$chr),
                    "strand"        = factor(readDat$strand),
                    "inter_base"    = vector("logical", nrow(readDat))),
                  sequence = as.character(readDat$seq),
                  reads    = as.integer(readDat$posfreq),
                  matches  = as.integer(readDat$matches),
                  id       = character(nrow(readDat)))
        return(GI)
})# setAs("AlignedRead", "AlignedGenomeIntervals")

### show method
setMethod("show",signature="AlignedGenomeIntervals",
          function(object){
            cat(formatC(nrow(object), big.mark=","),
                "genome intervals with",
                formatC(sum(object@reads), big.mark=","),
                "aligned reads ")
            cat("on",nlevels(annotation(object)$"seq_name"),"chromosome(s).")
            if (length(object@organism)>0)
              cat("\nOrganism:",object@organism, collapse=" ")
            cat("\n")
            invisible(NULL)}
) # setMethod("show",signature="AlignedGenomeIntervals")

### subsetting:
setMethod("[", signature( "AlignedGenomeIntervals" ),
          function( x, i, j, ..., drop=TRUE) {
            if ( missing(i) ) i <- rep( TRUE, nrow(x) )
            if ( missing(j) ) {
                # Note that both [i,] and [i] syntax subset rows.
                if (is.character(i)){
                  i <- match( i, rownames( x ) )}
                x@sequence   <- x@sequence[i]
                x@matches    <- x@matches[i]
                x@reads      <- x@reads[i]
                # also subset scores of intervals (if specified)
                if (length(x@score)==nrow(x))
                  x@score <- x@score[i]
                if (length(x@id)==nrow(x))
                  x@id <- x@id[i]
            }
            ## slots annotation and .Data are subset by methods for
            ###  super-classes
            callNextMethod( x, i, j, ..., drop=drop)
        }
)#setMethod("[", signature( "AlignedGenomeIntervals" )


###-----------------------------------------------------------------------
### visualization: ( a first attempt)
###----------------------------------------------------------------------
setMethod("plot", signature=c("AlignedGenomeIntervals", "missing"),
          function(x, y, ...){
            plotAligned(x, y=NULL, ...)
          }
) # setMethod("plot", signature=c("AlignedGenomeIntervals","missing"))


### see file plotAligned.R for the source code of the plotting functions
setMethod("plot", signature=c("AlignedGenomeIntervals", "Genome_intervals_stranded"), function(x, y, nameColumn="ID", typeColumn="type", ...)
          {
            stopifnot(nameColumn %in% names(annotation(y)))
            stopifnot(typeColumn %in% names(annotation(y)))
            this.gff <- cbind(annotation(y), start=y[,1], end=y[,2])
            #                 getGffAttribute(y, idAttributeName))
            #names(this.gff)[c(1,8)] <- c("seq_id", "attributes")
            plotAligned(x, y=NULL, gff=this.gff,
                        gffChrColumn="seq_name",
                        gffNameColumn=nameColumn,
                        gffTypeColumn=typeColumn, ...)
          }
) # setMethod("plot", signature=c("AlignedGenomeIntervals","genome_intervals_stranded"))

### accessing subsets by rows of the annotation data.frame:
setMethod("subset", signature(x="AlignedGenomeIntervals"),
          function(x, subset, drop=TRUE, ...){
            z <- eval(substitute(subset), annotation(x), parent.frame())
            x[z,,drop=drop]
})


### accessing the score:
setMethod("score", signature(x="AlignedGenomeIntervals"),
          function(x, ...){
            x@score  })

          
setReplaceMethod("score",
  signature(x="AlignedGenomeIntervals", value="numeric"),
  function(x, value){
    if (length(value)!=nrow(x))
      stop(paste("Number of scores must correspond to the number",
                 "of intervals (",nrow(x),").\n"))
    x@score <- value
    return(x) })

### accessing the strand:
setMethod("strand", signature("AlignedGenomeIntervals"),
          function( x ) x@annotation$strand )

### get the identifiers of aligned genome intervals
setMethod("id", signature("AlignedGenomeIntervals"),
          function(object,...) object@id )

setReplaceMethod("id", signature("AlignedGenomeIntervals","character"),
   function(object, value) {
      stopifnot(length(value)==nrow(object))
      object@id <- value
      return(object) } )

setReplaceMethod("strand", signature("AlignedGenomeIntervals","vector"),
   function(x, value) {
     value <- factor(value)
     strand(x) <- value
})

setReplaceMethod("strand", signature("AlignedGenomeIntervals","factor"),
   function(x, value) {
     if(nlevels(value)!=2) 
       stop(paste("The 'strand' argument should be a vector with",
                  "two distinct values or a factor with two levels.\n"))
     if(!( length( value ) %in% c(1,nrow(x))))
       stop(paste("The 'strand' argument should be a vector or a",
                  "factor of length equal to 1 or to the number of",
                  "rows of the end point matrix.\n"))
     if(length(value)==1) value <- rep(value, nrow(x))
     x@annotation$strand <- value
     return(x)
})

### organism assignment
setMethod("organism", signature("AlignedGenomeIntervals"),
          function( x ) x@organism
)

setReplaceMethod("organism", signature("AlignedGenomeIntervals","character"),
   function(x, value) {
     stopifnot(length(value)==1)
     orgpackage <- paste("org",value,"eg.db", sep=".")
     suppressWarnings(packageThere <- require(orgpackage, character.only=TRUE))
     if (!packageThere){
       warning(paste("An annotation package called '",orgpackage,"' ",
                     "could not be found.\n","Please check whether '",
                     value, "' is a valid organism name ",
                     "or install the package. Example organism identifiers ",
                     "are 'Hs' for Human and 'Mm' for Mouse. ",
                     "The annotation  packages can usually be obtained from ",
                     "the Bioconductor repositories. You may want to try:\n",
                     "> source('http://www.bioconductor.org/biocLite.R')\n",
                     "> biocLite('",orgpackage,"')\n", sep=""))
     }
     x@organism <- value
     x
})

### accessor and replacement methods for vector of chromosome lengths
setMethod("chrlengths", signature("AlignedGenomeIntervals"),
          function( x ) x@chrlengths
)

setReplaceMethod("chrlengths", signature("AlignedGenomeIntervals","numeric"),
   function(x, value) {
     if (is.null(names(value)))
       stop("Vector of chromosome lengths must have names!")
     unichrx <- unique(chromosome(x))
     areIn <- unichrx %in% names(value)
     if (!all(areIn))
       stop("Chromosomes '", paste(unichrx[!areIn], collapse=","),
            "' mentioned in object but not found in names of ",
            "vector of chromosome lengths!")
     mode(value) <- "integer"
     x@chrlengths <- value
     x
})

### get chromosome annotation
setMethod("chromosome", signature("Genome_intervals"),
          function(object, ...)
          as.character(object@annotation$"seq_name"))

setMethod("chromosome", signature("AlignedGenomeIntervals"),
          function(object, ...)
          callNextMethod(object, ...))

setMethod("seq_name", signature("AlignedGenomeIntervals"),
          function(x)
          callNextMethod(x))

### Other accessor methods
setMethod("matches", signature("AlignedGenomeIntervals"),
          function( x ) x@matches
)

setReplaceMethod("matches", signature("AlignedGenomeIntervals","integer"),
   function(x, value) {
     stopifnot(length(value)==nrow(x))
     x@matches <- value
     x
})

setMethod("reads", signature("AlignedGenomeIntervals"),
          function( x ) x@reads
)

setReplaceMethod("reads", signature("AlignedGenomeIntervals","character"),
   function(x, value) {
     stopifnot(length(value)==nrow(x))
     x@reads <- value
     x
})

### lenght of aligned genome intervals:
setMethod("width", signature(x="AlignedGenomeIntervals"),
          function(x) as.integer(x[,2]-x[,1]+1L))

### length of sequences:
setMethod("nchar", signature(x="AlignedGenomeIntervals"),
          function(x, type = "chars", allowNA = FALSE)
          nchar(x@reads))

### detailed information:
setMethod("detail", signature(x="AlignedGenomeIntervals"),
          function(x, ...){
            res <- cbind("start"=x@.Data[,1],
                         "end"=x@.Data[,2],
                         x@annotation,
                         reads=x@reads,
                         matches=x@matches,
                         sequence=x@sequence)
            if (any(nchar(x@id)>0))
              res$id <- x@id
            rownames(res) <- NULL
            res$"inter_base" <- NULL
            res
}) #detail
            
### export methods:
# see file export.R

### combining AlignedGenomeIntervals objects using base function 'c':
c.AlignedGenomeIntervals <- function( ... )
{
  args <- list(...)
  
  ## Drop NULL arguments
  if ( any( sapply( args, is.null ) ) )
    args <- args[ !sapply( args, is.null ) ]
  
  ## If anything else than AlignedGenomeIntervals returns a list
  classes <- sapply( args, class )
  if ( !all( classes %in% c( "AlignedGenomeIntervals")))
    return( list( ... ) )

  ## total number of intervals in result:
  ntotal <- sum(sapply(args, nrow))

  ## organism annotation
  orgs <- do.call("c", lapply(args, function(x) x@organism))
  if (length(unique(orgs))>1)
    stop("Cannot combine intervals aligned to different genomes.")

  ## chromosome lengths vector
  comb.chrlengths <- do.call("c", lapply(args, function(x) x@chrlengths))
  if (length(comb.chrlengths)>0 &&
      any(duplicated(names(comb.chrlengths)))){
    splitted <- split(as.vector(comb.chrlengths),
                      names(comb.chrlengths))
    comb.chrlengths <- sapply(splitted, max)
  }
 
  # rbinds the data frame if possible
  annot <-  try( do.call( rbind, lapply( args, function(x) annotation(x) ) ), silent = TRUE )
  if( class(annot) == "try-error")
    stop("failed to combine the annotation slots by 'rbind':\n", annot)

  ### combine scores:
  scores <- do.call("c", lapply(args, function(x) x@score))
  if (length(scores) != ntotal) ## in case of ununused scores
    scores <- as.numeric(NA)

  # append the rest
  result <- new("AlignedGenomeIntervals",
                do.call("c", lapply(args, as, "Intervals_full")),
                annotation = annot,
                sequence=do.call("c", lapply(args, function(x) x@sequence)),
                reads=do.call("c", lapply( args, function(x) x@reads)),
                matches=do.call("c", lapply(args, function(x) x@matches)),
                organism=unique(orgs),
                score=scores,
                id=do.call("c", lapply(args, function(x) x@id )),
                chrlengths=comb.chrlengths
                )
  return(result)
}# c.AlignedGenomeIntervals


### coercion from 'AlignedGenomeIntervals' class to
###  class 'RangedData' (package 'IRanges')
setAs("AlignedGenomeIntervals", "RangedData",
      function(from, to){
        coords <- IRanges(start=from[,1], end=from[,2])
        score  <- from@score
        org <- from@organism
        if (length(org)==0) org <- NULL
        if (any(from$"inter_base"))
          warning("Inter-base intervals contained in '",
                  deparse(substitute(from)),
                  "' cannot be properly converted.")
        RD <- RangedData(coords,
                        strand = as.character(strand(from)),
                        reads = from@reads,
                        matches = from@matches,
                        sequence = from@sequence,
                        score=score,
                        space = as.character(seq_name(from)),
                        universe = org[1]) # must be of length 1
        if (!validObject(RD))
          stop("Conversion failed.")
        return(RD)
})# setAs("AlignedGenomeIntervals", "RangedData",


### read extension, for example for ChIP-Seq data:
setMethod("extend", signature(x="AlignedGenomeIntervals"),
          function(x, threePrime=0L, fivePrime=0L, ...){
            ## left end:
            x[,1] <- pmax(1, x[,1] - ifelse(strand(x)=="-",
                                            threePrime, fivePrime))
            ## right end:
            x[,2] <- x[,2] + ifelse(strand(x)=="-",
                                    fivePrime, threePrime)
            ## do not exceed chromlength:
            chrlens <- getChromLengths(x)
            xchrom <- gsub("MT$", "M", chromosome(x))
            x[,2] <- pmin.int(x[,2], chrlens[xchrom], na.rm=TRUE)
            return(x)
}) #extend

### 'clusters' method
setMethod("clusters", signature("AlignedGenomeIntervals"),
          function (x, w=1, check_valid=TRUE, which=TRUE, ...){
            if ( check_valid && !( validObject(x)) )
              stop("The 'AlignedGenomeIntervals' object is invalid.")
            
            ## loop over chromosomes and call next method
            xLev <- factor(paste(chromosome(x),strand(x), sep=""))
            strdlev <- levels(xLev)
            
            ## which function to use for each iteration:
            if ("package:multicore" %in% search()){
              lFun <- mclapply
            } else {
              lFun <- lapply
            }
            ### first generate a list of overlap lists, one for each strdlev
            getOneStrandClusters <- function(s, w){
              xi <- which(xLev == s)
              res <-  intervals:::clusters(
                         genomeIntervals:::intervalsForOverlap(x[xi]),
                         w = w , which=TRUE)
               ## replace by correct indices:
              res <- lapply(res, function(j) xi[j])
              return(res)  } # getOneStrandClusters

            perStrand <- lFun(as.list(strdlev), getOneStrandClusters, w=w)
            rv <- do.call("c", perStrand)
            
            if (! which) ## do we want the actual intervals?
              rv <- lapply(rv, function(j) x[j])
            return(rv)
          }
) # setMethod("clusters", signature("AlignedGenomeIntervals"))

### histogram of widths of aligned reads
setMethod("hist", signature(x="AlignedGenomeIntervals"),
          function(x, plot=TRUE, ...){
            splitted <- split(reads(x), width(x))
            icounts <- sapply(splitted, sum)
            imids   <- as.integer(names(icounts))
            ## fill in missing lengths for homogenous histogram:
            missingLen <- setdiff(min(imids):max(imids), imids)
            imids <- c(imids, missingLen)
            icounts <- c(icounts, integer(length(missingLen)))
            ord <- order(imids)
            imids <- imids[ord]
            icounts <- icounts[ord]
            ibreaks <- c(imids-0.5, imids[length(imids)]+0.5)
            h <- list(breaks=ibreaks, counts=icounts,
                      mids=imids)
            class(h) <- c("histogram", class(h))
            if (plot) plot(h, ...)
            invisible(h) }
) # setMethod("hist", signature("AlignedGenomeIntervals"))

### sort intervals
setMethod("sort", signature(x="AlignedGenomeIntervals"),
          function(x, decreasing=FALSE, ...) {
            chr <- gsub("^[Cc]hr", "", chromosome(x))
            ## replace some chromosome names by numbers for sorting:
            chr <- gsub("X$","100", chr)
            chr <- gsub("Y$","200", chr)
            chr <- gsub("MT?$","300", chr)
            chr <- gsub("_random$","000", chr)
            suppressWarnings(chr <- as.numeric(chr))
            ## are there still non-numeric entries left in chr?
            if (any(is.na(chr))) chr <- chromosome(x)
            ord <- order(chr, x[,1], x[,2], decreasing=decreasing)
            return(x[ord])
} ) # sort

## subsample; draw n alignd reads from an AlignedGenomeIntervals
##  object and return the resulting AlignedGenomeIntervals object
setMethod("sample", signature(x="AlignedGenomeIntervals"),
 function(x, size, replace=TRUE, ...){
  stopifnot(inherits(x, "AlignedGenomeIntervals"),
            is.numeric(size), is.logical(replace))
  ## only supports sampling with replacement at the moment.
  #if (!replace && size>sum(reads(x)))
  #  stop("cannot take a sample larger than the population when 'replace = FALSE'\n")
  ## OLD VERSION: use Rle object during subsampling
  # -> but RLE problems with large objects "too large a range of values in 'x'"
  #R <- Rle(values=seq_len(nrow(x)), length=reads(x))
  #R2 <- sort(sample(R, size=size, replace=replace, ...))
  #x2 <- x[runValue(R2)]
  #x2@reads <- runLength(R2)

  ## NEW VERSION: using sample.int
  sel <- sample.int(nrow(x), size=size,
                    replace=TRUE, prob=reads(x))
  sel <- sort(sel)
  freq <- tabulate(sel, nbins=nrow(x))
  ## remove zero-read intervals
  freq <- freq[freq > 0L]
  x2 <- x[unique(sel)]
  x2@reads <- freq
  stopifnot(validObject(x2))
  return(x2)
} ) #sample


## summary method, currently display number of reads per chromosome,
##  optional limit to reads with specific number of matches
setMethod("summary", signature(object="AlignedGenomeIntervals"),
 function(object, nMatches=NA, ...){
   output <- "Number of aligned reads"
   if (!is.na(nMatches)){
     stopifnot(is.numeric(nMatches), length(nMatches)==1L)
     object <- object[matches(object)==nMatches]
     output <- paste(output, "(with", nMatches, "matches)")
   }
   splitted  <- split(reads(object), chromosome(object))
   chrcounts <- sapply(splitted, sum)
   output <- paste(output, "per chromosome:\n")
   cat(output)
   print(chrcounts)
   invisible(chrcounts)
} ) # summary
