### collapse/reduce aligned genome intervals by
## 1: combining intervals which are completely included in each other
## 2: combining overlapping intervals
setMethod("reduce", signature("AlignedGenomeIntervals"),
          function(x, method="standard", min.frac=0.0,
                   mem.friendly=FALSE, ...){
            ## mem.friendly: version that requires less RAM but is considerably slower
            stopifnot(is.logical(mem.friendly))
            method <- match.arg(method, c("standard", "exact",
                                          "same5", "same3"))
            ## which function to use for each iteration:
            if ("package:multicore" %in% search()){
              lFun <- mclapply
            } else {
              lFun <- lapply
            }
            if (mem.friendly){
              #lFun over individual chromosomes
              allChr <- unique(chromosome(x))
              res <- lFun(as.list(allChr), function(chr){
                reduceOne(x[x@annotation$seq_name==chr],
                          method=method, min.frac=min.frac)
              })
              res <- do.call("c", res)
            } else {
              res <- reduceOne(x, method=method, min.frac=min.frac)
            }
            res <- sort(res)
            stopifnot(validObject(res))
            return(res)
}) # setMethod reduce

## function for single iteration of reduce
reduceOne <- function(z, method=method, min.frac=0.0){
  ## separate method for reducing only intervals
  ###  at exactly the same position?
  method <- match.arg(method, c("standard", "exact",
                                "same5", "same3"))
  ## only reduce intervals with exactly same start and end?
  if (method == "exact") { #min.frac <- 1.0
    return(reduceOneExact(z))
  }
  ## only reduce intervals with exactly same 5' or 3' end?
  if (method %in% c("same5","same3"))
    return(reduceOneEnd(z, type=method))
  
  if (min.frac > 0.0) {
  ## check: use of new 'interval_included' method may
  ##  be faster here.
  ## basically generate another list 'ov' here
    ov <- as.list(1:nrow(z))
    fo <- fracOverlap(z,z, 0)
    fo <- subset(fo, fraction1 >= min.frac | fraction2 >= min.frac)
    perInd <- split(fo$Index2, fo$Index1)
    stopifnot(length(perInd)==nrow(z))
    for (j in 1:length(perInd))
      ov[[as.integer(names(perInd))[j]]] <- perInd[[j]]
  } else {
    # find overlapping and adjacent intervals:
    z2 <- z # for adjacent intervals
    z2[,1] <- z[,1]-1L; z2[,2] <- z[,2]+1L
    ov <- interval_overlap(z,z2)
  }# else
  ### intervals should only be combined if they are
  # overlapping/adjacent AND have the sampe number of matches
  # (i.e. the same match specificity [unique, 2 hits, 5 hits])
  ov2 <- lapply(as.list(seq(length(ov))),
              function(i){
                ovI <- ov[[i]]
                if (length(ovI)<2) return(ovI)
                return(ovI[z@matches[ovI] == z@matches[i]])
              })
  if (!any(listLen(ov2)>1)) return(z)
  ### prepare result of class AlignedGenomeIntervals:
  zr <- z # [1:length(CC)]
  hasIds <- length(z@id)==nrow(z)
  
  #treated <- vector("logical", length(ov))
  treated <- rep(TRUE, length(ov))
  removed <- vector("logical", length(ov))
  ofInt <- which(listLen(ov2)>1)
  treated[ofInt] <- FALSE
  done <- FALSE
  #if there are overlaps:
  while (!done) {
    i <- match(FALSE, treated)
    if (is.na(i)) { done <- TRUE }
    else {
      # repeat loop for single-linkage clustering-like
      ##  combination of consecutive match positions
      theseIdx <- ov2[[i]]
      repeat { ## works, but very slow
        newIdx <- unique(unlist(ov2[theseIdx], use.names=FALSE))
        if (all(newIdx %in% theseIdx)) break
        theseIdx <- union(theseIdx, newIdx)
      }
      treated[theseIdx] <- TRUE
      removed[theseIdx] <- TRUE
      removed[i] <- FALSE
      minx <- min(z[theseIdx,1])
      maxx <- max(z[theseIdx,2])
      zr@closed[i,] <- c(z@closed[which.min(z[theseIdx,1]),1],
                         z@closed[which.max(z[theseIdx,2]),2])
      zr@reads[i]   <- sum(z@reads[theseIdx])
      
      ## now for the sequence use 'shift' arg of consensusString
      # shift depends on strand of the reads:
      if (strand(z[i])=="-") {
        theseShifts <- maxx-z[theseIdx,2]
      } else { theseShifts <- z[theseIdx,1]-minx }
      ## get consensus matrix of read sequences:
      #iConsMat <- consensusMatrix(DNAStringSet(z@sequence[theseIdx]),
      #                            shift=theseShifts, baseOnly=TRUE)
      ## NEW VERSION with consensus weighted by read counts
      iConsMat <- weightedConsensusMatrix(z@sequence[theseIdx],
                                          z@reads[theseIdx],
                                          shift=theseShifts)
      zr@sequence[i] <- consensusString(iConsMat,  ambiguityMap="N")
      if (hasIds)
        zr@id[i] <- paste(sort(unique(z@id[theseIdx])),
                          collapse=",")
      zr[i,1] <- minx
      zr[i,2] <- maxx
    } # else
  }# while (!done)
  zr <- zr[!removed]
  return(zr)
}#reduceOne

reduceOneExact <- function(z){
  # simpler method if all intervals at exactly same position
  hasIds <- length(z@id)==nrow(z)
  readPos <- paste(chromosome(z), strand(z),z[,1], "-", z[,2],
                   matches(z), sep=".")
  splitted <- split(seq.int(nrow(z)), readPos)
  # prepare result: single aligned interval per set of overlapping intervals
  zr <- z[sapply(splitted, "[", 1L)]
  for (i in which(listLen(splitted)>1L)){
    ## iterate over each group of completely overlapping intervals
    theseIdx <- splitted[[i]]
    zr@reads[i] <- sum(z@reads[theseIdx])
    ## now for the sequence use consensusString
    ##  without shift since all reads on same position
    #zr@sequence[i] <-
    #  consensusString(DNAStringSet(z@sequence[theseIdx]),
    #                  ambiguityMap="N")
    ##New version: weighing consensus by read counts:
    iConsMat <- weightedConsensusMatrix(z@sequence[theseIdx],
                                        z@reads[theseIdx])
    zr@sequence[i] <- consensusString(iConsMat,  ambiguityMap="N")
    if (hasIds)
      zr@id[i] <- paste(sort(unique(z@id[theseIdx])),
                        collapse=",")
  }
  return(zr)  
}#reduceOneExact

reduceOneEnd <- function(z, type="same5"){
  # simpler method if all intervals have exactly same 5' or 3' end
  type <- match.arg(type, c("same5", "same3"))
  hasIds <- length(z@id)==nrow(z)
  if (type=="same5")
    readPos <- paste(chromosome(z), strand(z),
                     ifelse(strand(z)=="+", z[,1], z[,2]),
                     matches(z), sep=".")
  if (type=="same3")
    readPos <- paste(chromosome(z), strand(z),
                     ifelse(strand(z)=="+", z[,2], z[,1]),
                     matches(z), sep=".")
  splitted <- split(seq.int(nrow(z)), readPos)
  # prepare result: single aligned interval per set of overlapping intervals
  zr <- z[sapply(splitted, "[", 1L)]
  for (i in which(listLen(splitted)>1L)){
    ## iterate over each group of completely overlapping intervals
    theseIdx <- splitted[[i]]
    zr@reads[i] <- sum(z@reads[theseIdx])
    ## now for the sequence use consensusString
    ##  without shift since all reads on same position
    #zr@sequence[i] <-
    #  consensusString(DNAStringSet(z@sequence[theseIdx]),
    #                  ambiguityMap="N")
    ## NEW: weighted consensus sequence
    iConsMat <- weightedConsensusMatrix(z@sequence[theseIdx],
                                        z@reads[theseIdx])
    zr@sequence[i] <- consensusString(iConsMat,  ambiguityMap="N")
    if (hasIds)
      zr@id[i] <- paste(sort(unique(z@id[theseIdx])),
                        collapse=",")
  }
  return(zr)  
}#reduceOneEnd

### for Genome_intervals
setMethod("reduce", signature("Genome_intervals"),
          function(x, method="standard", min.frac=0.0, ...){
            method <- match.arg(method, c("standard", "exact"))
            ## separate method for reducing only intervals
            ###  at exactly the same position?
            ## which function to use for each iteration:
            if ("package:multicore" %in% search()){
              lFun <- mclapply
            } else {
              lFun <- lapply
            }
            if (method=="exact") min.frac <- 1.0
            if (min.frac > 0.0) {
            ## basically generate another list 'ov' here
              ov <- as.list(1:nrow(x))
              fo <- fracOverlap(x,x, 0)
              fo <- subset(fo, fraction1 >= min.frac | fraction2 >= min.frac)
              perInd <- split(fo$Index2, fo$Index1)
              stopifnot(length(perInd)==nrow(x))
              for (j in 1:length(perInd))
                ov[[as.integer(names(perInd))[j]]] <- perInd[[j]]
            } else {
              # find overlapping and adjacent intervals:
              x2 <- x # for adjacent intervals
              x2[,1] <- x[,1]-1L; x2[,2] <- x[,2]+1L
              ov <- interval_overlap(x,x2)
            }# else
            if (!any(listLen(ov)>1)) return(x)
            ### prepare result of class Genome_intervals:
            xr <- x # [1:length(CC)]
            ## which intervals to keep and which to remove
            treated <- vector("logical", length(ov))
            removed <- vector("logical", length(ov))
            #if there are overlaps:
            for (i in which(listLen(ov)>1)){
              if (treated[i]) next
              # repeat loop for single-linkage clustering-like
              ##  combination of consecutive match positions
              theseIdx <- ov[[i]]
              repeat { ## works, but very slow
                newIdx <- unique(unlist(ov[theseIdx], use.names=FALSE))
                if (all(newIdx %in% theseIdx)) break
                theseIdx <- union(theseIdx, newIdx)
              }
              treated[theseIdx] <- TRUE
              removed[theseIdx] <- TRUE
              removed[i] <- FALSE
              #for (i in 1:length(CC)){
              #theseIdx <- as.integer(gsub("^n", "", CC[[i]]))
              minx <- min(x[theseIdx,1])
              maxx <- max(x[theseIdx,2])
              xr@closed[i,] <- c(x@closed[which.min(x[theseIdx,1]),1],
                                x@closed[which.max(x[theseIdx,2]),2])
              xr[i,1] <- minx
              xr[i,2] <- maxx
              }# for i
            xr <- xr[!removed]
            stopifnot(validObject(xr))
            return(xr)
}) # setMethod reduce for Genome_intervals

### resurrect the reduce methods from package 'IRanges'
###  which were overwritten by later import of 'intervals'
# for class 'IRanges'
setMethod("reduce", signature("IRanges"),
          function(x, ...){
            getMethod("reduce", signature("IRanges"), where=match("package:IRanges",search())) (x, ...)  }  )
# for class 'Ranges'
setMethod("reduce", signature("Ranges"),
          function(x, ...){
            getMethod("reduce", signature("Ranges"), where=match("package:IRanges",search())) (x, ...)  }  )
# for class 'RangesList'
setMethod("reduce", signature("RangesList"),
          function(x, ...){
            getMethod("reduce", signature("RangesList"), where=match("package:IRanges",search())) (x, ...)  }  )
# for class 'CompressedIRangesList'
setMethod("reduce", signature("CompressedIRangesList"),
          function(x, ...){
            getMethod("reduce", signature("CompressedIRangesList"), where=match("package:IRanges",search())) (x, ...)  }  )
# for class 'RangedData'
setMethod("reduce", signature("RangedData"),
          function(x, by, ...){
            if (missing(by)) by <- seq(ncol(rd))
            getMethod("reduce", signature("RangedData"), where=match("package:IRanges",search())) (x, by=by, ...)  }  )

