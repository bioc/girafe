### collapse/reduce aligned genome intervals by
## 1: combining intervals which are completely included in each other
## 2: combining overlapping intervals
setMethod("reduce", signature("AlignedGenomeIntervals"),
          function(x, exact=FALSE, ...){
            stopifnot(is.logical(exact))

            ## which function to use for each iteration:
            if ("package:multicore" %in% search()){
              lFun <- mclapply
            } else {
              lFun <- lapply
            }
            ## separate method for reducing only intervals
            ###  at exactly the same position?
            if (exact) {
              ## check: use of new 'interval_included' method may
              ##  be faster here.
              ## basically generate another list 'ov' here
              ov <- as.list(1:nrow(x))
              fo <- fracOverlap(x,x, 1.0)
              fo <- subset(fo, frac1 == 1.0 & frac2 == 1.0)
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
            ### intervals should only be combined if they are
            # overlapping/adjacent AND have the sampe number of matches
            # (i.e. the same match specificity [unique, 2 hits, 5 hits])
            ov2 <- lFun(as.list(seq(length(ov))),
                   function(i){
                     ovI <- ov[[i]]
                     if (length(ovI)<2) return(ovI)
                     return(ovI[x@matches[ovI] == x@matches[i]])
            })
            if (!any(listLen(ov2)>1)) return(x)

            ### prepare result of class AlignedGenomeIntervals:
            xr <- x # [1:length(CC)]

            ### OLD STUFF:
            #treated <- vector("logical", length(ov))
            treated <- rep(TRUE, length(ov))
            removed <- vector("logical", length(ov))
            ofInt <- which(listLen(ov2)>1)
            treated[ofInt] <- FALSE
            done <- FALSE
            #if there are overlaps:
            #for (i in which(listLen(ov2)>1)){
            #  if (treated[i]) next
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
                minx <- min(x[theseIdx,1])
                maxx <- max(x[theseIdx,2])
                xr@closed[i,] <- c(x@closed[which.min(x[theseIdx,1]),1],
                                   x@closed[which.max(x[theseIdx,2]),2])
                xr@reads[i]   <- sum(x@reads[theseIdx])
                ## now for the sequence use 'shift' arg of consensusString
                # shift depends on strand of the reads:
                if (strand(x[i])=="-") {
                  theseShifts <- maxx-x[theseIdx,2]
                } else { theseShifts <- x[theseIdx,1]-minx }
                xr@sequence[i] <-
                  consensusString(DNAStringSet(x@sequence[theseIdx]),
                                  ambiguityMap="N", #still problem with IUPAC 
                                  shift=theseShifts)
                xr[i,1] <- minx
                xr[i,2] <- maxx
              } # else
            }# while (!done)
            xr <- xr[!removed]
            stopifnot(validObject(xr))
            return(xr)
}) # setMethod reduce


### for Genome_intervals
setMethod("reduce", signature("Genome_intervals"),
          function(x, exact=FALSE, ...){
            stopifnot(is.logical(exact))
            ## separate method for reducing only intervals
            ###  at exactly the same position?
            ## which function to use for each iteration:
            if ("package:multicore" %in% search()){
              lFun <- mclapply
            } else {
              lFun <- lapply
            }
            
            if (exact) {
            ## basically generate another list 'ov' here
              ov <- as.list(1:nrow(x))
              fo <- fracOverlap(x,x, 1.0)
              fo <- subset(fo, frac1 == 1.0 & frac2 == 1.0)
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
