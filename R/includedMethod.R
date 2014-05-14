
#### 5. interval_included methods for AlignedGenomeIntervals

## 5b. stranded case
setMethod("interval_included",
          signature("AlignedGenomeIntervals", "Genome_intervals_stranded" ),
          function(from, to, check_valid = TRUE){
            nextMethod <- getMethod("interval_included", signature("AlignedGenomeIntervals", "AlignedGenomeIntervals"))
            return(nextMethod(from, to))
          }
) # setMethod("interval_included")

setMethod("interval_included",
          signature("Genome_intervals_stranded" ,"AlignedGenomeIntervals"),
          function(from, to, check_valid = TRUE){
            nextMethod <- getMethod("interval_included", signature("AlignedGenomeIntervals", "AlignedGenomeIntervals"))
            return(nextMethod(from, to))
          }
) # setMethod("interval_included")

## 5c. AlignedGenomeIntervals, AlignedGenomeIntervals
setMethod("interval_included",
          signature("AlignedGenomeIntervals", "AlignedGenomeIntervals" ),
          function(from, to, check_valid = TRUE){
            if ( check_valid && !( validObject(to) && validObject(from) ) )
              stop( "The 'to' and/or 'from' Genome_intervals are invalid." )
            if( any( is.na(strand(to)) ) || any( is.na(strand(from)) ) )
              stop("NA(s) present in the strand of 'to' or 'from'.")
   
            ## loop over strands and call next method
            fromLev <- paste(chromosome(from),strand(from), sep="")
            toLev   <- paste(chromosome(to),strand(to), sep="")
            strdlev <- levels(factor(fromLev))
            
            ## which function to use for each iteration:
            if ("package:parallel" %in% search()){
              lFun <- mclapply
            } else {
              lFun <- lapply
            }
            ### first generate a list of overlap lists, one for each strdlev
            getOneStrandOverlap <- function(s){
               fi <- which(fromLev == s)
               ti <- which(toLev   == s)
               res <-  interval_included(
                         genomeIntervals:::intervalsForOverlap(from[fi]),
                         genomeIntervals:::intervalsForOverlap(to[ti]) )
              ## replace by correct indices:
              res <- lapply(res, function(j) ti[j])
              names(res) <- as.character(fi)
              return(res)  } # getOneStrandOverlap

            perStrand <- lFun(as.list(strdlev), getOneStrandOverlap)

            ### now prepare result:
            rv <- vector( mode="list", length=nrow(from) )
            ## fill in result:
            for (i in 1:length( perStrand)){
              oz <- perStrand[[i]] ## oz is a list
              sIndex <- as.integer(names(oz))
              rv[sIndex] <- oz
            }
            return(rv)
          }
) # setMethod("interval_included")
