
#### 5. interval_overlap methods for AlignedGenomeIntervals

## 5a. un-stranded case
setMethod("interval_overlap",
          signature("AlignedGenomeIntervals", "Genome_intervals"),
          function(from, to, check_valid = TRUE){
            if ( check_valid && !( validObject(to) && validObject(from) ) )
              stop( "The 'to' and/or 'from' Genome_intervals are invalid." )
            nextMethod <- getMethod("interval_overlap", signature("Genome_intervals", "Genome_intervals"))
            return(nextMethod(from, to))}
) # setMethod("interval_overlap")

setMethod("interval_overlap",
          signature("Genome_intervals", "AlignedGenomeIntervals"),
          function(from, to, check_valid = TRUE){
            if ( check_valid && !( validObject(to) && validObject(from) ) )
              stop( "The 'to' and/or 'from' Genome_intervals are invalid." )
            nextMethod <- getMethod("interval_overlap", signature("Genome_intervals", "Genome_intervals"))
            return(nextMethod(from, to))}
) # setMethod("interval_overlap")

## 5b. stranded case
setMethod("interval_overlap",
          signature("AlignedGenomeIntervals", "Genome_intervals_stranded" ),
          function(from, to, check_valid = TRUE){
            nextMethod <- getMethod("interval_overlap", signature("AlignedGenomeIntervals", "AlignedGenomeIntervals"))
            return(nextMethod(from, to))
          }
) # setMethod("interval_overlap")

setMethod("interval_overlap",
          signature("Genome_intervals_stranded" ,"AlignedGenomeIntervals"),
          function(from, to, check_valid = TRUE){
            nextMethod <- getMethod("interval_overlap", signature("AlignedGenomeIntervals", "AlignedGenomeIntervals"))
            return(nextMethod(from, to))
          }
) # setMethod("interval_overlap")

## 5c. AlignedGenomeIntervals, AlignedGenomeIntervals
## outdated, check below for new version
oldAGIoverlap <- function(from, to, check_valid = TRUE){
            if ( check_valid && !( validObject(to) && validObject(from) ) )
              stop( "The 'to' and/or 'from' Genome_intervals are invalid." )
            if( any( is.na(strand(to)) ) || any( is.na(strand(from)) ) )
              stop("NA(s) present in the strand of 'to' or 'from'.")

   
            ## loop over strands and call next method
            fromLev <- paste(chromosome(from),strand(from), sep="")
            toLev   <- paste(chromosome(to),strand(to), sep="")
            strdlev <- levels(factor(fromLev))
            
            ## which function to use for each iteration:
            if ("package:multicore" %in% search()){
              lFun <- mclapply
            } else {
              lFun <- lapply
            }
            ### first generate a list of overlap lists, one for each strdlev
            getOneStrandOverlap <- function(s){
               fi <- which(fromLev == s)
               ti <- which(toLev   == s)
               res <-  interval_overlap(
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
}  # oldAGIoverlap


### try new method without lapply:
## add to every position a scaling factor larger than the maximal
##  coordinate, negative for minus strand
setMethod("interval_overlap",
          signature("AlignedGenomeIntervals", "AlignedGenomeIntervals" ),
          function(from, to, check_valid = TRUE){
            if ( check_valid && !( validObject(to) && validObject(from) ) )
              stop( "The 'to' and/or 'from' Genome_intervals are invalid." )
            if( any( is.na(strand(to)) ) || any( is.na(strand(from)) ) )
              stop("NA(s) present in the strand of 'to' or 'from'.")
            ## which chromosomes are there:
            fromChrom <- as.character(chromosome(from))
            toChrom   <- as.character(chromosome(to))
            allChrom <- union(fromChrom, toChrom)
            ## get scaling factor from data:
            sf <- 10^ceiling(log10(max(max(from[,2]),max(to[,2]))))
            ## do we need a check whether this is greater than
            ##  the machine double range?
            stopifnot(sf > 0)
            adjustCoords <- function(x, fact=1e9){
              XI <- genomeIntervals:::intervalsForOverlap(x)
              shifts <- fact * match(chromosome(x), allChrom) *
                ifelse(strand(x)=="-", -1, +1)
              XI <- XI + shifts
              return(XI)
            }
            FI <- adjustCoords(from, sf)
            TI <- adjustCoords(to, sf)
            interval_overlap(FI, TI)
} ) # setMethod

#system.time(I1 <- interval_overlap(exAI, mm.gi))
#system.time(I2 <- oldAGIoverlap(exAI, mm.gi))
