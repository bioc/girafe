### New methods for classes Genome_intervals and Genome_intervals_stranded

## additional accessor:
setMethod("chromosome", signature("Genome_intervals"),
          function(object, ...)
          as.character(object@annotation$"seq_name"))

### read extension, for example for ChIP-Seq data:
setMethod("extend", signature(x="Genome_intervals_stranded"),
          function(x, threePrime=0L, fivePrime=0L, ...){
            ## left end:
            x[,1] <- pmin.int(1L, x[,1] - ifelse(strand(x)=="-",
                                                 threePrime, fivePrime))
            ## right end:
            x[,2] <- x[,2] + ifelse(strand(x)=="-",
                                    fivePrime, threePrime)
            return(x)
}) #extend

setMethod("extend", signature(x="Genome_intervals"),
          function(x, threePrime=0L, fivePrime=0L, ...){
            ## left end:
            x[,1] <- pmin.int(1L, x[,1] - fivePrime )
            ## right end:
            x[,2] <- x[,2] + threePrime
            return(x)
}) #extend

### fracOverlap:
## see file fracOverlap.R

### reduce
## see file reduceMethod.R

### 'clusters' method for Genome_intervals
setMethod("clusters", signature("Genome_intervals"),
          function (x, w=1, check_valid=TRUE, which=TRUE, ...){
            if ( check_valid && !( validObject(x)) )
              stop("The 'Genome_intervals' object is invalid.")
            
            ## loop over chromosomes and call next method
            #xLev <- paste(chromosome(from),strand(from), sep="")
            xLev <- factor(chromosome(x))
            strdlev <- levels(xLev)
            
            ## which function to use for each iteration:
            if ("package:parallel" %in% search()){
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
) # setMethod("clusters", signature("Genome_intervals"))

