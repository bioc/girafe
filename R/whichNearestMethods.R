## method which_nearest for Genome_intervals
setMethod("which_nearest",
          signature("Genome_intervals", "Genome_intervals"),
          function( from, to ) {
            ## return value
            n <- nrow(from)
            res <- data.frame(distance_to_nearest=rep(as.numeric(NA),n),
                              which_nearest=I(vector("list", n)),
                              which_overlap=I(vector("list", n)))
            ## adapt interval representation (basically moves to R
            ##  with appropriate changes for inter_base)
            ## empty intervals become NA
            fints <- genomeIntervals:::intervalsForDistance(from)
            tints <- genomeIntervals:::intervalsForDistance(to)

            ## loop over seqnames and call next method
            seqfrom <- as.character(seqnames(from))
            seqto   <- as.character(seqnames(to))

            ## loop over seqnames and call next method
            for(s in unique(seqfrom)){
                fi <- which(seqfrom == s)
                ti <- which(seqto   == s)
                if (length(ti)==0){
                  warning(paste("No intervals on chromosome",s,
                                "in genome_intervals object",
                                deparse(substitute(to)),"."))
                  next
                }
                W  <- which_nearest(fints[fi],
                                    tints[ti])
                res$"distance_to_nearest"[fi] <- W$"distance_to_nearest"
                res$"which_nearest"[fi] <- lapply(W$"which_nearest",
                                                  function(x) ti[x])
                res$"which_overlap"[fi] <- lapply(W$"which_overlap",
                                                  function(x) ti[x])
            }#for(s in seqlev)
            return(res)
        }
)# setMethod("which_nearest")

## if both are stranded, then take strands in account
## user must cast to unstranded if this is not wished.
setMethod("which_nearest",
          signature("Genome_intervals_stranded",
                    "Genome_intervals_stranded"),
          function( from, to ) {
            if( any( is.na(strand(to)) ) )
              stop("NA(s) present in the strand of 'to'.")
            ## return value
            n <- nrow(from)
            res <- data.frame(distance_to_nearest=rep(as.numeric(NA),n),
                              which_nearest=I(vector("list", n)),
                              which_overlap=I(vector("list", n)))
            strandFrom <- as.character(strand(from))
            strandTo   <- as.character(strand(to))
            
            strdlev <- unique(strandFrom)
            ## loop over strands and call next method
            nextMethod <- getMethod("which_nearest", signature("Genome_intervals", "Genome_intervals"))
            for(s in strdlev){
              si <- which(strandFrom == s)
              ti <- which(strandTo   == s)
              wni <- nextMethod( from[si], to[strand(to) == s] )
              res$"distance_to_nearest"[si] <- wni$"distance_to_nearest"
              res$"which_nearest"[si] <- lapply(wni$"which_nearest",
                                                function(x) ti[x])
              res$"which_overlap"[si] <- lapply(wni$"which_overlap",
                                                function(x) ti[x])
            }
            return(res)
          }
)# setMethod("which_nearest"), signature("Genome_intervals_stranded", "Genome_intervals_stranded"))

## for AlignedGenomeIntervals
setMethod("which_nearest",
          signature("AlignedGenomeIntervals", "Genome_intervals_stranded"),
          function( from, to ) {
            nextMethod <- getMethod("which_nearest", signature( "Genome_intervals_stranded",  "Genome_intervals_stranded" ))
            nextMethod(from, to)
          }
)# setMethod("which_nearest", signature("AlignedGenomeIntervals", "Genome_intervals_stranded"))
