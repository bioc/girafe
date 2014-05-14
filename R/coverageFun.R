
### generic method definition from IRanges
#setGeneric("coverage", signature="x",
#    function(x, shift=0L, width=NULL, weight=1L, ...)
#        standardGeneric("coverage")
#)

# R function for calling C++ coverage method:
coverageOneStrand <- function(x, chrmax=NA){
  stopifnot(inherits(x, "AlignedGenomeIntervals"))
  if (nrow(x) < 1){
    if (is.na(chrmax)) chrmax <- 0
    return(Rle(0, chrmax))
  }
  ## sorting:
  ord <- order(x[,1], x[,2])
  ## get chromosome length or use maximum interval end
  if (is.na(chrmax)) chrmax <- max(x[,2])
  cp <- .Call(girafeCoverage,
              as.integer(x[ord,1]), as.integer(x[ord,2]),
              as.integer(reads(x)), as.integer(chrmax))
  ## construct result and return
  res <- Rle(values=cp[,1], lengths=cp[,2])
  return(res)
}# coverageOneStrand


setMethod("coverage", signature("AlignedGenomeIntervals"),
         function(x, shift=NA, width=NULL, weight=NA,
                  byStrand=FALSE, mem.friendly=FALSE, ...){
           # shift, width, weight have to be there because
           #  the generic definition from IRanges mentions them, however
           #  they are not used here; only byStrand is relevant and
           #  says whether one wants strand-wise coverage or combined
           #  over + and - strand (default)

           ## changed on 23 march 2010 to use coverage method from IRanges
           ##  (previous implemented commented out with #o#

           ## 1. over which sets (=spaces) to compute coverage:
           spacex <- chromosome(x)
           allChr <- unique(spacex)

           ## 2. get chromosome lengths
           chrlens <- getChromLengths(x)
           if (byStrand){
             # adapt chromosome lengths and spacex vector
             oldchrlens <- chrlens
             chrlens <- rep(chrlens, each=2)
             names(chrlens) <- paste(rep(names(oldchrlens), each=2),
                                     c("+","-"),sep="")
             spacex <- paste(spacex, strand(x), sep="")
             allChr <- unique(spacex)
           }#byStrand

           if (!all(allChr %in% names(chrlens)))
             stop("No lengths for some chromosome names!")
     
           # new part, using coverage for RangedData from IRanges
           #  (thanks to Patrick for the suggestion)
           if (!mem.friendly){ # default
             coords <- RangedData(IRanges(start=x[,1], end=x[,2]),
                                  space=spacex, reads=reads(x))
             byChr <- coverage(coords, weight="reads",
                               width=as.list(chrlens[allChr]))
           } else { # mem.friendly==TRUE
             
             ## which function to use for each iteration:
             if ("package:parallel" %in% search())
               lFun <- mclapply
             else 
               lFun <- lapply
             byChr <-
               lFun(structure(as.list(allChr), names=allChr),
                    function(z){
                      on.z <- which(spacex == z)
                      coords <- RangedData(IRanges(start=x[on.z,1],
                                                   end=x[on.z,2]),
                                           space=z, reads=reads(x)[on.z])
                      coverage(coords, width=as.list(chrlens[z]),
                               weight="reads")[[1]]
                    })
           }# mem.friendly==TRUE

           #o#  old implementation, slower and not up to date:
           #byChr <- lFun(as.list(allChr),  function(z){
           #  on.z <- chromosome(x) == z
           #  if (!byStrand)
           #    coverageOneStrand(x[on.z], chrmax=chrlens[z])
           #  else 
           #    list("minus"=coverageOneStrand(x[on.z & on.minus],
           #           chrmax=chrlens[z]),
           #         "plus"=coverageOneStrand(x[on.z & !on.minus],
           #           chrmax=chrlens[z]))
           #})
           #names(byChr) <- allChr
           ## convert to SimpleRleList
           #byChr <- RleList(byChr)
           
           return(byChr)
} ) #setMethod("coverage", signature("AlignedGenomeIntervals")
