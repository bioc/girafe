
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
                  byStrand=FALSE, ...){
           # shift, width, weight have to be there because
           #  the generic definition from IRanges mentions them, however
           #  they are not used here; only byStrand is relevant and
           #  says whether one wants strand-wise coverage or combined
           #  over + and - strand (default)

           ## changed on 23 march 2010 to use coverage method from IRanges
           ##  (previous implemented commented out with #o#
           
           if (byStrand){
             #o# on.minus <- strand(x)=="-"
             altStrandCoding <-
               factor(as.character(strand(x)), levels=c("-", "+"),
                      labels=c("minus", "plus"))
           }
 
           ## which chrosomes are in the data:
           allChr <- gsub("MT","M",unique(chromosome(x)))
           allChr <- paste("chr",gsub("^chr","",allChr),sep="")
           ## if organism defined, chromomse length specifies
           ##  the length of the coverage vectors:
           if (length(x@organism)>0){
             chrlens <- getChromLengths(x)
           } else {
             chrlens <- rep.int(as.integer(NA), length(allChr))
             names(chrlens) <- allChr
           }
           ## which function to use for each iteration:
           if ("package:multicore" %in% search()){
             lFun <- mclapply
           } else {
             lFun <- lapply
           }

           # new part, using coverage for RangesList from IRanges
           #  (thanks to Patrick for the suggestion)
           coords <- IRanges(start=x[,1], end=x[,2])
           if (!byStrand){
             splitCoords <- split(coords, factor(chromosome(x), levels=allChr))
             byChr <- coverage(splitCoords,
                               width=as.list(chrlens)[match(allChr, names(chrlens))])
           } else {
             byChr <-
               lFun(structure(as.list(allChr), names=allChr),
                    function(z){
                      on.z <- chromosome(x) == z
                      splitCoords <- split(coords[on.z], altStrandCoding[on.z])
                      lens <- list("minus"=chrlens[z], "plus"=chrlens[z])
                      coverage(splitCoords, width=lens)
                    })
           }

           #o#  old implementation, slower
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
