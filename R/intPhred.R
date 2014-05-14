intPhred <- function(x, method="Sanger", returnType="list"){
  stopifnot(inherits(x, "ShortReadQ"))
  method <- match.arg(method, c("Sanger", "Solexa", "previousSolexa"))
  returnType <- match.arg(returnType, c("list","matrix", "Rle"))
  conv.fun <- switch(method,
                     "Sanger"=function(z){ utf8ToInt(z)-33},
                     "Solexa"=function(z){ utf8ToInt(z)-64},
                     "previousSolexa"=function(z){
                       Q <- utf8ToInt(z)-64
                       round(10*log10(1+10^(Q/10))) })
  ## which function to use for each iteration:
  if ("package:parallel" %in% search()){
    lFun <- mclapply
  } else {
    lFun <- lapply
  }
  xqualchar <- as.character(quality(quality(x)))
  if (returnType=="list"){
    xl <- lFun(as.list(xqualchar), conv.fun)
    names(xl) <- as.character(id(x))
  }
  if (returnType=="Rle"){
    xl <- lFun(as.list(xqualchar), function(z) Rle(conv.fun(z)))
    names(xl) <- as.character(id(x))
  }
  if (returnType=="matrix"){
    maxWidth <- max(width(x))
    xl <- lFun(as.list(xqualchar), function(z) conv.fun(z)[1:maxWidth])
    xl <- do.call("rbind", xl)
    rownames(xl) <- as.character(id(x))
  }
  return(xl)
}#intPhred


## function for iterative median computation over read positions
medianByPosition <- function(x, method="Sanger", batchSize=100000L){
  stopifnot(inherits(x, "ShortReadQ"))
  if (length(x) <= batchSize){
    ends <- c(0L, length(x))
  } else {
    ends <- seq.int(0L, length(x), by=batchSize)
    ## is there some rest?
    if (length(x) %% batchSize != 0)
      ends <- c(ends, length(x))
  }
  ## initialise result:
  qualsPerPos <- lapply(1:max(width(x)),
                        function(z) new("Rle", values=numeric(0)))
  for (i in 1:(length(ends)-1)){
    theseQuals <- intPhred(x[(ends[i]+1L):ends[i+1L]],
                           method=method, returnType="matrix")
    for (j in 1:ncol(theseQuals))
      ## in each iteration: attach new values as sorted Rle object
      ##  to save memory
      qualsPerPos[[j]] <- sort(c(qualsPerPos[[j]],
                                 Rle(theseQuals[,j])))
  }
  ## compute medians:
  medianPerPos <- sapply(qualsPerPos, median, na.rm=TRUE)
  return(medianPerPos)
}##medianByPosition
