intPhred <- function(x, method="Sanger", returnType="list"){
  stopifnot(inherits(x, "ShortReadQ"))
  method <- match.arg(method, c("Sanger", "Solexa", "previousSolexa"))
  returnType <- match.arg(returnType, c("list","matrix"))
  conv.fun <- switch(method,
                     "Sanger"=function(z){ utf8ToInt(z)-33},
                     "Solexa"=function(z){ utf8ToInt(z)-64},
                     "previousSolexa"=function(z){
                       Q <- utf8ToInt(z)-64
                       log10(1+10^(Q/10)) })
                       #round(-10*log10((10^(-Q/10))/(1+10^(-Q/10))))})
  ## which function to use for each iteration:
  if ("package:multicore" %in% search()){
    lFun <- mclapply
  } else {
    lFun <- lapply
  }
  xqualchar <- as.character(quality(quality(x)))
  if (returnType=="list"){
    xl <- lFun(as.list(xqualchar), conv.fun)
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
