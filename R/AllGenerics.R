### generic functions used in package girafe

setGeneric("id<-", function(object, value) standardGeneric("id<-"))

setGeneric("organism<-", function(x, value) standardGeneric("organism<-"))

setGeneric("export", function(object, con, format, ...) standardGeneric("export"))
# deliberately similar to the generic definition in rtracklayer

setGeneric("reads", function(x) standardGeneric("reads"))

setGeneric("reads<-", function(x, value) standardGeneric("reads<-"))

setGeneric("matches", function(x) standardGeneric("matches"))

setGeneric("matches<-", function(x, value) standardGeneric("matches<-"))

setGeneric("extend", function(x, ...) standardGeneric("extend"))

## new slot 'chrlengths' vector of chromosome lengths
setGeneric("chrlengths", function(x) standardGeneric("chrlengths"))

setGeneric("chrlengths<-", function(x, value) standardGeneric("chrlengths<-"))

## standard generics
setGeneric("plot")
setGeneric("summary")
