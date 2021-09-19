#' Main class
#'
#' @param genome character defining the genome.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)

setClass("TREGELDataSet",
  slots=list(
  genome="character",
  query="GRanges",
  queryFolder="character"),

  prototype=list(
  genome = "Human",
  queryFolder=NA
  )
)
  # Make a function that can test to see if the data is consistent.
# This is not called if you have an initialize function defined!
validity=function(object)
{
  if(object@genome%in%c("Human","Mouse")) {
    return("Only Human or Mouse genomes are supported.")
  }
  return(TRUE)
}

#constructors
TREGELDataSetfromFolder

setGeneric("query", function(x) standardGeneric("query"))
setGeneric("query<-", function(x, value) standardGeneric("query<-"))
setMethod("query", "TREGELDataSet", function(x) x@query)
setMethod("query<-", "TREGELDataSet", function(x, value) {
  x@query <- value
  x})

setGeneric("genome", function(x) standardGeneric("genome"))
setGeneric("genome<-", function(x, value) standardGeneric("genome<-"))
setMethod("genome", "TREGELDataSet", function(x) x@genome)
setMethod("genome<-", "TREGELDataSet", function(x, value) {
  x@genome <- value
  x})

setGeneric("queryFolder", function(x) standardGeneric("queryFolder"))
setGeneric("queryFolder<-", function(x, value) standardGeneric("queryFolder<-"))
setMethod("queryFolder", "TREGELDataSet", function(x) x@queryFolder)
setMethod("queryFolder<-", "TREGELDataSet", function(x, value) {
  x@queryFolder <- value
  x})

