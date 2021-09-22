#' BuildTREGELDataSet
#'
#' @import GenomicRanges S4Vectors methods
#' @export
TREGELDataSet <- function(queryFolder){
  TREGELDataSet <- GRangesList()
  metadata(TREGELDataSet)$queryFolder=queryFolder


  return(TREGELDataSet)
}



