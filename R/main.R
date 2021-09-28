#' BuildTREGELDataSet
#'
#' @import GenomicRanges S4Vectors methods data.table
#' @importFrom assertthat assert_that
#' @export
TREGELDataSet <- function(queryFolder,subjectFolder){
  TREGELDataSet <- GRangesList()
  metadata(TREGELDataSet)$queryFolder <- queryFolder
  metadata(TREGELDataSet)$subjectFolder <- subjectFolder
  metadata(TREGELDataSet)$detailDT <- NA
  metadata(TREGELDataSet)$sumDT <- NA

  metadata(TREGELDataSet)$itracks <- list()


  return(TREGELDataSet)
}



