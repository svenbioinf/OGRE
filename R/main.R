#' BuildTREGELDataSet
#'
#' @import GenomicRanges S4Vectors methods data.table ggplot2 Gviz IRanges
#' @importFrom assertthat assert_that
#' @export
TREGELDataSet <- function(queryFolder,subjectFolder){
  assertthat::assert_that(!is.null(queryFolder),msg="Please specify query folder!")
  assertthat::assert_that(!is.null(subjectFolder),msg="Please specify subject folder!")
  assertthat::assert_that(is(queryFolder,"character"),msg="queryFolder must be of type character!")
  assertthat::assert_that(is(subjectFolder,"character"),msg="subjectFolder must be of type character!")


  message("Initializing TREGELDataSet... ")
  TREGELDataSet <- GRangesList()
  metadata(TREGELDataSet)$queryFolder <- queryFolder
  metadata(TREGELDataSet)$subjectFolder <- subjectFolder
                                          dir.create(file.path(dirname(metadata(TREGELDataSet)$queryFolder), "output")) #move up one level
  metadata(TREGELDataSet)$outputFolder <- file.path(dirname(metadata(TREGELDataSet)$queryFolder), "output") #move up one level
                                          dir.create(file.path(dirname(metadata(TREGELDataSet)$outputFolder), "gvizPlots")) #move up one level
  metadata(TREGELDataSet)$gvizPlotsFolder <- file.path(dirname(metadata(TREGELDataSet)$outputFolder), "gvizPlots") #move up one level
  metadata(TREGELDataSet)$detailDT <- NA
  metadata(TREGELDataSet)$sumDT <- NA
  metadata(TREGELDataSet)$barplot_summary <- NA
  metadata(TREGELDataSet)$barplot_summary_dt <- NA
  metadata(TREGELDataSet)$itracks <- list()


  return(TREGELDataSet)
}



