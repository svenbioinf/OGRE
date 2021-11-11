#' BuildTREGELDataSetFromDir
#'
#' Builds a `TREGELDataset` from user specified directories containing datasets for which an overlap between query and subject is to be calculated.
#' A `TREGELDataset` is a `GenomicRangesList` with stores datasets in a list like structure and possible metadata information.
#' @param queryFolder A \code{character} path pointing to the directory where your query dataset is located.
#' @param subjectFolder A \code{character} path pointing to the directory where your subject dataset(s) are located.
#' @import GenomicRanges S4Vectors methods data.table ggplot2 Gviz IRanges
#' @importFrom assertthat assert_that
#' @export
TREGELDataSetFromDir <- function(queryFolder,subjectFolder){
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
  metadata(TREGELDataSet)$subjectNames <- NULL
  metadata(TREGELDataSet)$detailDT <- NULL
  metadata(TREGELDataSet)$sumDT <- NULL
  metadata(TREGELDataSet)$barplot_summary <- NULL
  metadata(TREGELDataSet)$barplot_summary_dt <- NULL
  metadata(TREGELDataSet)$itracks <- list()


  return(TREGELDataSet)
}



