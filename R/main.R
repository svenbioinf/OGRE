#' BuildOGREDataSetFromDir
#'
#' Builds a `OGREDataset` from user specified directories containing datasets for which an overlap between query and subject is to be calculated.
#' A `OGREDataset` is a `GenomicRangesList` which stores datasets in a list like structure and possible metadata information.
#' @param queryFolder A \code{character} path pointing to the directory where your query dataset is located.
#' @param subjectFolder A \code{character} path pointing to the directory where your subject dataset(s) are located.
#' @import GenomicRanges methods ggplot2 Gviz S4Vectors
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table
#' @return A OGREDataSet.
#' @examples
#' myQueryFolder <- file.path(system.file('extdata', package = 'OGRE'),"query")
#' mySubjectFolder <- file.path(system.file('extdata', package = 'OGRE'),"subject")
#' myOGRE <- OGREDataSetFromDir(queryFolder=myQueryFolder,subjectFolder=mySubjectFolder)
#' @export
OGREDataSetFromDir <- function(queryFolder,subjectFolder){
  assertthat::assert_that(!is.null(queryFolder),msg="Please specify query folder!")
  assertthat::assert_that(!is.null(subjectFolder),msg="Please specify subject folder!")
  assertthat::assert_that(is(queryFolder,"character"),msg="queryFolder must be of type character!")
  assertthat::assert_that(is(subjectFolder,"character"),msg="subjectFolder must be of type character!")
  message("Initializing OGREDataSet... ")
  OGREDataSet <- GRangesList()
  metadata(OGREDataSet)$queryFolder <- queryFolder
  metadata(OGREDataSet)$subjectFolder <- subjectFolder
  dir.create(file.path(dirname(metadata(OGREDataSet)$queryFolder),"output"),
             recursive=TRUE) #move up one level
  metadata(OGREDataSet)$outputFolder <- file.path(dirname(metadata(OGREDataSet)$queryFolder), "output") 
  dir.create(file.path(dirname(metadata(OGREDataSet)$outputFolder),"gvizPlots"),
             recursive=TRUE) #move up one level
  metadata(OGREDataSet)$gvizPlotsFolder <- file.path(dirname(metadata(OGREDataSet)$outputFolder), "gvizPlots") 
  metadata(OGREDataSet)$subjectNames <- NULL
  metadata(OGREDataSet)$detailDT <- NULL
  metadata(OGREDataSet)$sumDT <- NULL
  metadata(OGREDataSet)$barplot_summary <- NULL
  metadata(OGREDataSet)$barplot_summary_dt <- NULL
  metadata(OGREDataSet)$itracks <- list()
  metadata(OGREDataSet)$aH <- NULL #annotation hub
  return(OGREDataSet)
}

#' BuildOGREDataSet
#'
#' Builds a `OGREDataset` as a `GenomicRangesList` for storing and analysing
#' datasets which can be added by `addDataSetFromHub()` or `addGRanges()`. 
#' Use `BuildOGREDataSetFromDir` for adding dataSets stored as files.
#' @import GenomicRanges methods ggplot2 Gviz S4Vectors
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table
#' @return A OGREDataSet.
#' @examples
#' myOGRE <- OGREDataSet()
#' @export
OGREDataSet <- function(){
  message("Initializing OGREDataSet... ")
  OGREDataSet <- GRangesList()
  metadata(OGREDataSet)$queryFolder <- NULL
  metadata(OGREDataSet)$subjectFolder <- NULL
  dir.create(file.path(getwd(), "OGRE/output"), #create default dir
             recursive = TRUE)
  metadata(OGREDataSet)$outputFolder <- file.path(getwd(), "OGRE/output")
  dir.create(file.path(dirname(metadata(OGREDataSet)$outputFolder),"gvizPlots"),
             recursive = TRUE) 
  metadata(OGREDataSet)$gvizPlotsFolder <- file.path(dirname(metadata(OGREDataSet)$outputFolder), "gvizPlots") 
  metadata(OGREDataSet)$subjectNames <- NULL
  metadata(OGREDataSet)$detailDT <- NULL
  metadata(OGREDataSet)$sumDT <- NULL
  metadata(OGREDataSet)$barplot_summary <- NULL
  metadata(OGREDataSet)$barplot_summary_dt <- NULL
  metadata(OGREDataSet)$itracks <- list()
  metadata(OGREDataSet)$aH <- NULL #annotation hub
  return(OGREDataSet)
}


