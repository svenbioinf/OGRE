#' OGRE package to calculate, analyze and visiualize overlap between annotated 
#' genomic region datasets
#' 
#'OGRE calculates overlap between user defined annotated genomic region 
#'datasets. Any regions can be supplied such as public annotations (genes), 
#'genetic variation (SNPs, mutations), regulatory elements (TFBS, promoters, 
#'CpG islands) and basically all  types of NGS output from sequencing 
#'experiments. After overlap calculation, key numbers help analyse the extend 
#'of overlaps which can also be visualized at a genomic level.
#'
#' The main functions are:
#'
#' \code{\link{OGREDataSetFromDir}} - build an OGRE dataset from a user 
#' defined directory with GRanges annotation files.
#' \itemize{
#' \item \code{\link{loadAnnotations}} - Load dataset files containing genomic 
#' regions annotation information from hard drive
#' }
#' 
#' \code{\link{OGREDataSet}} - build an empty OGRE dataset to flexibly add
#' datasets from other sources like AnnotationHub or custom GRanges objects.
#'  \itemize{
#'  \item \code{\link{addDataSetFromHub}} - adds datasets from AnnotationHub
#'  \item \code{\link{addGRanges}} - adds user defined GenomicRanges datasets
#'  }
#' 
#' \code{\link{fOverlaps}} - Finds all overlaps between query and subject
#' datasets\cr
#' \code{\link{sumPlot}} - calculates key numbers, tables and plots\cr
#' \code{\link{gvizPlot}} - generates a genomic plot around query elements with 
#' overlapping subject hits.
#' 
#' For additional information, see the package vignette, by typing
#' \code{vignette("OGRE")}. Software-related questions or issues can be posted 
#' to the Bioconductor Support Site:
#' 
#' \url{https://support.bioconductor.org}
#'
#' or on github:
#'
#' \url{https://https://github.com/svenbioinf/OGRE}\cr
#' @author Sven Berres, Jörg Gromoll, Marius Wöste, Sarah Sandmann, Sandra
#' Laurentino
#' @docType package
#' @name OGRE-package
#' @aliases OGRE-package
#' @keywords package
NULL

#' BuildOGREDataSetFromDir
#'
#' Builds a `OGREDataset` from user specified directories containing datasets 
#' for which an overlap between query and subject is to be calculated.
#' A `OGREDataset` is a `GenomicRangesList` which stores datasets in a list like
#' structure and possible metadata information.
#' @param queryFolder A \code{character} path pointing to the directory where 
#' your query dataset is located.
#' @param subjectFolder A \code{character} path pointing to the directory where 
#' your subject dataset(s) are located.
#' @import GenomicRanges methods ggplot2 S4Vectors
#' @importFrom assertthat assert_that
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
  suppressWarnings(
  dir.create(file.path(dirname(metadata(OGREDataSet)$queryFolder),"output"),
             recursive=TRUE)) #dirname returns parentfolder
  metadata(OGREDataSet)$outputFolder <- file.path(dirname(metadata(OGREDataSet)$queryFolder), "output") 
  suppressWarnings(
  dir.create(file.path(dirname(metadata(OGREDataSet)$outputFolder),"gvizPlots"),
             recursive=TRUE))
  metadata(OGREDataSet)$gvizPlotsFolder <- file.path(dirname(metadata(OGREDataSet)$outputFolder), "gvizPlots") 
  metadata(OGREDataSet)$subjectNames <- NULL
  metadata(OGREDataSet)$detailDT <- NULL
  metadata(OGREDataSet)$sumDT <- NULL
  metadata(OGREDataSet)$barplot_summary <- NULL
  metadata(OGREDataSet)$barplot_summary_log <- NULL
  metadata(OGREDataSet)$barplot_summary_dt <- NULL
  metadata(OGREDataSet)$hist <- NULL
  metadata(OGREDataSet)$hist_dt <- NULL
  metadata(OGREDataSet)$quickDT <- NULL
  metadata(OGREDataSet)$summaryDT <- list()
  metadata(OGREDataSet)$itracks <- list()
  metadata(OGREDataSet)$aH <- NULL #annotation hub
  metadata(OGREDataSet)$covPlot <- NULL
  return(OGREDataSet)
}

#' BuildOGREDataSet
#'
#' Builds a `OGREDataset` as a `GenomicRangesList` for storing and analysing
#' datasets which can be added by `addDataSetFromHub()` or `addGRanges()`. 
#' Use `BuildOGREDataSetFromDir` for adding dataSets stored as files.
#' @import GenomicRanges methods ggplot2 S4Vectors
#' @importFrom assertthat assert_that
#' @return A OGREDataSet.
#' @examples
#' myOGRE <- OGREDataSet()
#' @export
OGREDataSet <- function(){
  message("Initializing OGREDataSet... ")
  OGREDataSet <- GRangesList()
  metadata(OGREDataSet)$queryFolder <- NULL
  metadata(OGREDataSet)$subjectFolder <- NULL
  suppressWarnings(
  dir.create(file.path(getwd(), "OGRE/output"), #create default dir
             recursive = TRUE))
  metadata(OGREDataSet)$outputFolder <- file.path(getwd(), "OGRE/output")
  suppressWarnings(
  dir.create(file.path(dirname(metadata(OGREDataSet)$outputFolder),"gvizPlots"),
             recursive = TRUE))
  metadata(OGREDataSet)$gvizPlotsFolder <- file.path(dirname(metadata(OGREDataSet)$outputFolder), "gvizPlots") 
  metadata(OGREDataSet)$subjectNames <- NULL
  metadata(OGREDataSet)$detailDT <- NULL
  metadata(OGREDataSet)$sumDT <- NULL
  metadata(OGREDataSet)$quickDT <- NULL
  metadata(OGREDataSet)$barplot_summary <- NULL
  metadata(OGREDataSet)$barplot_summary_log <- NULL
  metadata(OGREDataSet)$barplot_summary_dt <- NULL
  metadata(OGREDataSet)$hist <- NULL
  metadata(OGREDataSet)$hist_dt <- NULL
  metadata(OGREDataSet)$summaryDT <- list()
  metadata(OGREDataSet)$itracks <- list()
  metadata(OGREDataSet)$aH <- NULL #annotation hub
  metadata(OGREDataSet)$covPlot <- NULL
  return(OGREDataSet)
}


