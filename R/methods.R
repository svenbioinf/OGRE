#' Load annotation datasets
#'
#' Load dataset files containing genomic regions annotation information from 
#' hard drive. `loadAnnotations` calls `readQuery` and `readSubject`
#' to read in genomic regions as `GenomicRanges` objects stored as .RDS / .rds 
#' files. Each region needs chromosome, start, end and strand information.
#' A unique ID and a name column must be present in the `GenomicRanges` object 
#' metadata. OGRE searches for the query file in your query folder and any 
#' number of subject files in your subjects folder. Alternatively, .gff (v2&v3) files
#' in the query or subject folder with attribute columns containing "ID" and "name" 
#' information are read in by OGRE.
#' @param OGREDataSet A OGREDataSet.
#' @return A OGREDataSet.
#' @examples
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' @export
loadAnnotations <- function(OGREDataSet){
  OGREDataSet <- readQuery(OGREDataSet)
  OGREDataSet <- readSubject(OGREDataSet)

 return(OGREDataSet)
}




#' Read query dataset
#'
#' [readQuery()] scanns `queryFolder` for a `GRanges` object stored as .RDS/.rds
#' or .gff .GFF file and attaches it to the OGREDataSet.
#' @importFrom rtracklayer import.gff
#' @param OGREDataSet A OGREDataSet.
#' @return A OGREDataSet.
#' @keywords internal
readQuery=function(OGREDataSet){
  assertthat::assert_that(length(list.files(metadata(OGREDataSet)$queryFolder))>0,
                          msg="Query folder is empty!")
  message("Reading query dataset... ")
  if(is.null(metadata(OGREDataSet)$queryFolder)){ #if no folder supplied, check default folder

  }else{ #read queryFolder
    queryPath <- list.files(metadata(OGREDataSet)$queryFolder,full.names = TRUE)
    queryName <- tools::file_path_sans_ext(list.files(metadata(OGREDataSet)$queryFolder))
    if(grepl("rds",queryPath,ignore.case=TRUE)){
    OGREDataSet[[queryName]]<- readRDS(queryPath)
    }
    else if(grepl("gff",queryPath,ignore.case=TRUE)){ #read in .gff
      OGREDataSet[[queryName]] <- rtracklayer::import.gff(queryPath)%>%
      GenomeInfoDb::keepStandardChromosomes("Homo_sapiens",pruning.mode="coarse")
      GenomeInfoDb::seqlevelsStyle(OGREDataSet[[queryName]]) <- "Ensembl"
    }
  }
  assertthat::assert_that(c("ID")%in%names(mcols(OGREDataSet[[queryName]])),
                          msg="Query must contain ID column.")
  assertthat::assert_that(!any(duplicated(OGREDataSet[[queryName]]$ID)),
                          msg="ID column must be unique.")
  assertthat::assert_that(length(OGREDataSet[[queryName]])!=0,
                          msg=paste0("Dataset has no ranges: ",queryName))

  return(OGREDataSet)
}
#' Read subject datasets
#'
#' [readSubject()] scanns `SubjectFolder` for `GRanges` objects stored as .RDS/.rds
#' or .gff .GFF files and attaches them to the OGREDataSet.
#' @param OGREDataSet A OGREDataSet.
#' @return A OGREDataSet.
#' @keywords internal
readSubject=function(OGREDataSet){
  assertthat::assert_that(length(list.files(metadata(OGREDataSet)$subjectFolder))>0,
                          msg="Subject folder is empty!")
  message("Reading subject datasets... ")
  #read SubjectFolder
    subjectPath <- list.files(metadata(OGREDataSet)$subjectFolder,full.names = TRUE)
    subjectName <- tools::file_path_sans_ext(list.files(metadata(OGREDataSet)$subjectFolder))
    OGREDataSet <- c(OGREDataSet,mapply(function(x,y){
      if(grepl("rds",y,ignore.case=TRUE)){
        tmp=readRDS(y)
      }
      else if(grepl("gff",y,ignore.case=TRUE)){ #read in .gff
        tmp <- rtracklayer::import.gff(y)%>%
          GenomeInfoDb::keepStandardChromosomes("Homo_sapiens",pruning.mode="coarse")
        GenomeInfoDb::seqlevelsStyle(tmp) <- "Ensembl"
      }
      assertthat::assert_that(c("ID")%in%names(mcols(tmp)),msg="Subject must contain ID column.")
      assertthat::assert_that(!any(duplicated(tmp$ID)),msg="ID column must be unique.")
      assertthat::assert_that(length(tmp)!=0,msg=paste0("Dataset has no ranges: ",x))
      tmp
    },x=subjectName,y=subjectPath))
  return(OGREDataSet)
}

#' Read dataset(s) from folder
#'
#' [readDataSetFromFolder()] scanns queryFolder and subjectFolder for either
#' .RDS/.rds or .CSV/.csv files and adds them to a OGREDataSet. Each region 
#' needs chromosome, start, end and strand information. (tabular file columns
#' must be named accordingly). A unique ID and a name column must be present in 
#' the `GenomicRanges` object's metatdata and tabular file.
#' @param OGREDataSet A OGREDataSet.
#' @param type \code{character} and one of query/subject.
#' @return A OGREDataSet.
#' @export
#' @examples 
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- readDataSetFromFolder(myOGRE,type="query")
#' myOGRE <- readDataSetFromFolder(myOGRE,type="subject")
readDataSetFromFolder=function(OGREDataSet,type){
  message("Reading dataset(s)... ")
  if(type=="query"){
    files <- list.files(metadata(OGREDataSet)$queryFolder)
    assertthat::assert_that(length(files)>0,msg="Folder is empty!")
    assertthat::assert_that(grep("rds|RDS|csv|CSV",files)>0,msg="Wrong file type!")
    queryPath <- list.files(metadata(OGREDataSet)$queryFolder,full.names = TRUE)
    queryName <- tools::file_path_sans_ext(list.files(metadata(OGREDataSet)$queryFolder))
    if(grep("rds|RDS",files)>0){
    OGREDataSet[[queryName]]<- readRDS(queryPath)}
    else{
    OGREDataSet[[queryName]]<- makeGRangesFromDataFrame(fread(queryPath),
                               keep.extra.columns = TRUE)
    }
    assertthat::assert_that(c("ID")%in%names(mcols(OGREDataSet[[queryName]])),
                            msg="Query must contain ID column.")
    assertthat::assert_that(!any(duplicated(OGREDataSet[[queryName]]$ID)),
                            msg="ID column must be unique.")
    assertthat::assert_that(length(OGREDataSet[[queryName]])!=0,
                            msg=paste0("Dataset has no ranges: ",queryName))
  }
  else if(type=="subject"){
    files <- list.files(metadata(OGREDataSet)$subjectFolder)
    assertthat::assert_that(length(files)>0,msg="Folder is empty!")
    assertthat::assert_that(length(grep("rds|RDS|csv|CSV",files))==length(files),
                            msg="Wrong file type!")
    subjectPath <- list.files(metadata(OGREDataSet)$subjectFolder,full.names = TRUE)
    subjectName <- tools::file_path_sans_ext(list.files(metadata(OGREDataSet)$subjectFolder))
    if(length(grep("rds|RDS",files))==length(files)){
    OGREDataSet <- c(OGREDataSet,mapply(function(x,y){
      tmp=readRDS(y)
      assertthat::assert_that(c("ID")%in%names(mcols(tmp)),msg="Subject must contain ID column.")
      assertthat::assert_that(!any(duplicated(tmp$ID)),msg="ID column must be unique.")
      assertthat::assert_that(length(tmp)!=0,msg=paste0("Dataset has no ranges: ",x))
      tmp
    },x=subjectName,y=subjectPath))
    }else{
    OGREDataSet <- c(OGREDataSet,mapply(function(x,y){
      tmp <-makeGRangesFromDataFrame(fread(y),keep.extra.columns = TRUE)
      assertthat::assert_that(c("ID")%in%names(mcols(tmp)),msg="Subject must contain ID column.")
      assertthat::assert_that(!any(duplicated(tmp$ID)),msg="ID column must be unique.")
      assertthat::assert_that(length(tmp)!=0,msg=paste0("Dataset has no ranges: ",x))
      tmp
    },x=subjectName,y=subjectPath))
    }
  }
  return(OGREDataSet)
}



#' Find overlaps
#'
#' Finds all overlaps between query and subject(s) and stores each hit (overlap)
#' in data table `detailDT`. Data table `sumDT` shows all
#' overlaps of a certain subject type for all query elements. By default also 
#' partially overlaps are reported. Overlap calculation is done using 
#' \code{GenomicRanges::findOverlaps()} implementation.
#'
#' @param OGREDataSet A OGREDataSet.
#' @param selfHits \code{logical} if FALSE(default) ignores self hits of 
#' identical regions (with identical IDs) within datasets. 
#' @param ignoreStrand \code{logical} If TRUE (default) two regions with 
#' overlapping locations on different strands are considered an overlap hit.
#' @param ... Additional parameters, see \code{GenomicRanges::findOverlaps()}
#' @return OGREDataSet.
#' @examples
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' myOGRE <- fOverlaps(myOGRE)
#' @export
fOverlaps <- function(OGREDataSet,selfHits=FALSE,ignoreStrand=TRUE,...){
  queryID <- subjID <- subjType <- . <-  NULL
  detailDT <- data.table() #data table to store all overlaps for query vs all subjects
  metadata(OGREDataSet)$subjectNames <- names(
    OGREDataSet[seq(2,length(OGREDataSet))])
  for(subj in metadata(OGREDataSet)$subjectNames){
    ol <- findOverlaps(OGREDataSet[[1]],OGREDataSet[[subj]],ignore.strand=ignoreStrand,...)
    if(length(ol)==0){ #if no overlap hits skip to next subj
      message("No overlap found for: ",subj)
      next
    }
    else{
      overlapWidth <- width(pintersect(OGREDataSet[[1]][queryHits(ol)], 
              OGREDataSet[[subj]][subjectHits(ol)],ignore.strand=ignoreStrand))
      #overlap ratio with reference to query.(overLen=2,queryLen=6, ratio=0.33)
      #https://support.bioconductor.org/p/72656/
      overlapRatio <- overlapWidth / width(OGREDataSet[[1]][queryHits(ol)])
      ol <- data.table(q=queryHits(ol),s=mcols(OGREDataSet[[subj]])[subjectHits(ol),1])
      detailDT <- rbind(detailDT,data.table(
       queryID=mcols(OGREDataSet[[1]])[ol$q,"ID"],
       queryType=names(OGREDataSet)[1],
       subjID=ol$s,
       subjType=subj,
       queryChr=as.character(seqnames(OGREDataSet[[names(OGREDataSet)[1]]]))[ol$q],
       queryStart=start(OGREDataSet[[names(OGREDataSet)[1]]])[ol$q],
       queryEnd=end(OGREDataSet[[names(OGREDataSet)[1]]])[ol$q],
       queryStrand=as.character(strand(OGREDataSet[[names(OGREDataSet)[1]]]))[ol$q],
       subjChr=as.character(seqnames(OGREDataSet[[subj]]))
       [match(ol$s,mcols(OGREDataSet[[subj]])[,"ID"])],
       subjStart=start(OGREDataSet[[subj]])
       [match(ol$s,mcols(OGREDataSet[[subj]])[,"ID"])],
       subjEnd=end(OGREDataSet[[subj]])
       [match(ol$s,mcols(OGREDataSet[[subj]])[,"ID"])],
       subjStrand=as.character(strand(OGREDataSet[[subj]]))
       [match(ol$s,mcols(OGREDataSet[[subj]])[,"ID"])],
       overlapWidth=overlapWidth,
       overlapRatio=overlapRatio
      ))
      if(isFALSE(selfHits)){ #remove self hits
        detailDT <- detailDT[!(queryID==subjID)]
      }
    }
  }
  metadata(OGREDataSet)$sumDT <- detailDT[,
  list(subjID=paste(subjID,collapse = " ")),by=.(queryID,subjType)]
  assert_that(nrow(metadata(OGREDataSet)$sumDT)!=0,
              msg="No overlaps found for any subjects!")
  metadata(OGREDataSet)$detailDT <- detailDT
  metadata(OGREDataSet)$quickDT <- detailDT[,c("queryID","subjType")]
  metadata(OGREDataSet)$quickDT <- data.table::dcast(data=metadata(OGREDataSet)$quickDT, 
                          queryID ~ subjType,value.var = "subjType",fun=length)
  #add and set query elements without hits to 0
  withoutHits <- setdiff(mcols(OGREDataSet[[1]])$ID,metadata(OGREDataSet)$quickDT$queryID)
  metadata(OGREDataSet)$quickDT <- rbind(metadata(OGREDataSet)$quickDT,
                                    data.table(queryID=withoutHits),fill=TRUE)
  metadata(OGREDataSet)$quickDT[is.na(metadata(OGREDataSet)$quickDT)] <- 0
  #set order to match other tables
  metadata(OGREDataSet)$quickDT <- metadata(OGREDataSet)$quickDT[
    match(mcols(OGREDataSet[[1]])$ID,metadata(OGREDataSet)$quickDT$queryID),]
  
  return(OGREDataSet)
}

#' Generate summary plot
#'
#' `sumPlot()` calculates key numbers i.e. 
#' (total number of overlaps, number of overlaps per subject...) to help with an
#' exploratory data evaluation and displays them in an informative barplot.
#' @importFrom AnnotationHub query
#' @param OGREDataSet A OGREDataSet.
#' @return OGREDataSet.
#' @examples
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' myOGRE <- fOverlaps(myOGRE)
#' myOGRE <- sumPlot(myOGRE)
#' @export
sumPlot <- function(OGREDataSet){
  assert_that(!is.null(metadata(OGREDataSet)$sumDT),
              msg="No data to generate sumPlot!")
  xtext <- lab <- NULL
  #summary barplot
  query_N <- length(unique(mcols(OGREDataSet[[1]])$ID))
  query_ol <- length(unique(metadata(OGREDataSet)$sumDT$queryID))
  query_No_ol <- query_N-query_ol
  query_full <- sum(table(metadata(OGREDataSet)$sumDT$queryID)>=length(OGREDataSet)-1)
  dtbarplot<-data.table(x=c("query_N","query_w_minOne_overlap",
        "genes_w_overlap_allSubjTypes"),query=c(query_N,query_ol,query_full))
  dtbarplot$x<-factor(x = dtbarplot$x,levels = c("query_N","query_w_minOne_overlap",
                                                 "genes_w_overlap_allSubjTypes"))
  dtbarplot$lab<-c(paste0("Total number of submitted queries: ",query_N," (100%)"),
                  paste0("Queries with >=1 subject overlaps: ",query_ol," (",
                         round(query_ol/query_N*100),"%)"),
                  paste0("Queries with overlaps in all subject types: ",
                         query_full," (",round(query_full/query_N*100),"%)"))
  dtbarplot$col<-"steelblue2"
  temp<-vapply(metadata(OGREDataSet)$subjectNames,function(x){
    length(unique(metadata(OGREDataSet)$detailDT[
      metadata(OGREDataSet)$detailDT$subjType==x,][["queryID"]]))},integer(1))
  temp<-data.table(x=names(temp),
   Subjects=temp,
   lab=paste0("Queries with ",names(temp),": ",temp," (",round(temp/query_N*100),"%)"),
   col="steelblue3")
  subjPerQuery<-round(vapply(metadata(OGREDataSet)$subjectNames,function(x){
  sum(metadata(OGREDataSet)$detailDT$subjType==x)},numeric(1))/query_N,digits=1)
  subjPerQuery<-data.table(x=names(subjPerQuery),Subjects=subjPerQuery,
  lab=paste0("Average ",names(subjPerQuery)," per query: ",subjPerQuery),col="steelblue4")
  dtbarplotDetailed<-rbind(dtbarplot,temp,subjPerQuery, use.names=FALSE)
  dtbarplotDetailed$sequence<-seq(dim(dtbarplotDetailed)[1],1)
  dtbarplotDetailed$xtext<-dtbarplotDetailed$sequence+0.55
  metadata(OGREDataSet)$barplot_summary <-
    ggplot(dtbarplotDetailed, aes(x = sequence, y =query,fill=col)) +
      geom_bar(stat = "identity", width = 0.3) +
      scale_fill_manual(guide="none",values = c("steelblue2" = "steelblue2",
                                               "steelblue3" = "steelblue3",
                                               "steelblue4" = "steelblue4"))+
      coord_flip() +labs(x = "", y = "N") +theme_bw() +  theme_classic() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      geom_text(aes(x = xtext, y = 0, label = lab),size=7,hjust = 0, vjust = 1)+
    guides(size = "none")+theme(axis.text = element_text(size = 14)) 
    ggsave(plot=metadata(OGREDataSet)$barplot_summary,
         path=metadata(OGREDataSet)$outputFolder,width = 30,
         filename="Barplot_Summary.png",height = 20,dpi = 400,units = "cm")
  dtbarplotDetailed$queryLog<-log2(dtbarplotDetailed$query+1)#add +1 to avoid log2(1)=0 and log2(0)=-inf
  metadata(OGREDataSet)$barplot_summary_dt <- dtbarplotDetailed
  metadata(OGREDataSet)$barplot_summary_log2 <-
    ggplot(dtbarplotDetailed, aes(x = sequence, y =queryLog,fill=col)) +
    geom_bar(stat = "identity", width = 0.3) +
    scale_fill_manual(guide="none",values = c("steelblue2" = "steelblue2",
                                              "steelblue3" = "steelblue3",
                                              "steelblue4" = "steelblue4"))+
    coord_flip() +labs(x = "", y = "Log2(N)") +theme_bw() +  theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    geom_text(aes(x = xtext, y = 0, label = lab),size=7,hjust = 0, vjust = 1)+
    guides(size = "none")+theme(axis.text = element_text(size = 14))
    ggsave(plot=metadata(OGREDataSet)$barplot_summary_log,
           path=metadata(OGREDataSet)$outputFolder,width = 30,
           filename="Barplot_Summary_log.png",height = 20,dpi = 400,units = "cm")


    return(OGREDataSet)
}

#' Generate Gviz plot
#'
#'`gvizPlot` generates a plot around one or many given query elements with all 
#'overlapping subject hits. In addition, each generated plot can be stored in 
#'the `gvizPlots` folder get or set by \code{gvizPlotsFolder}. A maximum of 25
#'elements can be plotted per track.
#' @importFrom grDevices dev.off pdf
#' @importFrom stats setNames
#' @rawNamespace import(Gviz, except = tags)
#' @param OGREDataSet A OGREDataSet.
#' @param query A character vector of one or many query elements ID's (i.e. Gene ID's).
#' @param gvizPlotsFolder A character pointing to the plot(s) output directory. 
#' If not supplied a folder is automatically generated and can be accessed by 
#' \code{metatdata(OGREDataSet)$gvizPlotsFolder}.
#' @param showPlot \code{logical} If \code{FALSE}(default) plots are only saved 
#' to `gvizPlotsFolder`. If \code{TRUE}
#' plots are additionally send to the plotting window.
#' @param trackRegionLabels A labeled character vector that defines the type of 
#' label that is displayed for query and subject elements during plotting. 
#' Vector values represent the type of label and vector labels define the type 
#' of subject element. In the following example 
#' \code{setNames(c("ID","name"),c("genes","CGI"))}
#' Value "ID" and label "genes" would annotate your genes with IDs taken from the
#' ID column of your dataset. Datasets not defined in this vector are 
#' plotted without track labels.
#' @param trackShapes A labeled character vector that defines the type of 
#' shape in which every dataset's elements are displayed. 
#' Vector values represent the type of shape and vector labels define the type 
#' of subject element. In the following example 
#' \code{setNames(c("fixedArrow","box"),c("genes","CGI"))}
#' Value "fixedArrow" and label "genes" would display your genes in fixedArrow 
#' and CGI as box shape. Possible values:
#' (box, arrow, fixedArrow, ellipse, and smallArrow) Default="fixedArrow"
#' @param extendPlot \code{int vector} Integer vector of length two that extends
#' the plot window to the left or right by adding the first value to query start
#' and the second value to query end coordinates(bp). e.g. \code{c(-1000,1000)}
#' zooms out, \code{c(1000,-1000)} zooms in and \code{c(-1000,0)} shifts the plot
#' window to the left.
#' @param nElements \code{integer} Number of elements that are displayed
#' in each track (Default=25). High n.elements can lead to overplotting. Use
#' \code{nElements=FALSE} to display all elements.
#' @return OGREDataSet.
#' @examples
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' myOGRE <- fOverlaps(myOGRE)
#' myOGRE <- gvizPlot(myOGRE,query="ENSG00000142168")
#' @export
gvizPlot <- function(OGREDataSet,query,
 gvizPlotsFolder = metadata(OGREDataSet)$gvizPlotsFolder,
 trackRegionLabels =setNames(rep("ID",length(OGREDataSet)),names(OGREDataSet)),
 trackShapes=setNames(rep("fixedArrow",length(OGREDataSet)),names(OGREDataSet)),
 showPlot=FALSE,extendPlot=c(-300,300),nElements=25){
  queryID <- subjType <-  NULL
  for(q in query){
    message(paste0("Plotting query: ",q))
    #Subject tracks
    GvizSubjTracks<-lapply(metadata(OGREDataSet)$subjectNames,function(x){
      tmp<-metadata(OGREDataSet)$detailDT[subjType==x & queryID==q]
      if(is.numeric(nElements&dim(tmp)[1]>nElements)){tmp<-tmp[sample
                            (seq_len(dim(tmp)[1]),size = nElements),]}
      regionLabels<-mcols(OGREDataSet[[x]])[[trackRegionLabels[x]]][match(
        tmp$subjID,mcols(OGREDataSet[[x]])[["ID"]])]
      if(dim(tmp)[1]!=0){# if subjectType overlaps with query create track
        if(is.null(regionLabels)){
  AnnotationTrack(start=tmp$subjStart,end = tmp$subjEnd,chromosome=tmp$subjChr,
              name=x,id = regionLabels,fontsize.title=24,shape=trackShapes[x],
              featureAnnotation=NULL,fontcolor.title="black",fontcolor="black",
              fontcolor.group="black",fontcolor.item="black",rotation.item=20,
              arrowHeadWidth=30,strand=tmp$subjStrand)
        }else{
  AnnotationTrack(start=tmp$subjStart,end = tmp$subjEnd,chromosome=tmp$subjChr,
              name=x,id = regionLabels,fontsize.title=24,shape=trackShapes[x],
              featureAnnotation="id",fontcolor.title="black",fontcolor="black",
              fontcolor.group="black",fontcolor.item="black",rotation.item=20,
              arrowHeadWidth=30,strand=tmp$subjStrand)          
        }
      }else{
  AnnotationTrack(GRanges(),name = x,fontcolor.title="black",fontsize.title=24)
      }
    })
    #Gviz helper tracks
    chr<-as.character(seqnames(OGREDataSet[[1]])[mcols(OGREDataSet[[1]])$ID==q])
    itrack <-IdeogramTrack(genome="hg38",chromosome=chr)#ideogram track (chromosome bands etc)
    gtrack<-GenomeAxisTrack()
    queryGR<-OGREDataSet[[1]][mcols(OGREDataSet[[1]])$ID==q]
    from <- start(queryGR)+extendPlot[1]  #add 300bp left and righ as plotWindow
    to <- end(queryGR)+extendPlot[2]
    #Gviz query track
    regionLabels<-mcols(queryGR)[[trackRegionLabels[1]]]
    queryTrack<-AnnotationTrack(range=queryGR,name=names(OGREDataSet)[1],fill="red",
       arrowHeadWidth=30,shape=trackShapes[1],featureAnnotation="id",
       id=regionLabels,fontsize.group=20,fontsize.title=24,
       fontcolor.title="black")#,stacking = "dense"
    allTracks<-c(itrack,gtrack,queryTrack,GvizSubjTracks)#;names(temp)=make.unique(names(temp))
    #Gviz plotting
    pdf(file.path(gvizPlotsFolder,paste0(q,".pdf")),width = 30/2.54,height = 20/2.54)
    plotTracks(allTracks,showOverplotting=TRUE,from = from, to = to,title.width=6,
               rotation.title=0,background.title="white",just.group="above")
    dev.off()
    if(showPlot){
      plotTracks(allTracks,showOverplotting=TRUE,from = from, to = to,title.width=6,
                 rotation.title=0,background.title="white",just.group="above")
    }
  }
  return(OGREDataSet)
}
#' Plot histogram
#'
#' Plots overlap histograms of all subject datasets and stores them as a list,
#' that can be accessed by \code{metadata(myOGRE)$hist}                                                                                                                                                                                           
#' @importFrom ggplot2 ggplot 
#' @param OGREDataSet An OGREDataSet
#' @param plot0 plot0=FALSE(default) plots a histogram of all dataset elements
#' with overlaps, excluding elements without overlaps. plot0=FALSE also includes
#' elements without overlaps.
#' @return OGREDataSet.
#' @examples 
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' myOGRE <- fOverlaps(myOGRE)
#' myOGRE <- plotHist(myOGRE)
#' metadata(myOGRE)$hist
#' @export
plotHist <- function(OGREDataSet,plot0=FALSE){
  for(i in colnames(metadata(OGREDataSet)$quickDT[,-c(1)])){
    dt <- metadata(OGREDataSet)$quickDT
    if(isFALSE(plot0))dt <- dt[dt[[i]]!=0,]
    p <- ggplot(dt, aes(x=.data[[i]])) +
    geom_histogram(fill="steelblue2", position="dodge")+
    geom_vline(xintercept=stats::median(dt[[i]]),
               linetype="dashed",size=1.5)+
    scale_x_continuous(trans = 'log10')+
    theme(legend.position="top")+theme_classic()+
    theme(axis.title.x=element_blank(),text = element_text(size = 20))+
    ylab(i)
    metadata(OGREDataSet)$hist[[i]] <- p
  }
  return(OGREDataSet)
}

.datatable.aware=TRUE #to use "[" without importing whole data.table
#' Coverage plot
#'
#' Generates coverage plots of all subject datasets and stores them as a list,
#' that can be accessed by \code{metadata(OGREDataSet)$covPlot}                                                                                                                                                                                           
#' @import ggplot2 
#' @rawNamespace import(data.table, except = c(shift,second,first))
#' @importFrom tidyr %>%
#' @param OGREDataSet An OGREDataSet
#' @param datasets \code{character vector} of subject dataset names. Default: 
#' Generates a coverage plots for all subjects
#' @param nbin Number of bins
#' @return OGREDataSet.
#' @examples 
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' myOGRE <- fOverlaps(myOGRE)
#' myOGRE <- covPlot(myOGRE)
#' metadata(myOGRE)$covPlot
#' @export
covPlot <- function(OGREDataSet,
                     datasets=names(OGREDataSet)[seq(2,length(OGREDataSet))],
                     nbin=100){
  message("Generating coverage plot(s), this might take a while...")
  #filter out queries<nbin
  #regions that are at least nbin nucleotides long 
  message("Excluding regions with nucleotides<nbin")
  regions <-mcols(OGREDataSet[[1]][width(OGREDataSet[[1]])>=nbin])$ID
  #& have overlaps
  regions <- regions[regions%in%metadata(OGREDataSet)$detailDT$queryID]
  for(d in datasets){
    regions <- regions[regions%in%metadata(OGREDataSet)$detailDT[subjType==d][["queryID"]]]
    cov <- GenomicRanges::coverage(OGREDataSet[[d]])
    covDT <- sapply(unique(regions),function(r){
      cor <- metadata(OGREDataSet)$detailDT[queryID==r&subjType==d,]
      dt <- data.table(rCov=cov[[cor$queryChr[1]]][seq(cor$queryStart[1],cor$queryEnd[1])]%>%as.numeric()) 
      dt[,bins:=ggplot2::cut_interval(seq(1,length(rCov)), n = nbin)]#set bins)
      dt[,sum(rCov),by=bins][["V1"]] #bin overlaps are summed up
    })
    covDT <- as.data.table(t(covDT),keep.rownames = "ID")
    covDT[,strand:=as.vector(strand(OGREDataSet[[1]]))[match(ID,mcols(OGREDataSet[[1]])$ID)]] 
    covDT_for <- subset(covDT[strand=="+",],select=-c(strand,ID))
    covDT_rev <- subset(covDT[strand=="-",],select=-c(strand,ID))
    covDT_rev <- rev(covDT_rev) #reverse minus strand
    covDT_both <- subset(covDT[strand=="*",],select=-c(strand,ID))
    covDT <- rbind(covDT_for,covDT_rev,covDT_both)
    covDT <- colSums(covDT) #column summarization for now
    #covDT <- colSums(covDT_for) #add support for forward/reverse plots
    #covDT <- colSums(covDT_rev)
    dt <- data.table(x=seq(1,length(covDT)),Coverage=covDT)
    p <- ggplot(dt, aes(x=x, y=Coverage)) + xlab("Region bins")+theme_classic()+
      theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
      geom_smooth(se=FALSE)+ylab("Overlap coverage")
    metadata(OGREDataSet)$covPlot[[d]] <- (list(plot=p,data=dt))
  }
  return(OGREDataSet)
}

#' Calculates min/max/average overlap
#' 
#' Calculates min/max/average overlap for all datasets using \code{summary()}. 
#' Results can be accessed by \code{metadata(OGREDataSet)$summaryDT} which is a 
#' \code{list()} of two \code{data.table} objects. The first one includes 
#' elements without any overlap at all and the second provides summary numbers
#' for all elements that have at least one overlap.
#' @param OGREDataSet An OGREDataSet
#' @return OGREDataSet.
#' @examples 
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' myOGRE <- fOverlaps(myOGRE)
#' myOGRE <- summarizeOverlap(myOGRE)
#' metadata(myOGRE)$summaryDT
#' @export
summarizeOverlap <- function(OGREDataSet){
  dt <- metadata(OGREDataSet)$quickDT[,-c(1)]
  metadata(OGREDataSet)$summaryDT[["includes0"]] <- do.call(cbind, lapply(dt, summary))
  dt[dt==0] <- NA
  metadata(OGREDataSet)$summaryDT[["excludes0"]] <- do.call(cbind, lapply(dt, summary))
  return(OGREDataSet)
}

#' List predefined datasets
#' 
#' Use `listPredefinedDataSets()` to receive a vector of names for predefined
#' datasets that can be aquired from AnnotationHub that are already correctly 
#' parsed and formatted. Each of the listed names can be used as input for
#' `addDataSetFromHub()`. Currently supported:
#' \itemize{
#'   \item protCodingGenes - Protein coding genes from HG19 (GRCh37) Ensembl
#'   For additional information use:
#'   `getInfoOnIds(AnnotationHub(), "AH10684")` 
#'   \item CGI - CpG islands from HG19 UCSC
#'   For additional information use:
#'   `getInfoOnIds(AnnotationHub(), "AH5086")`
#'   \item SNP - Common Single Nucleotide Polymorphism from HG19 UCSC
#'   For additional information use:
#'   `getInfoOnIds(AnnotationHub(), "AH5105")`
#'   \item TFBS - Transcription Factor Binding Sites conserved from HG19 UCSC
#'   For additional information use:
#'   `getInfoOnIds(AnnotationHub(), "AH5090")`
#'   \item Promoters - Promoter and flanking regions from HG19 Ensembl (Note: 
#'   This annotation is currently not included in AnnotationHub and is therefore
#'   downloaded from Ensembl's ftp site)
#' }
#' @return \code{character} vector.
#' @examples
#' listPredefinedDataSets()
#' @export
listPredefinedDataSets <- function(){
  return(c("protCodingGenes","CGI","SNP","TFBS","Promoters"))
}

#' Add dataSet from AnnotationHub
#' 
#' AnnotationHub offers a wide range of annotated datasets which can be manually
#' aquired but need some parsing to work with OGRE as detailed in vignette
#' section "Load datasets from AnnotationHub". 
#' For convienence `addDataSetFromHub()` adds one of the predefined human 
#' dataSets of `listPredefinedDataSets()` to a OGREDataSet.Those are taken from 
#' AnnotationHub and are ready to use for OGRE. Additional information on 
#' dataSets can be found here \code{\link{listPredefinedDataSets}}. 
#' @rawNamespace import(data.table, except = c(shift,second,first))
#' @importFrom assertthat assert_that
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom GenomeInfoDb keepStandardChromosomes seqlevelsStyle dropSeqlevels
#' @importFrom rtracklayer import.gff
#' @param OGREDataSet OGREDataSet
#' @param dataSet \code{character} Name of one predefined dataSets to add as
#' query or subject to a OGREDataSet. Possible dataSets can be show with
#' `listPredefinedDataSets()`.
#' @param type Type of dataSet, must be either query or subject. If query the
#' dataSet will be added as query and at the first position of OGREDataSet.
#' @return OGREDataSet.
#' @examples 
#' myOGRE <- OGREDataSet() 
#' myOGRE <- addDataSetFromHub(myOGRE,"protCodingGenes","query")
#' @export
addDataSetFromHub <- function(OGREDataSet,dataSet,type){
  assertthat::assert_that(type%in%c("query","subject"),
                          msg="Parameter type must be query or subject.")
  assertthat::assert_that(dataSet%in%listPredefinedDataSets(),
                          msg=paste("Parameter dataSet type must be among:",
                          paste(listPredefinedDataSets(),collapse = " ")))
  if(is.null(metadata(OGREDataSet)$aH)){aH <- AnnotationHub()}
  switch(dataSet,
         protCodingGenes={x <- aH[["AH10684"]]
         x <- GenomeInfoDb::keepStandardChromosomes(x,"Homo_sapiens",
                                                    pruning.mode="coarse")
         x <- x[mcols(x)$type=="gene"&mcols(x)$gene_biotype=="protein_coding"]
         mcols(x) <-mcols(x)[,c("gene_id","gene_name")]
         colnames(mcols(x)) <- c("ID","name")
         },
         CGI={x <- aH[["AH5086"]]
         x <- GenomeInfoDb::keepStandardChromosomes(x,"Homo_sapiens",
                                                    pruning.mode="coarse")
         GenomeInfoDb::seqlevelsStyle(x) <- "Ensembl"
         mcols(x) <- data.frame(ID=seq(1,length(x)),name=mcols(x)$name)
         },
         SNP={x <- aH[["AH5105"]]
         x <- GenomeInfoDb::keepStandardChromosomes(x,"Homo_sapiens",
                                                    pruning.mode="coarse")
         GenomeInfoDb::seqlevelsStyle(x) <- "Ensembl"
         mcols(x) <-data.frame(ID=make.unique(mcols(x)$name),name=mcols(x)$name)
         },
         TFBS={x <- aH[["AH5090"]]
         x <- GenomeInfoDb::keepStandardChromosomes(x,"Homo_sapiens",
                                                    pruning.mode="coarse")
         GenomeInfoDb::seqlevelsStyle(x) <- "Ensembl"
         mcols(x)$name <- gsub("V\\$","",mcols(x)$name)
         mcols(x)$ID <- mcols(x)$name
         mcols(x)$ID <- make.unique(gsub("_.*","",mcols(x)$ID))
         },
         Promoters={
           x <- as.data.table(rtracklayer::import.gff("http://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz"))
           x <- x[type%in%c("promoter","promoter_flanking_region"),c("seqnames","start","end","ID")]
           x[,name:=sub(".*:", "", ID)  ]
           x <- makeGRangesFromDataFrame(x,keep.extra.columns = TRUE)
           x <- GenomeInfoDb::keepStandardChromosomes(x,"Homo_sapiens",
                                                      pruning.mode="coarse")
           GenomeInfoDb::seqlevelsStyle(x) <- "Ensembl"
           }
  )
  #MT removed for compatibility with hg19,MT(hg19) differs from MT(GRCH37) see:
  #https://gatk.broadinstitute.org/hc/en-us/articles/360035890711?id=23390
  x <- GenomeInfoDb::dropSeqlevels(x,"MT","coarse") 
  genome(x) <- "hg19"
  mData <- metadata(OGREDataSet)
  if(dataSet%in%names(OGREDataSet)){
    dataSet <- tail(make.unique(c(names(OGREDataSet),dataSet)),n=1)
    warning("Renamed datasets cause of duplicated labels")
    }
  if(type=="query"){ #add dataSet as query or subject
    OGREDataSet <- c(GRangesList(x),OGREDataSet) #add x at first position
    names(OGREDataSet)[1] <- dataSet
    metadata(OGREDataSet) <- mData
  }else{
    OGREDataSet[[dataSet]] <- x
    metadata(OGREDataSet) <- mData
  if(length(OGREDataSet)>1){ #update subject names
  metadata(OGREDataSet)$subjectNames<-names(
    OGREDataSet[seq(2,length(OGREDataSet))])
  }}
  message(paste0("Added ",dataSet))
  return(OGREDataSet)
}

#' Add GenomicRanges 
#'
#' Add a GenomicRanges dataset to OGREDataSet
#' @importFrom assertthat assert_that 
#' @param OGREDataSet An OGREDataSet
#' @param dataSet A GRanges object. Each region needs chromosome, start, end and
#' strand information. A unique ID and a name column must be present in the 
#' `GenomicRanges` object metadata. Avoid different chromosome naming 
#' conventions i.e. (chr1, CHR1, 1, I) among all datasets
#' @param type Type of dataSet, must be either query or subject. If query the
#' dataSet will be added as query and at the first position of OGREDataSet. 
#' @param label A \code{character} that will label your GRanges object. If
#' not supplied, the label will be guessed from the dataset parameter.
#' @return OGREDataSet.
#' @examples 
#' myOGRE <- OGREDataSet()
#' myGRanges <- makeExampleGRanges()
#' myOGRE <- addGRanges(myOGRE,myGRanges,"query")
#' @export
addGRanges <- function(OGREDataSet,dataSet,type,label=NULL){
  assertthat::assert_that(type%in%c("query","subject"),
                          msg="Parameter type must be query or subject.")
  assertthat::assert_that(is(dataSet,"GRanges"),
                          msg="dataSet class must be GRanges.")
  assertthat::assert_that(c("ID")%in%names(mcols(dataSet)),
                          msg="DataSet must contain ID column.")
  assertthat::assert_that(!any(duplicated(dataSet$ID)),
                          msg="ID column must be unique.")
  assertthat::assert_that(length(dataSet)!=0,
                          msg="DataSet has no ranges.")
  mData <- metadata(OGREDataSet)
  if(is.null(label))label <- deparse(substitute(dataSet)) #if label not supplied
  if(label%in%names(OGREDataSet)){
    label <- tail(make.unique(c(names(OGREDataSet()),label)),n=1)
    warning("Renamed datasets cause of duplicated labels")
  }
  
  if(type=="query"){
    OGREDataSet <- c(GRangesList(dataSet),OGREDataSet)
    names(OGREDataSet)[1] <- label
    metadata(OGREDataSet) <- mData
  }else{
    OGREDataSet[[label]] <- dataSet
    metadata(OGREDataSet) <- mData
    if(length(OGREDataSet)>1){ #update subject names
      metadata(OGREDataSet)$subjectNames <-names(OGREDataSet
      [seq(2,length(OGREDataSet))])
    }  
  }
  return(OGREDataSet)
}

#' Make a example OGRE dataset
#' 
#' `makeExampleOGREDataSet` generates a example OGREDataSet from dataset files
#' stored in OGRE's extdata directory.
#' @return OGREDataSet.
#' @examples
#' myOGRE <- makeExampleOGREDataSet()
#' @export
makeExampleOGREDataSet <- function()
{
  myQueryFolder <- file.path(system.file('extdata', package = 'OGRE'),"query")
  mySubjectFolder <- file.path(system.file('extdata', package = 'OGRE'),"subject")
  myOGRE <- OGREDataSetFromDir(queryFolder=myQueryFolder,subjectFolder=mySubjectFolder)
  return(myOGRE)
}
#' Make an example GRanges dataset
#'
#' `makeExampleGRanges` generates an example GRanges dataset.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @return OGREDataSet.
#' @examples
#' myGRanges <- makeExampleGRanges()
#' @export
makeExampleGRanges <- function()
{
  myGRanges <- GRanges(Rle(c("2", "2", "1", "3"), c(1, 3, 2, 4)),
             IRanges(seq_len(10), width=seq(10,1), names=head(letters, 10)),
             Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
             ID=seq_len(10), name=paste0("gene",seq_len(10)))
  return(myGRanges)
}
#' Subset a GRanges object
#'
#' Subsets a GRanges object with reference to it's ID column using a ID vector.
#' @importFrom GenomicRanges GRanges
#' @param OGREDataSet An OGREDataSet
#' @param IDs \code{character vector} with IDs used to subset the GRanges object
#' defined in \code{name}
#' @param name \code{character} Name of the GRanges object for subsetting. One of 
#' the GRanges objects in a \code{OGREDataSet}
#' @return OGREDataSet.
#' @examples
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' myOGRE <- subsetGRanges(myOGRE,c("ENSG00000142168","ENSG00000256715"),"genes")
#' @export
subsetGRanges <- function(OGREDataSet,IDs,name)
{
  OGREDataSet[[name]] <- OGREDataSet[[name]][mcols(OGREDataSet[[name]])$ID%in%IDs,]
  return(OGREDataSet)
}

#' Extend a GRanges object
#'
#' Extend(shrink) ranges of a GRanges object.
#' @importFrom GenomicRanges GRanges
#' @param OGREDataSet  An OGREDataSet
#' @param name \code{character} Name of the GRanges object for extending 
#' @param upstream \code{int} (positive or negative number)
#' @param downstream \code{int} (positive or negative number)
#' @return OGREDataSet
#' @examples
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' #extend range by shifting start 100 bp in upstream direction
#' myOGRE <- extendGRanges(myOGRE,"genes",upstream=100)
#' #shrinking range by shifting end 100 bp in upstream direction
#' myOGRE <- extendGRanges(myOGRE,"genes",downstream=-100)
#' #shrinking range by shifting from both sides to the center
#' myOGRE <- extendGRanges(myOGRE,"genes",upstream=-10,downstream=-10)
#' @export
extendGRanges <- function(OGREDataSet,name, upstream=0, downstream=0)     
  #taken from https://support.bioconductor.org/p/78652/
{
  if (any(strand(OGREDataSet[[name]]) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(OGREDataSet[[name]]) == "+" | strand(OGREDataSet[[name]]) == "*"
  new_start <- start(OGREDataSet[[name]]) - ifelse(on_plus, upstream, downstream)
  new_end <- end(OGREDataSet[[name]]) + ifelse(on_plus, downstream, upstream)
  ranges(OGREDataSet[[name]]) <- IRanges(new_start, new_end)
  OGREDataSet[[name]] <- trim(OGREDataSet[[name]])
  return(OGREDataSet)
}

#' Extract promoter
#'
#' A wrapper of \code{GenomicRanges::promoters()} to extract promoter regions of 
#' a GRanges object stored in a OGREDataSet
#' @importFrom GenomicRanges GRanges
#' @param OGREDataSet  An OGREDataSet
#' @param name \code{character} Name of the GRanges object 
#' @param upstream \code{int} (positive) upstream=2000(default)
#' @param downstream \code{int} (positive) downstream=200(default)
#' @return OGREDataSet
#' @examples
#' myOGRE <- makeExampleOGREDataSet()
#' myOGRE <- loadAnnotations(myOGRE)
#' myOGRE <- extractPromoters(myOGRE,"genes", upstream=2000, downstream=200)
#' @export
extractPromoters <- function(OGREDataSet,name, upstream=2000, downstream=200)     
  #taken from https://support.bioconductor.org/p/78652/
{
  if (any(strand(OGREDataSet[[name]]) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(OGREDataSet[[name]]) == "+" | strand(OGREDataSet[[name]]) == "*"
  new_start <- start(OGREDataSet[[name]]) - ifelse(on_plus, upstream, downstream)
  new_end <- end(OGREDataSet[[name]]) + ifelse(on_plus, downstream, upstream)
  ranges(OGREDataSet[[name]]) <- IRanges(new_start, new_end)
  OGREDataSet[[name]] <- trim(OGREDataSet[[name]])
  return(OGREDataSet)
}
