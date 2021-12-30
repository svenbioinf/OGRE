#' Load annotation datasets
#'
#' Load dataset files containing genomic regions annotation information from hard drive. `loadAnnotations` calls `readQuery` and `readSubject`
#' to read in genomic regions as `GenomicRanges` objects stored as .RDS / .rds files. Each region needs chromosome, start, end and strand information.
#' A unique ID and a name column must be present in the `GenomicRanges` object metadata. OGRE searches for the query file in your query folder and any number of subject files in your subjects folder.
#' @param OGREDataSet A OGREDataSet.
#' @return A OGREDataSet.
#' @examples
#' myOGRE=makeExampleOGREDataSet()
#' myOGRE=loadAnnotations(myOGRE)
#' @export

loadAnnotations <- function(OGREDataSet){
  OGREDataSet <- readQuery(OGREDataSet)
  OGREDataSet <- readSubject(OGREDataSet)

 return(OGREDataSet)
}




#' Read query dataset
#'
#'[readQuery()] scanns `queryFolder` for a `GRanges` object stored as .RDS / .rds file and attaches it to the OGREDataSet.
#' @param OGREDataSet A OGREDataSet.
#' @return A OGREDataSet.
readQuery=function(OGREDataSet){
  assertthat::assert_that(length(list.files(metadata(OGREDataSet)$queryFolder))>0,
                          msg="Query folder is empty!")
  message("Reading query dataset... ")
  if(is.null(metadata(OGREDataSet)$queryFolder)){ #if no folder supplied, check default folder

  }else{ #read queryFolder
    queryPath <- list.files(metadata(OGREDataSet)$queryFolder,full.names = TRUE)
    queryName <- tools::file_path_sans_ext(list.files(metadata(OGREDataSet)$queryFolder))
    OGREDataSet[[queryName]]<- readRDS(queryPath)
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
#'[readSubject()] scanns `SubjectFolder` for any `GRanges` objects stored as .RDS / .rds files and attaches them to the OGREDataSet.
#' @param OGREDataSet A OGREDataSet.
#' @return A OGREDataSet.
readSubject=function(OGREDataSet){
  assertthat::assert_that(length(list.files(metadata(OGREDataSet)$subjectFolder))>0,
                          msg="Subject folder is empty!")
  message("Reading subject datasets... ")
  #read queryFolder
    subjectPath <- list.files(metadata(OGREDataSet)$subjectFolder,full.names = TRUE)
    subjectName <- tools::file_path_sans_ext(list.files(metadata(OGREDataSet)$subjectFolder))
    OGREDataSet <- c(OGREDataSet,mapply(function(x,y){
      tmp=readRDS(y)
      assertthat::assert_that(c("ID")%in%names(mcols(tmp)),msg="Subject must contain ID column.")
      assertthat::assert_that(!any(duplicated(tmp$ID)),msg="ID column must be unique.")
      assertthat::assert_that(length(tmp)!=0,msg=paste0("Dataset has no ranges: ",x))

      tmp
    },x=subjectName,y=subjectPath))

  return(OGREDataSet)

}




#' Find overlaps
#'
#' Finds all overlaps between query and subject(s) and stores each hit (overlap) in data table `detailDT`. Data table `sumDT` shows all
#' overlaps of a certain subject type for all query elements. By default also partially overlaps are reported.
#'
#' @param OGREDataSet A OGREDataSet.
#' @param selfHits \code{logical} if FALSE(default) ignores self hits of identical regions within datasets. 
#' @return OGREDataSet.
#' @examples
#' myOGRE=makeExampleOGREDataSet()
#' myOGRE=loadAnnotations(myOGRE)
#' myOGRE=fOverlaps(myOGRE)
#' @export
fOverlaps <- function(OGREDataSet,selfHits=FALSE){
  detailDT <- data.table() #data table to store all overlaps for query vs all subjects
  metadata(OGREDataSet)$subjectNames <- names(OGREDataSet[2:length(OGREDataSet)])
  for(subj in metadata(OGREDataSet)$subjectNames){
    ol <- findOverlaps(OGREDataSet[[1]],OGREDataSet[[subj]])
    ol <- data.table(q=queryHits(ol),s=mcols(OGREDataSet[[subj]])[subjectHits(ol),1])
    if(length(ol)==0){ #if no overlap hits skip to next subj
      message("No overlap found for: ",subj)
      next
    }
    else{
      detailDT <- rbind(detailDT,data.table(
       queryID=mcols(OGREDataSet[[1]])[ol$q,"ID"],
       queryType=names(OGREDataSet)[1],
       subjID=ol$s,
       subjType=subj,
       queryChr=as.character(seqnames(OGREDataSet[[names(OGREDataSet)[1]]]))[ol$q],
       queryStart=start(OGREDataSet[[names(OGREDataSet)[1]]])[ol$q],
       queryEnd=end(OGREDataSet[[names(OGREDataSet)[1]]])[ol$q],
       subjChr=as.character(seqnames(OGREDataSet[[subj]]))
       [match(ol$s,mcols(OGREDataSet[[subj]])[,"ID"])],
       subjStart=start(OGREDataSet[[subj]])
       [match(ol$s,mcols(OGREDataSet[[subj]])[,"ID"])],
       subjEnd=end(OGREDataSet[[subj]])
       [match(ol$s,mcols(OGREDataSet[[subj]])[,"ID"])]
      ))
      if(isFALSE(selfHits)){ #remove self hits
        detailDT <- detailDT[!(queryID==subjID)]
      }
    }
  }
  metadata(OGREDataSet)$sumDT <- detailDT[,list(
    subjType=subjType,
    subjID=paste(subjID,collapse = " ")),by="queryID"]
  assert_that(nrow(metadata(OGREDataSet)$sumDT)!=0,
              msg="No overlaps found for any subjects!")
  metadata(OGREDataSet)$detailDT <- detailDT
  return(OGREDataSet)
}

#' Generate summary plot
#'
#' `sumPlot()` calculates key numbers i.e. (total number of overlaps, number of overlaps per subject...) to help with an
#' exploratory data evaluation and displays them in an informative barplot.
#' @importFrom AnnotationHub query
#' @param OGREDataSet A OGREDataSet.
#' @return OGREDataSet.
#' @examples
#' myOGRE=makeExampleOGREDataSet()
#' myOGRE=loadAnnotations(myOGRE)
#' myOGRE=fOverlaps(myOGRE)
#' myOGRE=sumPlot(myOGRE)
#' @export
sumPlot <- function(OGREDataSet){
  assert_that(!is.null(metadata(OGREDataSet)$sumDT),
              msg="No data to generate sumPlot!")
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
  metadata(OGREDataSet)$barplot_summary_dt <- dtbarplot

  temp<-sapply(metadata(OGREDataSet)$subjectNames,function(x){
    length(unique(metadata(OGREDataSet)$detailDT[
      metadata(OGREDataSet)$detailDT$subjType==x,][["queryID"]]))})
  temp<-data.table(x=names(temp),
   Subjects=temp,
   lab=paste0("Queries with ",names(temp),": ",temp," (",round(temp/query_N*100),"%)"),
   col="steelblue3")
  subjPerQuery<-round(sapply(metadata(OGREDataSet)$subjectNames,function(x){
    sum(metadata(OGREDataSet)$detailDT$subjType==x)})/query_N,digits=1)
  subjPerQuery<-data.table(x=names(subjPerQuery),
                           Subjects=subjPerQuery,
  lab=paste0("Average ",names(subjPerQuery)," per query: ",subjPerQuery),col="steelblue4")

  dtbarplotDetailed<-rbind(dtbarplot,temp,subjPerQuery, use.names=FALSE)
  dtbarplotDetailed$sequence<-dim(dtbarplotDetailed)[1]:1
  dtbarplotDetailed$xtext<-dtbarplotDetailed$sequence+0.55
  dtbarplotDetailed$query<-log2(dtbarplotDetailed$query+1)#add +1 to avoid log2(1)=0 and log2(0)=-inf
  metadata(OGREDataSet)$barplot_summary_dt <- dtbarplotDetailed
  metadata(OGREDataSet)$barplot_summary <-
    ggplot(dtbarplotDetailed, aes(x = sequence, y =query,fill=col)) +
      geom_bar(stat = "identity", width = 0.3) +
      scale_fill_manual(guide=FALSE,values = c("steelblue2" = "steelblue2",
                                               "steelblue3" = "steelblue3",
                                               "steelblue4" = "steelblue4"))+
      coord_flip() +labs(x = "", y = "Log2(N)") +theme_bw() +  theme_classic() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      geom_text(aes(x = xtext, y = 0, label = lab),size=7,hjust = 0, vjust = 1)+
    guides(size = FALSE)
    ggsave(plot=metadata(OGREDataSet)$barplot_summary,
           path=metadata(OGREDataSet)$outputFolder,width = 30,
           filename="Barplot_Summary.png",height = 20,dpi = 400,units = "cm")
    return(OGREDataSet)
}

#' Generate Gviz plot
#'
#'`gvizPlot` generates a plot around one or many given query elements with all overlapping subject hits. In addition,
#'each generated plot is stored in the `gvizPlots` folder next to you
#' @importFrom grDevices dev.off pdf
#' @importFrom stats setNames
#' @param OGREDataSet A OGREDataSet.
#' @param query A character vector of one or many query elements ID's (i.e. Gene ID's).
#' @param gvizPlotsFolder A character pointing to the plot(s) output directory. If not supplied a folder is
#' generated next to the dataset folder.
#' @param showPlot \code{logical} If \code{FALSE}(default) plots are only saved to `gvizPlotsFolder`. If \code{TRUE}
#' plots are additionally send to the plotting window.
#' @param trackRegionLabels A labeled character vector that defines the type of label that is displayed for subject elements.
#' Values represent subject types and and vector labels define the type of label of each subject element. Value "ID" and label "genes" would
#' use the "ID" column of your gene dataset as labels.
#' @return OGREDataSet.
#' @examples
#' myOGRE=makeExampleOGREDataSet()
#' myOGRE=loadAnnotations(myOGRE)
#' myOGRE=fOverlaps(myOGRE)
#' myOGRE=gvizPlot(myOGRE,query="ENSMUSG00000068196")
#' @export
gvizPlot <- function(OGREDataSet,query,
 gvizPlotsFolder = metadata(OGREDataSet)$gvizPlotsFolder,
 trackRegionLabels =setNames(rep("ID",length(OGREDataSet)),names(OGREDataSet)),
 showPlot=FALSE){
  for(q in query){
    #Subject tracks
    GvizSubjTracks<-lapply(metadata(OGREDataSet)$subjectNames,function(x){
      tmp<-metadata(OGREDataSet)$detailDT[subjType==x & queryID==q]
      if(dim(tmp)[1]>25){tmp<-tmp[sample(1:dim(tmp)[1],size = 25),]}
      regionLabels<-mcols(OGREDataSet[[x]])[[trackRegionLabels[x]]][match(
        tmp$subjID,mcols(OGREDataSet[[x]])[["ID"]])]
        if(dim(tmp)[1]!=0){# if subjectType overlaps with query create track
          AnnotationTrack(start=tmp$subjStart,end = tmp$subjEnd,chr=tmp$subjChr,
                        name=x,id = regionLabels,fontsize.title=24,
                        featureAnnotation="id",fontcolor.title="black",fontcolor="black",
                        fontcolor.group="black",fontcolor.item="black",rotation.item=20)
        }else{
          AnnotationTrack(GRanges(),name = x,fontcolor.title="black",fontsize.title=24)
        }
    })
    #Gviz helper tracks
    chr<-as.character(seqnames(OGREDataSet[[1]])[mcols(OGREDataSet[[1]])$ID==q])
    itrack <-IdeogramTrack(genome="hg38",chromosome=chr)#ideogram track (chromosome bands etc)
    gtrack<-GenomeAxisTrack()
    queryGR<-OGREDataSet[[1]][mcols(OGREDataSet[[1]])$ID==q]
    from <- start(queryGR)-300  #add 300bp left and righ as plotWindow
    to <- end(queryGR)+300
    plotWindow<-GRanges(seqnames=chr,ranges=IRanges::IRanges(start = from, end = to))
    #Gviz query track
    regionLabels<-mcols(queryGR)[[trackRegionLabels[1]]]
    queryTrack<-AnnotationTrack(range=queryGR,name="Query",fill="red",arrowHeadWidth=30,
                               shape="fixedArrow",featureAnnotation="id",
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

#' List predefined datasets
#' 
#' Use `listPredifinedDataSets()` to receive a vector of names for predefined
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
#' }
#' @return \code{character} vector.
#' @examples
#' listPredefinedDataSets()
#' @export
listPredefinedDataSets <- function(){
  return(c("protCodingGenes","CGI"))
}

#' Add dataSet from AnnotationHub
#' AnnotationHub offers a wide range of annotated datasets which can be manually
#' aquired but need some parsing to work with OGRE as detailed in vignette
#' section "Access to annotation data". 
#' For convienence `addDataSetFromHub()` adds one of the predefined human 
#' dataSets of `listPredefinedDataSets()` to a OGREDataSet.Those are taken from 
#' AnnotationHub and are ready to use for OGRE. Additional information on 
#' dataSets can be found here \code{\link{listPredefinedDataSets}}. 
#' @importFrom assertthat assert_that
#' @importFrom AnnotationHub AnnotationHub
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
         x <- x[mcols(x)$type=="gene"&mcols(x)$gene_biotype=="protein_coding"]
         mcols(x) <-mcols(x)[,c("gene_id","gene_name")]
         colnames(mcols(x)) <- c("ID","name")
         },
         CGI={x <- aH[["AH5086"]]
         x <- keepStandardChromosomes(x,"Homo_sapiens",pruning.mode="coarse")
         seqlevelsStyle(x) <- "Ensembl"
         mcols(x) <- data.frame(ID=seq(1,length(x)),name=mcols(x)$name)
         },
         SNP={x <- aH[["AH5105"]]
         x <- keepStandardChromosomes(x,"Homo_sapiens",pruning.mode="coarse")
         seqlevelsStyle(x) <- "Ensembl"
         mcols(x) <- data.frame(ID=make.unique(mcols(x)$name),name=mcols(x)$name)
         }       
         
  )
  if(type=="query"){ #add dataSet as query or subject
    OGREDataSet <- c(GRangesList(x),OGREDataSet) #add x at first position
    names(OGREDataSet)[1] <- dataSet
  }else{
    OGREDataSet[[dataSet]] <- x
  }
  if(length(OGREDataSet)>1){ #update subject names
  metadata(OGREDataSet)$subjectNames<-names(OGREDataSet[2:length(OGREDataSet)])
  }
  return(OGREDataSet)
}

#' Add a GenomicRanges dataset to OGREDataSet                                                                                                                                                                                             
#' @importFrom assertthat assert_that 
#' @param OGREDataSet An OGREDataSet
#' @param dataSet A GRanges object. Each region needs chromosome, start, end and
#' strand information. A unique ID and a name column must be present in the 
#' `GenomicRanges` object metadata.
#' @param type Type of dataSet, must be either query or subject. If query the
#' dataSet will be added as query and at the first position of OGREDataSet. 
#' @return OGREDataSet.
#' @examples 
#' myOGRE <- OGREDataSet()
#' myGRanges <- makeExampleGRanges()
#' myOGRE <- addGRanges(myOGRE,myGRanges,"query")
#' @export
addGRanges <- function(OGREDataSet,dataSet,type){
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
  if(type=="query"){
    #OGREDataSet[[1]] <- c(GRangesList(dataSet),OGREDataSet)
    OGREDataSet[[1]] <- dataSet
    names(OGREDataSet)[1] <- deparse(substitute(dataSet))
    
  }else{
    OGREDataSet[[deparse(substitute(dataSet))]] <- dataSet
    metadata(OGREDataSet)$subjectNames <-names(OGREDataSet
    [2:length(OGREDataSet)])
  }
  return(OGREDataSet)
}

#' Make a example OGRE dataset
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
  myGRanges <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
                       IRanges(1:10, width=10:1, names=head(letters, 10)),
                       Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
                       ID=1:10, name=paste0("gene",1:10))
  return(myGRanges)
}