#' Load annotation datasets
#'
#' Load dataset files containing genomic regions annotation information from hard drive. `loadAnnotations` calls `readQuery` and `readSubject`
#' to read in genomic regions as `GenomicRanges` objects stored as .RDS / .rds files. Each region needs chromosome, start, end and strand information.
#' A unique ID and a name column must be present in the `GenomicRanges` object metadata. TREGEL searches for the query file in your query folder and any number of subject files in your subjects folder.
#' @param TREGELDataSet A TREGELDataSet.
#' @return A TREGELDataSet.
#' @export

loadAnnotations <- function(TREGELDataSet){
  TREGELDataSet <- readQuery(TREGELDataSet)
  TREGELDataSet <- readSubject(TREGELDataSet)

 return(TREGELDataSet)
}




#' Read query dataset
#'
#'[readQuery()] scanns `queryFolder` for a `GRanges` object stored as .RDS / .rds file and attaches it to the TREGELDataSet.
#' @param TREGELDataSet A TREGELDataSet.
#' @return A TREGELDataSet.
readQuery=function(TREGELDataSet){
  assertthat::assert_that(length(list.files(metadata(TREGELDataSet)$queryFolder))>0,msg="Query folder is empty!")
  message("Reading query dataset... ")
  if(is.null(metadata(TREGELDataSet)$queryFolder)){ #if no folder supplied, check default folder

  }else{ #read queryFolder
    queryPath <- list.files(metadata(TREGELDataSet)$queryFolder,full.names = TRUE)
    queryName <- tools::file_path_sans_ext(list.files(metadata(TREGELDataSet)$queryFolder))
    TREGELDataSet[[queryName]]<- readRDS(queryPath)
  }
  assertthat::assert_that(c("ID")%in%names(mcols(TREGELDataSet[[queryName]])),msg="Query must contain ID column.")
  assertthat::assert_that(!any(duplicated(TREGELDataSet[[queryName]]$ID)),msg="ID column must be unique.")
  assertthat::assert_that(length(TREGELDataSet[[queryName]])!=0,msg=paste0("Dataset has no ranges: ",queryName))

  return(TREGELDataSet)
}
#' Read subject datasets
#'
#'[readSubject()] scanns `SubjectFolder` for any `GRanges` objects stored as .RDS / .rds files and attaches them to the TREGELDataSet.
#' @param TREGELDataSet A TREGELDataSet.
#' @return A TREGELDataSet.
readSubject=function(TREGELDataSet){
  assertthat::assert_that(length(list.files(metadata(TREGELDataSet)$subjectFolder))>0,msg="Subject folder is empty!")
  message("Reading subject datasets... ")
  #read queryFolder
    subjectPath <- list.files(metadata(TREGELDataSet)$subjectFolder,full.names = TRUE)
    subjectName <- tools::file_path_sans_ext(list.files(metadata(TREGELDataSet)$subjectFolder))
    TREGELDataSet <- c(TREGELDataSet,mapply(function(x,y){
      tmp=readRDS(y)
      assertthat::assert_that(c("ID")%in%names(mcols(tmp)),msg="Query must contain ID column.")
      assertthat::assert_that(!any(duplicated(tmp$ID)),msg="ID column must be unique.")
      assertthat::assert_that(length(tmp)!=0,msg=paste0("Dataset has no ranges: ",x))

      tmp
    },x=subjectName,y=subjectPath))

  return(TREGELDataSet)

}




#' Find overlaps
#'
#' Finds all overlaps between query and subject(s) and stores each hit (overlap) in data table `detailDT`. Data table `sumDT` shows all
#' overlaps of a certain subject type for all query elements. By default also partially overlaps are reported.
#'
#' @param TREGELDataSet A TREGELDataSet.
#' @return TREGELDataSet.
#' @export
fOverlaps <- function(TREGELDataSet){
  detailDT <- data.table() #data table to store all overlaps for query vs all subjects
  sumDT <- data.table()
  s <- NULL
  metadata(TREGELDataSet)$subjectNames <- names(TREGELDataSet[2:length(TREGELDataSet)])
  for(subj in metadata(TREGELDataSet)$subjectNames){
    ol <- findOverlaps(TREGELDataSet[[1]],TREGELDataSet[[subj]])
    ol <- data.table(q=queryHits(ol),s=mcols(TREGELDataSet[[subj]])[subjectHits(ol),1])
    if(length(ol)==0){ #if no overlap hits skip to next subj
      message("No overlap found for: ",subj)
      next
    }
    else{
      detailDT <- rbind(detailDT,data.table(queryID=mcols(TREGELDataSet[[1]])[ol$q,"ID"],
                                     subjID=ol$s,
                                     subjType=subj,
                                     subjChr=as.character(seqnames(TREGELDataSet[[subj]]))[match(ol$s,mcols(TREGELDataSet[[subj]])[,"ID"])],
                                     subjStart=start(TREGELDataSet[[subj]])[match(ol$s,mcols(TREGELDataSet[[subj]])[,"ID"])],
                                     subjEnd=end(TREGELDataSet[[subj]])[match(ol$s,mcols(TREGELDataSet[[subj]])[,"ID"])]
      ))
      ol <- ol[,list(s=paste(s,collapse = " ")),by="q"]
      ol$q <- mcols(TREGELDataSet[[1]])[ol$q,"ID"]
      ol$subjType <- subj
      names(ol) <- c("queryID","subjID","subjType")
      sumDT <- rbind(sumDT,ol)
    }
  }
  assert_that(nrow(sumDT)!=0,msg="No overlaps found for any subjects!")
  metadata(TREGELDataSet)$detailDT <- detailDT
  metadata(TREGELDataSet)$sumDT <- sumDT
  return(TREGELDataSet)
}

#' Generate summary plot
#'
#' `sumPlot()` calculates key numbers i.e. [total number of overlaps, number of overlaps per subject...] to help with an
#' exploratory data evaluation and displays them in an informative barplot.
#' @param TREGELDataSet A TREGELDataSet.
#' @return TREGELDataSet.
#' @export
sumPlot <- function(TREGELDataSet){
  assert_that(!is.null(metadata(TREGELDataSet)$sumDT),msg="No data to generate sumPlot!")
  #summary barplot
  query_N <- length(unique(mcols(TREGELDataSet[[1]])$ID))
  query_ol <- length(unique(metadata(TREGELDataSet)$sumDT$queryID))
  query_No_ol <- query_N-query_ol
  query_full <- sum(table(metadata(TREGELDataSet)$sumDT$queryID)>=length(TREGELDataSet)-1)

  dtbarplot<-data.table(x=c("query_N","query_w_minOne_overlap","genes_w_overlap_allSubjTypes"),query=c(query_N,query_ol,query_full))
  dtbarplot$x<-factor(x = dtbarplot$x,levels = c("query_N","query_w_minOne_overlap","genes_w_overlap_allSubjTypes"))
  dtbarplot$lab<-c(paste0("Total number of submitted queries: ",query_N," (100%)"),
                  paste0("Queries with >=1 subject overlaps: ",query_ol," (",round(query_ol/query_N*100),"%)"),
                  paste0("Queries with overlaps in all subject types: ",query_full," (",round(query_full/query_N*100),"%)"))
  dtbarplot$col<-"steelblue2"
  metadata(TREGELDataSet)$barplot_summary_dt <- dtbarplot

  temp<-sapply(metadata(TREGELDataSet)$subjectNames,function(x){length(unique(metadata(TREGELDataSet)$detailDT[metadata(TREGELDataSet)$detailDT$subjType==x,][["queryID"]]))})
  temp<-data.table(x=names(temp),Subjects=temp,lab=paste0("Queries with ",names(temp),": ",temp," (",round(temp/query_N*100),"%)"),col="steelblue3")
  subjPerQuery<-round(sapply(metadata(TREGELDataSet)$subjectNames,function(x){sum(metadata(TREGELDataSet)$detailDT$subjType==x)})/query_N,digits=1)
  subjPerQuery<-data.table(x=names(subjPerQuery),Subjects=subjPerQuery,lab=paste0("Average ",names(subjPerQuery)," per subject: ",subjPerQuery),col="steelblue4")

  dtbarplotDetailed<-rbind(dtbarplot,temp,subjPerQuery, use.names=FALSE)
  dtbarplotDetailed$sequence<-dim(dtbarplotDetailed)[1]:1
  dtbarplotDetailed$xtext<-dtbarplotDetailed$sequence+0.55
  dtbarplotDetailed$query<-log2(dtbarplotDetailed$query+1)#add +1 to avoid log2(1)=0 and log2(0)=-inf
  metadata(TREGELDataSet)$barplot_summary_dt <- dtbarplotDetailed
  metadata(TREGELDataSet)$barplot_summary <- ggplot(dtbarplotDetailed, aes(x = sequence, y =query,fill=col)) +
      geom_bar(stat = "identity", width = 0.3) +
      scale_fill_manual(guide=FALSE,values = c("steelblue2" = "steelblue2", "steelblue3" = "steelblue3", "steelblue4" = "steelblue4"))+
      coord_flip() +labs(x = "", y = "Log2(N)") +theme_bw() +  theme_classic() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      geom_text(aes(x = xtext, y = 0, label = lab),size=7,hjust = 0, vjust = 1) +guides(size = FALSE)
    ggsave(plot=metadata(TREGELDataSet)$barplot_summary,path=metadata(TREGELDataSet)$outputFolder,width = 30,filename="Barplot_Summary.png",height = 20,dpi = 400,units = "cm")
    return(TREGELDataSet)
}

#' Generate Gviz plot
#'
#'`gvizPlot` generates a plot around one or many given query elements with all overlapping subject hits. In addition,
#'each generated plot is stored in the `gvizPlots` folder next to you
#'
#' @param TREGELDataSet A TREGELDataSet.
#' @param query A character vector of one or many query elements ID's (i.e. Gene ID's).
#' @param gvizPlotsFolder A character pointing to the plot(s) output directory. If not supplied a folder is
#' generated next to the dataset folder.
#' @param showPlot \code{logical} If \code{FALSE}(default) plots are only saved to `gvizPlotsFolder`. If \code{TRUE}
#' plots are additionally send to the plotting window.
#' @param trackRegionLabels A labeled character vector that defines the type of label that is displayed for subject elements.
#' Values represent subject types and and vector labels define the type of label of each subject element. Value "ID" and label "genes" would
#' use the "ID" column of your gene dataset as labels.
#' @return TREGELDataSet.
#' @export
gvizPlot <- function(TREGELDataSet,query,
                     gvizPlotsFolder = metadata(TREGELDataSet)$gvizPlotsFolder,
                     trackRegionLabels =setNames(rep("ID",length(TREGELDataSet)),names(TREGELDataSet)),
                     showPlot=FALSE){
  for(q in query){
    #Subject tracks
    GvizSubjTracks<-lapply(metadata(TREGELDataSet)$subjectNames,function(x){
      tmp<-metadata(TREGELDataSet)$detailDT[subjType==x & queryID==q]
      if(dim(tmp)[1]>25){tmp<-tmp[sample(1:dim(tmp)[1],size = 25),]}
      regionLabels<-mcols(TREGELDataSet[[x]])[[trackRegionLabels[x]]][match(tmp$subjID,mcols(TREGELDataSet[[x]])[["ID"]])]
        if(dim(tmp)[1]!=0){# if subjectType overlaps with query create track
          AnnotationTrack(start=tmp$subjStart,end = tmp$subjEnd,chr=tmp$subjChr,name=x,id = regionLabels,
                        fontsize.title=24,featureAnnotation="id",fontcolor.title="black",fontcolor="black",
                        fontcolor.group="black",fontcolor.item="black",rotation.item=20)
        }else{
          AnnotationTrack(GRanges(),name = x,fontcolor.title="black",fontsize.title=24)
        }
    })
    #Gviz helper tracks
    chr<-as.character(seqnames(TREGELDataSet[[1]])[mcols(TREGELDataSet[[1]])$ID==q])
    itrack <-IdeogramTrack(genome="hg38",chromosome=chr)#ideogram track (chromosome bands etc)
    gtrack<-GenomeAxisTrack()
    queryGR<-TREGELDataSet[[1]][mcols(TREGELDataSet[[1]])$ID==q]
    from <- start(queryGR)-300  #add 300bp left and righ as plotWindow
    to <- end(queryGR)+300
    plotWindow<-GRanges(seqnames=chr,ranges=IRanges::IRanges(start = from, end = to))
    #Gviz query track
    regionLabels<-mcols(queryGR)[[trackRegionLabels[1]]]
    queryTrack<-AnnotationTrack(range=queryGR,name="Query",fill="red",arrowHeadWidth=30,shape="fixedArrow",featureAnnotation="id",
                               id=regionLabels,fontsize.group=20,fontsize.title=24,fontcolor.title="black")#,stacking = "dense"
    allTracks<-c(itrack,gtrack,queryTrack,GvizSubjTracks)#;names(temp)=make.unique(names(temp))
    #Gviz plotting
    pdf(file.path(gvizPlotsFolder,paste0(q,".pdf")),width = 30/2.54,height = 20/2.54)
    plotTracks(allTracks,showOverplotting=TRUE,from = from, to = to,title.width=6,rotation.title=0,background.title="white",just.group="above")
    dev.off()
    if(showPlot){
      plotTracks(allTracks,showOverplotting=TRUE,from = from, to = to,title.width=6,rotation.title=0,background.title="white",just.group="above")
    }

  }
  return(TREGELDataSet)


}


