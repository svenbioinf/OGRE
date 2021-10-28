#' loadAnnotations
#'
#' @param TREGELDataSet A TREGELDataSet.
#' @return TREGELDataSet.
#' @export

loadAnnotations <- function(TREGELDataSet){
  TREGELDataSet <- readQuery(TREGELDataSet)
  TREGELDataSet <- readSubject(TREGELDataSet)

 return(TREGELDataSet)
}




#' Read query
readQuery=function(TREGELDataSet){
  message("Reading query dataset... ")
  if(is.na(metadata(TREGELDataSet)$queryFolder)){ #if no folder supplied, check default folder
    #query(obj) <- readRDS()
    #colnames(mcols(query(obj)))=c("name","id")
    #query(obj)<- makeGRangesFromDataFrame(query(obj),keep.extra.columns = TRUE)
    #obj <- readQuery(obj)

  }else{ #read queryFolder
    queryPath <- list.files(metadata(TREGELDataSet)$queryFolder,full.names = TRUE)
    queryName <- tools::file_path_sans_ext(list.files(metadata(TREGELDataSet)$queryFolder))
    TREGELDataSet[[queryName]]<- readRDS(queryPath)
    #query(obj)<- makeGRangesFromDataFrame(query(obj),keep.extra.columns = TRUE)
    #colnames(mcols(query(obj)))=c("name","id")
  }
  assertthat::assert_that(c("ID")%in%names(mcols(TREGELDataSet[[queryName]])),msg="Query must contain ID column.")
  assertthat::assert_that(!any(duplicated(TREGELDataSet[[queryName]]$ID)),msg="ID column must be unique.")

  return(TREGELDataSet)
}

#' Read subjects
readSubject=function(TREGELDataSet){
  message("Reading subjet datasets... ")
  if(is.na(metadata(TREGELDataSet)$subjectFolder)){ #if no folder supplied, check default folder

    #query(obj) <- readRDS()
    #colnames(mcols(query(obj)))=c("name","id")
    #query(obj)<- makeGRangesFromDataFrame(query(obj),keep.extra.columns = TRUE)
    #obj <- readQuery(obj)

  }else{ #read queryFolder
    subjectPath <- list.files(metadata(TREGELDataSet)$subjectFolder,full.names = TRUE)
    subjectName <- tools::file_path_sans_ext(list.files(metadata(TREGELDataSet)$subjectFolder))
    TREGELDataSet <- c(TREGELDataSet,mapply(function(x,y){
      tmp=readRDS(y)
      assertthat::assert_that(c("ID")%in%names(mcols(tmp)),msg="Query must contain ID column.")
      assertthat::assert_that(!any(duplicated(tmp$ID)),msg="ID column must be unique.")
      tmp
    },x=subjectName,y=subjectPath))
  }
  return(TREGELDataSet)

}





#' Find overlaps for query and subject
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
  metadata(TREGELDataSet)$detailDT <- detailDT
  metadata(TREGELDataSet)$sumDT <- sumDT
  return(TREGELDataSet)
}

#' Summary barplot
#' @param TREGELDataSet A TREGELDataSet.
#' @return TREGELDataSet.
#' @export
sumPlot <- function(TREGELDataSet){
  #summary barplot
  query_N <- length(unique(mcols(TREGELDataSet[[1]])$ID))
  query_ol <- length(unique(metadata(TREGELDataSet)$sumDT$queryID))
  query_No_ol <- query_N-query_ol
  query_full <- sum(table(metadata(TREGELDataSet)$sumDT$queryID)>=length(TREGELDataSet)-1)

  dtbarplot=data.table(x=c("query_N","query_w_minOne_overlap","genes_w_overlap_allSubjTypes"),query=c(query_N,query_ol,query_full))
  dtbarplot$x=factor(x = dtbarplot$x,levels = c("query_N","query_w_minOne_overlap","genes_w_overlap_allSubjTypes"))
  dtbarplot$lab=c(paste0("Total number of submitted queries: ",query_N," (100%)"),
                  paste0("Queries with >=1 subject overlaps: ",query_ol," (",round(query_ol/query_N*100),"%)"),
                  paste0("Queries with overlaps in all subject types: ",query_full," (",round(query_full/query_N*100),"%)"))
  dtbarplot$col="steelblue2"
  metadata(TREGELDataSet)$barplot_summary_dt <- dtbarplot

  temp=sapply(metadata(TREGELDataSet)$subjectNames,function(x){length(unique(metadata(TREGELDataSet)$detailDT[metadata(TREGELDataSet)$detailDT$subjType==x,][["queryID"]]))})
  temp=data.table(x=names(temp),Subjects=temp,lab=paste0("Queries with ",names(temp),": ",temp," (",round(temp/query_N*100),"%)"),col="steelblue3")
  subjPerQuery=round(sapply(metadata(TREGELDataSet)$subjectNames,function(x){sum(metadata(TREGELDataSet)$detailDT$subjType==x)})/query_N,digits=1)
  subjPerQuery=data.table(x=names(subjPerQuery),Subjects=subjPerQuery,lab=paste0("Average ",names(subjPerQuery)," per subject: ",subjPerQuery),col="steelblue4")

  dtbarplotDetailed=rbind(dtbarplot,temp,subjPerQuery, use.names=FALSE)
  dtbarplotDetailed$sequence=dim(dtbarplotDetailed)[1]:1
  dtbarplotDetailed$xtext=dtbarplotDetailed$sequence+0.55
  dtbarplotDetailed$query=log2(dtbarplotDetailed$query+1)#add +1 to avoid log2(1)=0 and log2(0)=-inf
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

#'Gviz plotting
#' @param TREGELDataSet A TREGELDataSet.
#' @return TREGELDataSet.
#' @export
gvizPlot <- function(TREGELDataSet,query,trackRegionLabels=setNames(rep("ID",length(TREGELDataSet)),names(TREGELDataSet))){
  for(q in query){
    #Subject tracks
    GvizSubjTracks=lapply(metadata(TREGELDataSet)$subjectNames,function(x){
      tmp=metadata(TREGELDataSet)$detailDT[subjType==x & queryID==q]
      if(dim(tmp)[1]>25){tmp=tmp[sample(1:dim(tmp)[1],size = 25),]}
      regionLabels=mcols(TREGELDataSet[[x]])[[trackRegionLabels[x]]][match(tmp$subjID,mcols(TREGELDataSet[[x]])[["ID"]])]
        if(dim(tmp)[1]!=0){# if subjectType overlaps with query create track
          AnnotationTrack(start=tmp$subjStart,end = tmp$subjEnd,chr=tmp$subjChr,name=x,id = regionLabels,
                        fontsize.title=24,featureAnnotation="id",fontcolor.title="black",fontcolor="black",
                        fontcolor.group="black",fontcolor.item="black",rotation.item=20)
        }else{
          AnnotationTrack(GRanges(),name = x,fontcolor.title="black",fontsize.title=24)
        }
    })
    #Gviz helper tracks
    chr=GvizSubjTracks[[1]]@chromosome
    itrack =IdeogramTrack(genome="hg38",chromosome=chr)#ideogram track (chromosome bands etc)
    gtrack=GenomeAxisTrack()
    queryGR=TREGELDataSet[[1]][mcols(TREGELDataSet[[1]])$ID==q]
    from = start(queryGR)-300  #add 300bp left and righ as plotWindow
    to = end(queryGR)+300
    plotWindow=GRanges(seqnames=chr,ranges=IRanges::IRanges(start = from, end = to))
    #Gviz query track
    regionLabels=mcols(queryGR)[[trackRegionLabels[1]]]
    queryTrack=AnnotationTrack(range=queryGR,name="Query",fill="red",arrowHeadWidth=30,shape="fixedArrow",featureAnnotation="id",
                               id=regionLabels,fontsize.group=20,fontsize.title=24,fontcolor.title="black")#,stacking = "dense"
    allTracks=c(itrack,gtrack,queryTrack,GvizSubjTracks)#;names(temp)=make.unique(names(temp))
     #Gviz plotting
    pdf(file.path(metadata(TREGELDataSet)$gvizPlotsFolder,paste0(q,".pdf")),width = 30/2.54,height = 20/2.54)
    plotTracks(allTracks,showOverplotting=TRUE,from = from, to = to,title.width=6,rotation.title=0,background.title="white",just.group="above") #20%BP up/down
    dev.off()
  return(TREGELDataSet)
  }


}


#' Add together two numbers

# extend <- function(x, upstream=0, downstream=0)     #https://support.bioconductor.org/p/78652/
# {
#   if (any(strand(x) == "*"))
#     warning("'*' ranges were treated as '+'")
#   on_plus <- strand(x) == "+" | strand(x) == "*"
#   new_start <- start(x) - ifelse(on_plus, upstream, downstream)
#   new_end <- end(x) + ifelse(on_plus, downstream, upstream)
#   ranges(x) <- IRanges(new_start, new_end)
#   trim(x)
# }
# subsetTFBS400=function(x,chrOI){
#   tmp=GRanges()
#   for(f in list.files(file.path(currDir,"data/anno/TFBShg"),full.names = TRUE)){
#     if(isFALSE(sub(".bed","",tail(unlist(strsplit(f,"/")),n=1))%in%chrOI)){next}
#     message(print(f))
#     TFBSonChr=readRDS(f)
#     tmp=c(tmp,TFBSonChr[TFBSonChr$TFBSID%in%x])
#   }
#   tmp$Name=tmp$TFBSID
#   tmp$TFBSID=make.unique(tmp$TFBSID)
#   return(tmp)
# }
# getTFBS=function(x){
#   TFBS=data.table()
#   for(i in 1:length(x)){ #for every gene take TFBS
#     chr=as.character(seqnames(x))[i]; if(chr=="MT"){chr="M"}
#     start=start(x)[i]
#     end=end(x)[i]
#     tmp=fread(cmd=paste0("./scripts/bigBedToBed http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/JASPAR2020_hg38.bb -chrom=chr",chr," -start=",start, " -end=",end, " stdout"),
#               col.names = c("chr","start","end","Name","score"),select=1:5)
#     if(length(tmp)==0)next #skip if there is no TFBS for gene (tmp is null)
#     tmp[,chr:=sub("chr","",tmp$chr,fixed = TRUE)];if(tmp[1,1]=="M"){tmp[,chr:="MT"]}
#     tmp[,Name:=toupper(tmp$Name)]
#     tmp[,TFBSID:=make.unique(tmp$Name)]
#     TFBS=rbindlist(list(TFBS,tmp))
#   }
#   TFBS=makeGRangesFromDataFrame(TFBS,keep.extra.columns=TRUE)
#   return(TFBS)
# }
