#' loadAnnotations
#'
#' @param TREGELDataSet A TREGELDataSet.
#' @return TREGELDataSet.

loadAnnotations <- function(TREGELDataSet){
  TREGELDataSet <- readQuery(TREGELDataSet)
  TREGELDataSet <- readSubject(TREGELDataSet)

 return(TREGELDataSet)
}




#' Read query
readQuery=function(TREGELDataSet){
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
fOverlaps <- function(TREGELDataSet){
  detailDT <- data.table() #data table to store all overlaps for query vs all subjects
  sumDT <- data.table()
  s <- NULL
  for(subj in names(TREGELDataSet[2:length(TREGELDataSet)])){
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
