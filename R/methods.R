#' Add together two numbers
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
loadAnnotations <- function(obj){
 obj <- readQuery(obj)
 return(obj)
}

  if (FALSE){
    query(obj)=readRDS(list.files(queryFolder(obj),full.names =TRUE )[[1]])
  colnames(mcols(query(obj)))=c("name","id")

  v$geneGR=makeGRangesFromDataFrame(hggeneCoord,keep.extra.columns = TRUE)
    v$itracks=hgitracks
    v$regEle=hgregEle
    v$exon=hgExon

    v$geneCoord=mmgeneCoord
    colnames(mmgeneCoord)=c("chr","start","end","strand","symbol","id");v$geneGR=makeGRangesFromDataFrame(mmgeneCoord,keep.extra.columns = TRUE)
    v$itracks=mmitracks
    v$regEle=mmregEle
    v$exon=mmExon
  return(obj)
  }



#' Add together two numbers
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
readQuery=function(obj){
  if(is.na(queryFolder(obj))){ #no folder supplied, check default folder
    #query(obj) <- readRDS()
    #colnames(mcols(query(obj)))=c("name","id")
    #query(obj)<- makeGRangesFromDataFrame(query(obj),keep.extra.columns = TRUE)
    #obj <- readQuery(obj)

  }else{ #read queryFolder
    queryPath <- list.files(queryFolder(obj),full.names = TRUE)
    query(obj) <- readRDS(queryPath)
    query(obj)<- makeGRangesFromDataFrame(query(obj),keep.extra.columns = TRUE)
    colnames(mcols(query(obj)))=c("name","id")
  }
  return(obj)
}
#' Add together two numbers
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
readRubject=function(obj){

}

#' Add together two numbers
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
extend <- function(x, upstream=0, downstream=0)     #https://support.bioconductor.org/p/78652/
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}
subsetTFBS400=function(x,chrOI){
  tmp=GRanges()
  for(f in list.files(file.path(currDir,"data/anno/TFBShg"),full.names = TRUE)){
    if(isFALSE(sub(".bed","",tail(unlist(strsplit(f,"/")),n=1))%in%chrOI)){next}
    message(print(f))
    TFBSonChr=readRDS(f)
    tmp=c(tmp,TFBSonChr[TFBSonChr$TFBSID%in%x])
  }
  tmp$Name=tmp$TFBSID
  tmp$TFBSID=make.unique(tmp$TFBSID)
  return(tmp)
}
getTFBS=function(x){
  TFBS=data.table()
  for(i in 1:length(x)){ #for every gene take TFBS
    chr=as.character(seqnames(x))[i]; if(chr=="MT"){chr="M"}
    start=start(x)[i]
    end=end(x)[i]
    tmp=fread(cmd=paste0("./scripts/bigBedToBed http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/JASPAR2020_hg38.bb -chrom=chr",chr," -start=",start, " -end=",end, " stdout"),
              col.names = c("chr","start","end","Name","score"),select=1:5)
    if(length(tmp)==0)next #skip if there is no TFBS for gene (tmp is null)
    tmp[,chr:=sub("chr","",tmp$chr,fixed = TRUE)];if(tmp[1,1]=="M"){tmp[,chr:="MT"]}
    tmp[,Name:=toupper(tmp$Name)]
    tmp[,TFBSID:=make.unique(tmp$Name)]
    TFBS=rbindlist(list(TFBS,tmp))
  }
  TFBS=makeGRangesFromDataFrame(TFBS,keep.extra.columns=TRUE)
  return(TFBS)
}
