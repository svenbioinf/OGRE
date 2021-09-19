## app.R ##
#RegulationFinder, RegElements, Tr eg el TREGEL Transcription Regulatory Elements ReTrEl
#-order detailed table by summary table
#-3D data include! or H3K4me3 marks
#check for global variables
#1000bp enhancer?
#check what isoforms are used for TP53 human
#check fshb cgi
#highlight start button
#user defined tfbs
#ucsc genome browser for human
#take only one ccds ?
#todo:
# - TFBS is big and a lot of overlapping ranges, optimize?
library(shiny)
library(Gviz)
library(DT)
library(data.table)
library(biomaRt)
library(GenomicRanges)
library(shinydashboard)
library(shinycssloaders)
library(ggplot2)
library(reshape2)
library(stringr)
library(shinyBS)

#load functions
source("methods.R")
source("classes.R")

#for testing
if(FALSE){
  geneSymbols=c("ACSL4","EBF3","FAM57B","CXCR4")
  geneSymbols=read.table(file=file.path("/mnt/idata/LinaSeq/new/output/DE","DEG_padj0.01_FC2_padjSorted_TPM.csv"))$gene #read DEG of Lina's experiment
  geneSymbols=toupper(geneSymbols)
}

TREGELtest <- new("TREGELDataSet",queryFolder="/mnt/idata/TREGEL/tool/TREGEL3/data/human/query")

TREGELtest <- loadAnnotations(TREGELtest)

#preload annotations
currDir=getwd()
annoDir=file.path(currDir,"data/anno")

mmgeneCoord=readRDS(file.path(currDir,"data/genes/biomartGenesOkt2020mm38.rds"))
TFBSexpressed=readRDS(file.path(currDir,"data/TFBS","TFBSexpressed.RDS"))
mmitracks=readRDS(file.path(currDir,"data/gviz/mmitracks.RDS"))
mmregEle=list()
mmregEle[["Promoter"]]=readRDS(file.path(annoDir,"mmPromoters.rds"))
mmregEle[["CGI"]]=readRDS(file.path(annoDir,"mmCGI.rds"))
mmregEle[["TFBS"]]=readRDS(file.path(annoDir,"mmTFBS.rds"))
mmExon=readRDS(file.path(annoDir,"exonMM.rds"))

hggeneCoord=readRDS(file.path(currDir,"data/genes/biomartGenesOkt2020hg38.rds"))
TFBSexpressed=readRDS(file.path(currDir,"data/TFBS","TFBSexpressed.RDS"))
hgitracks=readRDS(file.path(currDir,"data/gviz/hgitracks.RDS"))
hgregEle=list()
hgregEle[["Promoter"]]=readRDS(file.path(annoDir,"hgPromoters.rds"))
hgregEle[["CGI"]]=readRDS(file.path(annoDir,"hgCGI.rds"))
hgregEle[["TFBS"]]=readRDS(file.path(annoDir,"hgTFBS.rds"))
hgExon=readRDS(file.path(annoDir,"exonHG.rds"))








  GvizGeneID=NULL


  #file upload
  #"chr,start,end,strand,ID" sep=",", NO header
  customData=list()
  #CustomDataLabels=make.unique(rep("CustomData",dim(input$files)[1]))
    i=0
    for(inputf in input$files$datapath){
      i=i+1
      temp=fread(inputf,select = c(1,2,3,4,5),col.names = c("Chr","Start","End","Strand","ID"))
      temp=makeGRangesFromDataFrame(temp,keep.extra.columns=TRUE)
      customData[[CustomDataLabels[i]]]=temp
    }




  #Main Start Script
    #' Main function.
    #'
    #' @param x A number.
    #' @param y A number.
    #' @return The sum of \code{x} and \code{y}.
    #' @examples
    #' add(1, 1)
    #' add(10, 1)
    findRegEle(Genome="human"){
    setProgress(message = "Analysis started\nthis might take a moment...");Sys.sleep(0.5)

    #Loading annotations
      loadAnnotations()
    }



    #Custom Data & custom selection
    if(input$CustomData==TRUE){v$regEle=c(v$regEle,customDataR())}
    if(input$checkboxPromoter!=TRUE){v$regEle["Promoter"]=NULL}
    if(input$checkboxCGI!=TRUE){v$regEle["CGI"]=NULL}
    if(input$checkboxTFBS==FALSE&input$checkboxTFBSexpressed==FALSE){v$regEle[["TFBS"]]=NULL}
    #get gene symbols
    geneSymbols=strsplit(input$genes,split = "\n")[[1]]
    geneSymbols=toupper(geneSymbols) #change to upper case names for later compatibility
    #get coordinates for DE genes
    if(TRUE){
      setProgress(message = "Accessing Ensembl Biomart");Sys.sleep(0.5)
      v$geneCoord=v$geneCoord[v$geneCoord$symbol%in%geneSymbols,]
      symbolsNotFound=geneSymbols[!(geneSymbols%in%v$geneCoord$symbol)]
      if(is.null(input$regionFile)){updateTextInput(session, "symbolsNotFound", value = paste(symbolsNotFound,collapse ="\n"))}
      #v$geneCoord=v$geneCoord[v$geneCoord$ccds!="",]#Use  ccds consensus transcript https://www.ensembl.org/info/genome/genebuild/ccds.html, Sometimes there is multiple CCDS for one gene symbol, it seems like they have the same start/end coords though. So I take the first one
      v$geneCoord=v$geneCoord[!duplicated(v$geneCoord$symbol),]#take only one ccds
      v$geneCoord$Chr=v$geneCoord$chromosome_name
      colnames(v$geneCoord)=c("chromosome_name", "start_position", "end_position", "strand",    "Gene",    "GeneID","Chr")
      v$geneCoord=v$geneCoord[,c(1:5,7,6)]
      v$geneCoord=makeGRangesFromDataFrame(v$geneCoord,seqnames.field = "chromosome_name",start.field = "start_position",end.field = "end_position",strand.field="strand",keep.extra.columns = TRUE)
    }
    #use regions file if given
    if(!is.null(input$regionFile)){
      temp=fread(input$regionFile$datapath,select = c(1,2,3,4,5),col.names = c("Chr","Start","End","Strand","ID"))
      temp=makeGRangesFromDataFrame(temp,keep.extra.columns=TRUE)
      temp$Gene=temp$ID;temp$Chr=as.character(seqnames(temp));temp$GeneID=temp$ID;temp$ID=NULL
      v$geneCoord=temp
    }
    v$geneCoord=extend(v$geneCoord,input$start,input$end) #add user supplied bp offset
    if(input$sliderPromoter!=0){v$geneCoord=promoters(v$geneCoord,upstream =input$sliderPromoter,downstream = 0);v$geneCoord$Gene=paste0(v$geneCoord$Gene,"_",input$sliderPromoter,"bp.up")}
    #TFBSScore
    if(input$TFBSScore!=600){v$regEle[["TFBS"]]=getTFBS(v$geneCoord)}
    if(input$checkboxTFBSexpressed==TRUE){v$regEle["TFBS"]=v$regEle[["TFBS"]][v$regEle[["TFBS"]]$Name%in%toupper(strsplit(input$TFBSexpressed,split = "\n")[[1]])]}

    #analysis Find overlaps
    setProgress(message = "Searching for Regulatory Elements")
    summt=v$geneCoord #summary table
    summtD= data.table(GeneID=NULL,RegEleID=NULL,Type=NULL) #detailed table
    for(regEleN in names(v$regEle)){
      sbjct=v$regEle[[regEleN]]
      ol=findOverlaps(summt,sbjct); ol=data.table(q=queryHits(ol),s=mcols(sbjct)[subjectHits(ol),1])
      summtD=rbind(summtD,data.table(Gene=v$geneCoord$Gene[ol$q],
                                      GeneID=v$geneCoord$GeneID[ol$q],
                                      RegEleID=ol$s,
                                      Type=regEleN,
                                      Chr=as.character(seqnames(sbjct))[match(ol$s,mcols(sbjct)[,1])],
                                      Start=start(sbjct)[match(ol$s,mcols(sbjct)[,1])],
                                      End=end(sbjct)[match(ol$s,mcols(sbjct)[,1])]
                                      ))
      ol=ol[,list(s=paste(s,collapse = " ")),by="q"]
      mcols(summt)[ol$q,regEleN]=ol$s
      ol=summt
    }
    summt=mcols(ol)
    summt=as.data.table(mcols(ol))
    summtD=summtD[order(summtD$Gene),]
    summtD=na.omit(summtD)
    #Add additional information (Chr, Positions, etc)
    summt$Chr=v$geneCoord$Chr[match(summt$GeneID,v$geneCoord$GeneID)]
    summt$GeneStart=start(v$geneCoord)[match(summt$GeneID,v$geneCoord$GeneID)]
    summt$GeneEnd=end(v$geneCoord)[match(summt$GeneID,v$geneCoord$GeneID)]
    v$summt=summt
    summtC=summt[complete.cases(summt),]#Genes with all regEle
    summtInC=summt[!(complete.cases(summt)),]#Genes with some regEle
    for(regEleN in names(v$regEle)){
      i=which(str_count(summtC[[regEleN]]," ")>5)#reduce elements to fit in table column
      i2=gregexpr(" ",summtC[[regEleN]][i])#locate " "
      cut_Pos=unlist(lapply(i2, `[[`, 5))#locate 9th " "
      cutted=substr(summtC[[regEleN]][i],1,cut_Pos)#cutting
      summtC[i,regEleN]=paste0(cutted,"...")
    }



      #summary barplot
      if(is.null(input$regionFile)){genes_t=length(unique(geneSymbols))}else{genes_t=length(v$geneCoord)}
      genes_identified=length(unique(mcols(v$geneCoord)$GeneID))
      genes_w_minOne_ele=length(unique(summtD$GeneID))
      genes_w_all_ele=dim(summtC)[1]
      dtbarplot=data.table(x=c("genes_t","genes_identified","genes_w_minOne_ele","genes_w_all_ele"),Genes=c(genes_t,genes_identified,genes_w_minOne_ele,genes_w_all_ele))
      dtbarplot$x=factor(x = dtbarplot$x,levels = c("genes_w_all_ele","genes_w_minOne_ele","genes_identified","genes_t"))
      dtbarplot$lab=c(paste0("Total number of submitted genes: ",genes_t," (100%)"),paste0("Genes found in database: ",genes_identified," (",round(genes_identified/genes_t*100),"%)"),paste0("Genes with >=1 regulatory elements: ",genes_w_minOne_ele," (",round(genes_w_minOne_ele/genes_t*100),"%)"),paste0("Genes with all regulatory elements: ",genes_w_all_ele," (",round(genes_w_all_ele/genes_t*100),"%)"))
      dtbarplot$sequence=4:1
      dtbarplot$xtext=dtbarplot$sequence+0.55
      output$summtBarPlot <- renderPlot({
        ggplot(dtbarplot, aes(x = sequence, y =Genes)) +
          geom_bar(stat = "identity", width = 0.3,color="lightblue",fill="lightblue") +
          coord_flip() +labs(x = "", y = "") +theme_bw() +  theme_classic() +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
          geom_text(aes(x = xtext, y = 0, label = lab),size=6,hjust = 0, vjust = 1) +guides(size = FALSE)
      })
      #Summary barplot (detailed)
      temp=sapply(unique(summtD$Type),function(x){length(unique(summtD[summtD$Type==x,][["GeneID"]]))})
      temp=data.table(x=names(temp),Genes=temp,lab=paste0("Genes with ",names(temp),": ",temp," (",round(temp/genes_t*100),"%)"),sequence="NA",xtext="NA")
      regElePerRegion=round(sapply(unique(summtD$Type),function(x){sum(summtD$Type==x)})/temp$Genes,digits=1)
      regElePerRegion=data.table(x=names(regElePerRegion),Genes=regElePerRegion,lab=paste0("Average ",names(regElePerRegion)," per gene: ",regElePerRegion),sequence="NA",xtext="NA")

      dtbarplotDetailed=rbind(dtbarplot,temp,regElePerRegion)
      dtbarplotDetailed$sequence=dim(dtbarplotDetailed)[1]:1
      dtbarplotDetailed$xtext=dtbarplotDetailed$sequence+0.55
      output$summtBarPlotDetailed <- renderPlot(height=dim(dtbarplotDetailed)[1]*50+50,{
        ggplot(dtbarplotDetailed, aes(x = sequence, y =Genes)) +
          geom_bar(stat = "identity", width = 0.3,color="lightblue",fill="lightblue") +
          coord_flip() +labs(x = "", y = "") +theme_bw() +  theme_classic() +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
          geom_text(aes(x = xtext, y = 0, label = lab),size=7,hjust = 0, vjust = 1) +guides(size = FALSE)
      })

      #summary table
      output$BestHits =DT::renderDT(server=TRUE,{
          datatable(summtC,extensions = 'Buttons',callback = callback2,
            options = list(
            dom = 'B<"dwnldb">frtip',
            autoWidth=FALSE,
            scrollX = TRUE,
            pageLength=5,
            buttons = list(list(extend='csv',filename="BestHits"),list(extend='excel',filename="BestHits")))
          )
        })
      output$downloadb <- downloadHandler(
        filename = function() {
          paste("BestHitsFull.csv")
        },
        content = function(file) {
          fwrite(summtC, file)
        })

      #detailed table
      output$summtD = DT::renderDT(server=TRUE,{
        datatable(summtD,extensions = 'Buttons',selection=list(mode = 'multiple', selected = c(1)),callback = callback,
                  options = list(
                    dom = 'B<"dwnld">frtip',
                    autoWidth=FALSE,
                    scrollX = TRUE,
                    pageLength=5,
                    buttons = list(list(extend='csv',filename="DetailedTable"),list(extend='excel',filename="DetailedTable")))
        )
      })
     #download button
      output$download1 <- downloadHandler(
        filename = function() {
          paste("DetailedTableFull.csv")
        },
        content = function(file) {
          fwrite(summtD, file)
        }
      )

    #get CustomDataDisplay for gviz
    if(!is.null(input$filesDisplay)){v$CustomDataDisplay=getCustomDataDisplay()}
    #make global avail
    v$summtC=summtC
    v$summtD=summtD
    setProgress(message = "Done");Sys.sleep(0.5)



  #get data ready for gTable
  observe({
    req(v$summt)
    x=v$summt[,!c("Chr","GeneID","GeneStart","GeneEnd")]
    to.replace=names(x)[!(names(x)%in%c("Gene"))]
    for (var in to.replace) x[[var]]= as.integer(!is.na(x[,..var]))
    x[x=="1"]=as.character(icon("check",lib = "glyphicon"))
    x[x=="0"]=as.character(icon("minus",lib = "glyphicon"))
    v$summtG=x
    #
  })
  #show gene in genome browser
  output$my_test = renderUI({
    req(v$GvizGeneID)
    v$GvizGeneID
    geneGR=v$geneCoord[v$geneCoord$GeneID==v$GvizGeneID]
    geneStart=start(geneGR)
    geneEnd=end(geneGR)
    geneChr=as.character(seqnames(geneGR))
    if(input$Genome=="mouse"){
      urlx=  paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr",
                    geneChr,"%3A",
                    geneStart,"-",
                    geneEnd,#"%2D",
                    "&hgsid=898780277_NxNuEi2rA7EmfrSkCF9c9ujVsSa0"
      )
    }else if(input$Genome=="human"){
      urlx=  paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr",
                    geneChr,"%3A",
                    geneStart,"-",
                    geneEnd,#"%2D",
                    "&hgsid=898780277_NxNuEi2rA7EmfrSkCF9c9ujVsSa0"
      )
    }
    my_test=shiny::tags$iframe(src= urlx ,seamless = "seamless", style='width:100vw;height:100vh;')
    my_test
  })
  #Manage which selected Gene to plot
  observe({message("GeneSelecter")
      proxy_BestHits = dataTableProxy('BestHits')
      proxy_summtD = dataTableProxy('summtD')
  if(is.null(input$BestHits_rows_selected)&is.null(input$summtD_rows_selected)) {v$GvizGeneID=NULL}
      else if(!is.null(input$BestHits_rows_selected)&!is.null(input$summtD_rows_selected)){v$GvizGeneID=NULL}
      else if(is.null(input$BestHits_rows_selected)&!is.null(input$summtD_rows_selected)){
        selectRows(proxy = proxy_BestHits,selected = NULL)
        v$GvizGeneID=as.character(v$summtD[input$summtD_rows_selected[[1]],"GeneID"])
        v$chr <- as.character(v$summtD[v$summtD$GeneID==v$GvizGeneID,"Chr"][1])
      }
      else if(!is.null(input$BestHits_rows_selected)){
        selectRows(proxy = proxy_summtD,selected = NULL)
        v$GvizGeneID=as.character(v$summtC[input$BestHits_rows_selected[[1]],"GeneID"])
        v$chr <- as.character(v$summtC[v$summtC$GeneID==v$GvizGeneID,"Chr"][1])
      }
      else if(!is.null(input$summtD_rows_selected)){
        selectRows(proxy = proxy_BestHits,selected = NULL)
        v$GvizGeneID=as.character(v$summtD[input$summtD_rows_selected[[1]],"GeneID"])
        v$chr <- as.character(v$summtD[v$summtD$GeneID==v$GvizGeneID,"Chr"][1])
      }
  })

  calcGvizPlotData=observe({
    req(v$GvizGeneID,v$chr)
    message("gvizPlotdata")
    print(system.time({
    #setProgress(message = "Print");Sys.sleep(0.5)
    summtD=v$summtD[v$summtD$GeneID==v$GvizGeneID,]
    #Regele
    maxRegEleToDisplay=100
    GvizRegEleTracks=lapply(unique(summtD$Type),function(x){
      tmp=summtD[summtD$Type==x&summtD$GeneID==v$GvizGeneID,]
      if(dim(tmp)[1]>maxRegEleToDisplay){tmp=tmp[sample(1:dim(tmp)[1],size = maxRegEleToDisplay),]}
      if(input$checkboxPlotDisplayPredefiendTrackIDs==TRUE){
        AnnotationTrack(start=tmp$Start,end = tmp$End,chr=tmp$Chr,name=x,id = as.character(tmp$RegEleID),fontsize.title=24,featureAnnotation="id",fontcolor.title="black",fontcolor="black",fontcolor.group="black",fontcolor.item="black",rotation.item=20)
      }else{
        AnnotationTrack(start=tmp$Start,end = tmp$End,chr=tmp$Chr,name=x,stacking="dense",id = as.character(tmp$RegEleID),fontsize.title=24,fontcolor.title="black",fontcolor="black",rotation.item=20)
      }
    })
    #Additional Tracks to display
    displayTracks=list()
    returnList=list()
    itrack =v$itracks[[v$chr]]#ideogram track (chromosome bands etc)
    gtrack=GenomeAxisTrack()
    geneGR=v$geneCoord[v$geneCoord$GeneID==v$GvizGeneID]

    v$from = start(geneGR)-round(width(geneGR)*0.2)+1; v$to = end(geneGR)+round(width(geneGR)*0.2)+1 #20%BP up/down
    plotWindow=GRanges(seqnames=v$chr,ranges=IRanges(start = v$from, end = v$to))
    exon=subsetByOverlaps(v$exon, plotWindow);exon=exon[exon$ensembl_gene_id==v$GvizGeneID];mcols(exon)=NULL
    geneWexon=c(range(geneGR),exon);geneWexon$group=c(as.character(summtD$Gene[1]),rep("Exon",times=length(geneWexon)-1))
    #geneTrack=AnnotationTrack(range=geneGR,name="Gene",fill="red",stacking = "dense",arrowHeadWidth=30,shape="fixedArrow",id = as.character(summtD$Gene[1]),featureAnnotation="id",fontsize=20,fontsize.title=24,fontcolor.title="black")
    geneTrack=AnnotationTrack(range=geneWexon,name="Gene",fill="red",arrowHeadWidth=30,shape="fixedArrow",groupAnnotation = "group",fontsize.group=20,fontsize.title=24,fontcolor.title="black")#,stacking = "dense"

    if(input$Genome=="mouse"){
      if(input$checkboxGeneDisplay){displayTracks[["Genes.MM"]]=subsetByOverlaps(v$geneGR, plotWindow)}
      if(input$checkboxPromoterDisplay){displayTracks[["Promoter.MM"]]=subsetByOverlaps(hgregEle[["Promoter"]],plotWindow)}
      if(input$checkboxCGIDisplay){displayTracks[["CGI.MM"]]=subsetByOverlaps(hgregEle[["CGI"]],plotWindow)}
      if(isTRUE(input$checkboxTFBSDisplay)){displayTracks[["TFBS.MM"]]=subsetByOverlaps(hgregEle[["TFBS"]],plotWindow)}
      if(input$checkboxTFBSexpressedDisplay){displayTracks[["TFBSexpressed.MM"]]=subsetByOverlaps(hgregEle[["TFBS"]],plotWindow)}
      CustomDataDisplay=lapply(v$CustomDataDisplay,function(x){subsetByOverlaps(x,plotWindow)})
      displayTracks=c(displayTracks,CustomDataDisplay)
      for(i in names(displayTracks) ){
        name=i;x=displayTracks[[i]]
        if(input$checkboxPlotDisplayCustomTrackIDs==TRUE){
          returnList=c(returnList,AnnotationTrack(start=start(x),end = end(x),chr=as.character(seqnames(x)),strand = as.character(strand(x)),name=name,id=mcols(x)[,1],featureAnnotation="id",stacking="dense",fontsize.title=24,fontcolor.title="black",fontcolor="black",fontcolor.group="black",fontcolor.item="black"))
        }else{
          returnList=c(returnList,AnnotationTrack(start=start(x),end = end(x),chr=as.character(seqnames(x)),strand = as.character(strand(x)),name=name,id=mcols(x)[,1],stacking="dense",fontsize.title=24,fontcolor.title="black",fontcolor="black"))
        }
      }
    }
    else if(input$Genome=="human"){
      if(input$checkboxGeneDisplay){displayTracks[["Genes.HG"]]=subsetByOverlaps(v$geneGR, plotWindow)}
      if(input$checkboxPromoterDisplay){displayTracks[["PromoterHG"]]=subsetByOverlaps(hgregEle[["Promoter"]],plotWindow)}
      if(input$checkboxCGIDisplay){displayTracks[["CGI.HG"]]=subsetByOverlaps(hgregEle[["CGI"]],plotWindow)}
      if(isTRUE(input$checkboxTFBSDisplay)){displayTracks[["TFBS.HG"]]=subsetByOverlaps(hgregEle[["TFBS"]],plotWindow)}
      if(input$checkboxTFBSexpressedDisplay){displayTracks[["TFBSexpressed.HG"]]=subsetByOverlaps(hgregEle[["TFBS"]],plotWindow)}
      CustomDataDisplay=lapply(v$CustomDataDisplay,function(x){subsetByOverlaps(x,plotWindow)})
      displayTracks=c(displayTracks,CustomDataDisplay)
      for(i in names(displayTracks)){ #list to store named annotationTracks
        name=i;x=displayTracks[[i]]
        if(length(x)==0)next# do not put track in list for plotting if its empty
        if(input$checkboxPlotDisplayCustomTrackIDs==TRUE){
          if(name=="Genes.GH"){returnList=c(returnList,AnnotationTrack(start=start(x),end = end(x),chr=as.character(seqnames(x)),strand = as.character(strand(x)),name=name,id=mcols(x)[,1],featureAnnotation="id",fontsize.title=24,fontcolor.title="black",fontcolor="black",fontcolor.group="black",fontcolor.item="black"))}
          else{returnList=c(returnList,AnnotationTrack(start=start(x),end = end(x),chr=as.character(seqnames(x)),strand = as.character(strand(x)),name=name,id=mcols(x)[,1],featureAnnotation="id",fontsize.title=24,fontcolor.title="black",fontcolor="black",fontcolor.group="black",fontcolor.item="black"))}
        }else{
          if(name=="Genes.GH"){returnList=c(returnList,AnnotationTrack(start=start(x),end = end(x),chr=as.character(seqnames(x)),strand = as.character(strand(x)),name=name,id=mcols(x)[,1],fontsize.title=24,fontcolor.title="black",fontcolor="black"))}#
          else{returnList=c(returnList,AnnotationTrack(start=start(x),end = end(x),chr=as.character(seqnames(x)),strand = as.character(strand(x)),name=name,id=mcols(x)[,1],fontsize.title=24,fontcolor.title="black",fontcolor="black"))}
        }
      }
    }

    dummy=GRanges(Rle(c(geneTrack@chromosome), c(1)),IRanges(1, width=10:1))
    dummy=AnnotationTrack(start=start(dummy),end = end(dummy),chr=as.character(seqnames(dummy)),strand = as.character(strand(dummy)),name="123456789101112",fontsize.title=24)
    v$allTracks=c(itrack,gtrack,geneTrack,GvizRegEleTracks,returnList,dummy)#;names(temp)=make.unique(names(temp))
    #v$allTracks=c(v$allTracks,dummy)
    v$allTracksLength=100+80*(2+length(v$allTracks))
    }))
  })

  #Gviz gene plot
  output$plot1= renderPlot({
    req(v$allTracks)
    message("plotting")
    print(system.time({p=plotTracks(v$allTracks,showOverplotting=TRUE,from = v$from, to = v$to,title.width=6,rotation.title=0,background.title="white",just.group="above")})) #20%BP up/down
    p
  })#,height =function(){v$allTracksLength}
  output$Gtable = DT::renderDT({
   datatable(v$summtG,extensions = 'Buttons',escape = FALSE,
              options = list(
                dom = 'Bfrtip',
                autoWidth=FALSE,
                scrollX = TRUE,
                scrollY="40vh",
                paging = FALSE,
                #columnDefs = list(list(width = '50', targets = c(1,2,3))),
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
              ) #options = list(scrollY="300px",scrollX="300px", pageLength = 100, autoWidth = TRUE))
    )
  },server=FALSE)

#server end

if(FALSE){
  save(summtD, geneCoord, file = "tempData.RData")

  library(Gviz)
  maxRegEleToDisplay=100
  GvizGene="ACSL4"
  chr <- summtD$Chr[summtD$Gene==GvizGene][1]
  itrack <- IdeogramTrack(genome = "mm10", chromosome = chr)
  gtrack=GenomeAxisTrack()
  #geneTrack=AnnotationTrack(start=summt$geneStart[1],end = summt$geneEnd[1],chromosome = chr,name="Gene",col="red",fill="red")
  geneGR=geneCoord[geneCoord$Gene==GvizGene];colnames(mcols(geneGR))=c("group","Chr","GeneID")
  #geneTrack=GeneRegionTrack(range=geneGR)
  #geneTrack=GeneRegionTrack(range=geneGR,transcriptAnnotation = "Gene")
  #geneTrack=GeneRegionTrack(range=geneGR,groupAnnotation = "Gene")
  #if(as.character(strand(geneGR))[1]=="-"){offset5=6000;offset3=0}else{offset5=0;offset3=6000}
  #s=seq(start(geneGR)+3000,end(geneGR)-3000,3000);tempGR=GRanges(seqnames = seqnames(geneGR),IRanges(s,s+1000),strand = strand(geneGR))
  #mcols(tempGR)=mcols(geneGR);geneGR=c(geneGR,tempGR);geneGR$group=1:length(geneGR)
  geneTrack=AnnotationTrack(range=geneGR,name="Gene",fill="red",stacking = "dense",arrowHeadWidth=30,shape="fixedArrow",id = GvizGene,featureAnnotation="id")
  GvizRegEleTracks=lapply(names(regEle),function(x){
    tmp=summtD[summtD$Type==x&summtD$Gene==GvizGene,]
    if(dim(tmp)[1]>maxRegEleToDisplay){tmp=tmp[sample(1:dim(tmp)[1],size = maxRegEleToDisplay),]}
    AnnotationTrack(start=tmp$Start,end = tmp$End,chr=tmp$Chr,name=x,stacking="dense")

  })

  promoterTrack=AnnotationTrack(start=summtD$Start[1],end = summtD$End[1],chromosome = chr,name="Promoter")
  cpgTrack=AnnotationTrack(start=summtD$Start[6],end = summtD$End[1],chromosome = chr,name="CpG",fill="green")
  #TFBSGRlist=TFBSGR[queryHits(findOverlaps(TFBSGR,GRanges(seqnames = chr,ranges = IRanges(start = summt$promoterStart[1],end = summt$promoterEnd[1]))))];colnames(mcols(TFBSGRlist))=c("group")
  #TFBStrack=AnnotationTrack(range=TFBSGRlist,name = "TFBS",showId=TRUE)
  #plotTracks(list(itrack, gtrack, geneTrack,promoterTrack,cpgTrack),from = summt$promoterStart[1]-1000 ,to = summt$promoterEnd[1]+1000)
  plotTracks(list(itrack, gtrack, geneTrack,promoterTrack,cpgTrack))
  plotTracks(c(itrack,gtrack,geneTrack,GvizRegEleTracks),showOverplotting=TRUE)
  availableDisplayPars(geneTrack)
}



#initialization







# junk
if(FALSE){
  summt=v$summt[,!c("Chr","GeneID","GeneStart","GeneEnd")]

  to.replace=names(summt)[!(names(summt)%in%c("Gene"))]
  for (var in to.replace) summt[[var]]= as.integer(!is.na(summt[,..var]))
  #https://stackoverflow.com/questions/44762015/header-direction-in-shiny-data-table
  # https://stackoverflow.com/questions/30671958/how-to-embed-an-image-in-a-cell-a-table-using-dt-r-and-shiny
  #https://community.rstudio.com/t/icons-show-up-as-text-in-datatable/25495
  #https://gist.github.com/cecilialee/e848068cfd5862fd63961b03b42ead78
  #https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/


  if(!exists("urlx")){ #if url not set
    urlx=  paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr",
                  "1%3A","1%2D","10000&hgsid=898780277_NxNuEi2rA7EmfrSkCF9c9ujVsSa0"
    )
  }else{
    urlx=  paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr",
                  summt[input$BestHits_rows_selected[[1]],"geneChr"],"%3A",
                  summt[input$BestHits_rows_selected[[1]],"geneStart"],"%2D",
                  summt[input$BestHits_rows_selected[[1]],"geneEnd"],
                  "&hgsid=898780277_NxNuEi2rA7EmfrSkCF9c9ujVsSa0"
    )
  }
}

#1000bp distance enhancer https://www.researchgate.net/publication/49637582_Integrative_analysis_of_genomic_functional_and_protein_interaction_data_predicts_long-range_enhancer-target_gene_interactions
shinyApp(ui, server)
#library(rsconnect)
#setRepositories(ind=c(1,2,3,4))
#deployApp(appDir = "/mnt/idata/PCTReg/scripts/git/PCTReg/test")

# server commands
# scp -r /mnt/idata/PCTReg/scripts/git/PCTReg/test root@cera-p4:/home/shiny/TREGEL1.0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        # scp -r /mnt/idata/TREGEL/tool/TREGEL2 root@cera-p4:/home/shiny/TREGEL1.0
#(without sudo)
#danach noch test.sh auf neue version umstellen und eventuell packages installieren
#in /home/shiny/TREGEL1.0/packages.R
#danach shiny server neu starten und warten
#systemctl restart shiny.service

# if(FALSE){
#  list_obj_sizes = function(list_obj=ls(envir=.GlobalEnv)){
#   sizes = sapply(list_obj, function(n) object.size(get(n)), simplify = FALSE)
#   print(sapply(sizes[order(-as.integer(sizes))], function(s) format(s, unit = 'auto'))) }
# list_obj_sizes()
#
# list_obj_sizes <- function(list_obj=ls(envir=.GlobalEnv)){
#   sizes <- sapply(list_obj, function(n) object.size(get(n)), simplify = FALSE)
#   print(sapply(sizes[order(-as.integer(sizes))], function(s) format(s, unit = 'auto'))) }
#
# list_obj_sizes <- function(list_obj=names(v)){
#   sizes <- sapply(list_obj, function(n) object.size(n), simplify = FALSE)
#   print(sapply(sizes[order(-as.integer(sizes))], function(s) format(s, unit = 'auto'))) }
# list_obj_sizes()
# }
#

# system.time({plotTracks(v$allTracks,showOverplotting=TRUE,from = v$from, to = v$to,title.width=6,rot.title=0,background.title="white")})
#
# saveRDS(v,"v.RDS")
# v=readRDS("v.RDS")


#tmp=fread(cmd="./bigBedToBed http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2020/JASPAR2020_hg38.bb -chrom=chr21 -start=25900000 -end=30000000 stdout",
#          col.names = c("chr","start","end","name","score"),select=1:5)
#tmp[,chr:=sub("chr","",tmp$chr,fixed = TRUE)]
#tmp[,name:=toupper(tmp$name)]

#package created with usethis::create_package("TREGEL3")
#check if package is valid with devtools::check()
#create documentation with devtools::document()

#
