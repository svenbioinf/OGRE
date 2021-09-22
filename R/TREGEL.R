#' @author Sven Berres

#library(shiny)
#library(Gviz)
#library(DT)
#library(data.table)
#library(biomaRt)
#library(GenomicRanges)
#library(shinydashboard)
#library(shinycssloaders)
#library(ggplot2)
#library(reshape2)
#library(stringr)
#library(shinyBS)

#load functions



#TREGELtest <- TREGELDataSet(queryFolder="/mnt/idata/TREGEL/tool/TREGEL3/data/human/query")
#TREGELtest <- loadAnnotations(TREGELtest)
#TREGELtest <- loadAnnotations(TREGELtest)




#1000bp distance enhancer https://www.researchgate.net/publication/49637582_Integrative_analysis_of_genomic_functional_and_protein_interaction_data_predicts_long-range_enhancer-target_gene_interactions
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
#check if package is valid with
#create documentation with devtools::document()

#
