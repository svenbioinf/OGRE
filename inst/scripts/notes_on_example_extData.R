# Example data consists of human protein coding genes,CGI, SNP and TFBS (hg19), 
# taken from AnnotationHub using OGRE and reduced to Chr21. The aquired GRanges
# data is then stored under
# the extdata folder as .RDS files. It can be easily reproduced using below 
# provided R code.
# 
# List of example datasets:
#
# -protCodingGenes - Protein coding genes from HG19 (GRCh37) Ensembl
#   For additional information use:
#   `getInfoOnIds(AnnotationHub(), "AH10684")` 
# -CGI - CpG islands from HG19 UCSC
#   For additional information use:
#   `getInfoOnIds(AnnotationHub(), "AH5086")`
# -SNP - Common Single Nucleotide Polymorphism from HG19 UCSC
#   For additional information use:
#   `getInfoOnIds(AnnotationHub(), "AH5105")`
# -TFBS - Transcription Factor Binding Sites conserved from HG19 UCSC
#   For additional information use:
#   `getInfoOnIds(AnnotationHub(), "AH5090")`

library(OGRE)
myOGRE <- OGREDataSet()
myOGRE <- addDataSetFromHub(myOGRE,"protCodingGenes","query")
myOGRE$protCodingGenes <- keepSeqlevels(myOGRE$protCodingGene,"21","coarse")
saveRDS(myOGRE$protCodingGenes,"/path/genes.RDS")
myOGRE <- addDataSetFromHub(myOGRE,"CGI","subject")
myOGRE$CGI <- keepSeqlevels(myOGRE$CGI,"21","coarse")
saveRDS(myOGRE$CGI,"/path/CGI.RDS")
myOGRE <- addDataSetFromHub(myOGRE,"TFBS","subject")
myOGRE$TFBS <- keepSeqlevels(myOGRE$TFBS,"21","coarse")
saveRDS(myOGRE$TFBS,"/path/TFBS.RDS")

