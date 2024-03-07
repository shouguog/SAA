rm(list = ls())
#################################
## variance partition analysis ##
#################################
library(variancePartition)
library(dplyr)
library(tidyverse)
# library(BiocParallel)
## in the manuscript, we provided all celltypes regardless of how terminally defined they were.
## in this example, we provide only the terminally defined populations as template code for the approach, the results are the same.
load("../count_meta.RData")
rownames(metaDataMerge)<-metaDataMerge$V1
colnames(metaDataMerge)<-c("cellID", "libID", "celltype")

subjects<-c(); disease<-c(); tmt<-c()
for(ii in 1:dim(metaDataMerge)[1]){
  items<-strsplit(metaDataMerge$libID[ii], "_")[[1]]
  subjects<-c(subjects, items[3])
  if(grepl("HD", items[6])){
    disease<-c(disease, "HD"); tmt<-c(tmt, "NO")
  }
  if(items[6]=="baseline"){
    disease<-c(disease, "SAA"); tmt<-c(tmt, "NO")
  }
  if(items[6]=="3M"||items[6]=="6M"){
    disease<-c(disease, "SAA"); tmt<-c(tmt, "YES")
  }
}
metaDataMerge$disease<-disease
metaDataMerge$subject<-subjects
metaDataMerge$tmt<-tmt
metaDataMerge$sevirity<-"VSAA"
#TAYC	SAA
#BELC	SAA
#QUET	SAA
#LOPA	SAA
#SIMT	SAA
#WELR	SAA
#AGUD	SAA
metaDataMerge$sevirity[grepl("TAYC", metaDataMerge$libID)]<-"SAA"
metaDataMerge$sevirity[grepl("BELC", metaDataMerge$libID)]<-"SAA"
metaDataMerge$sevirity[grepl("QUET", metaDataMerge$libID)]<-"SAA"
metaDataMerge$sevirity[grepl("LOPA", metaDataMerge$libID)]<-"SAA"
metaDataMerge$sevirity[grepl("SIMT", metaDataMerge$libID)]<-"SAA"
metaDataMerge$sevirity[grepl("WELR", metaDataMerge$libID)]<-"SAA"
metaDataMerge$sevirity[grepl("AGUD", metaDataMerge$libID)]<-"SAA"
###remove HD
#metaDataMerge<-metaDataMerge[metaDataMerge$disease=="SAA",]
geneExpression<-geneExpression[, rownames(metaDataMerge)]


form <- as.formula("~ (1|disease) + (1|sevirity) + (1|tmt)")
varPart_disease_sevirity_tmt <- fitExtractVarPartModel(geneExpression, form, metaDataMerge)


library(ggtern)
png("varPart_disease_sevirity_tmt.png", width = 1000, height = 1000, res = 200)
ggtern(data=varPart_disease_sevirity_tmt, aes(disease, sevirity, tmt)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  geom_point(size= 0.8, alpha = 0.5, color="black")  +
  geom_density_tern(
    bins=10,
    bdl = 0.02,
    color="black",
    base = "identity", n = 20,
    bdl.val = 0.1) +
  Llab("disease") +
  Tlab("sevirity") +
  Rlab("tmt") 
dev.off()


