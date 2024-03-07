#R4.0.4 on PC
##R version 4.1.3 (2022-03-10) -- "One Push-Up"
setwd("/data/gaos2/2022H0150SAAProject/AnalysisNew/BMRPCA/")
rm(list = ls())
library(Seurat)
load("seuratObjList_BM_step5.RData")
metaData<-sAA.combined@meta.data
#sAA.combined<-SetIdent(sAA.combined, "integrated_snn_res.0.8")
#Idents(sAA.combined)<-"integrated_snn_res.0.5"
#for(cluster in unique(metaData$integrated_snn_res.0.5)){
#  markers <- FindMarkers(sAA.combined, ident.1 = cluster,print.bar = TRUE, only.pos = TRUE, logfc.threshold = 0.2)
#  write.csv(markers, file=paste("BM.0.5_cluster_", cluster, ".csv",sep=""))
#}

#Idents(sAA.combined)<-"integrated_snn_res.0.8"
#for(cluster in unique(metaData$integrated_snn_res.0.8)){
#  markers <- FindMarkers(sAA.combined, ident.1 = cluster,print.bar = TRUE, only.pos = TRUE, logfc.threshold = 0.2)
#  write.csv(markers, file=paste("BM.0.8_cluster_", cluster, ".csv",sep=""))
#}

Idents(sAA.combined)<-"integrated_snn_res.1.5"
for(cluster in unique(metaData$integrated_snn_res.1.5)){
  markers <- FindMarkers(sAA.combined, ident.1 = cluster,print.bar = TRUE, only.pos = TRUE, logfc.threshold = 0.2)
  write.csv(markers, file=paste("BM.1.5_cluster_", cluster, ".csv",sep=""))
}



