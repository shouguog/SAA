library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(tidyverse)

load("../seuratObjList_BM_step5.RData")
load("../plotResultAliceAnnotation/metaData_AliceAnnotation_CD4CD8.RData")
sAA.combined@meta.data<-metaData[rownames(sAA.combined@meta.data),]
metaData<-metaData[grepl("baseline|_HD", metaData$sampleID),]
sAA.combined<-subset(sAA.combined, cells = rownames(metaData))

##--- SCE conversion
sce <- as.SingleCellExperiment(sAA.combined)
altExps(sce) <- NULL                           

reducedDim(sce, "PCA", withDimnames=TRUE) <- sAA.combined@reductions$pca@cell.embeddings   
reducedDim(sce, "UMAP", withDimnames=TRUE) <- sAA.combined@reductions$umap@cell.embeddings
save(sce, file = "sce.RData")
