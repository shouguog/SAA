#R4.0.4 on PC
##R version 4.1.3 (2022-03-10) -- "One Push-Up"
setwd("/data/gaos2/2022H0150SAAProject/AnalysisNew/BMRPCA/")
rm(list = ls())
library(Seurat)
load("seuratObjList_BM_step1.RData")

sAA.combined <- ScaleData(sAA.combined, verbose = FALSE)
sAA.combined <- RunPCA(sAA.combined, verbose = FALSE)
sAA.combined <- RunUMAP(sAA.combined, dims = 1:50)
#find neighbor
sAA.combined <- FindNeighbors(sAA.combined, dims = 1:50)
sAA.combined <- FindClusters(sAA.combined, resolution = 0.8)
save(sAA.combined, file = "seuratObjList_BM_step2.RData")
sAA.combined <- FindClusters(sAA.combined, resolution = 0.5)
save(sAA.combined, file = "seuratObjList_BM_step3.RData")
sAA.combined <- FindClusters(sAA.combined, resolution = 1)
save(sAA.combined, file = "seuratObjList_BM_step4.RData")
sAA.combined <- FindClusters(sAA.combined, resolution = 1.5)
save(sAA.combined, file = "seuratObjList_BM_step5.RData")


#https://satijalab.org/seurat/articles/integration_large_datasets.html
#bm280k.integrated <- ScaleData(bm280k.integrated, verbose = FALSE)
#bm280k.integrated <- RunPCA(bm280k.integrated, verbose = FALSE)
#bm280k.integrated <- RunUMAP(bm280k.integrated, dims = 1:50)