#R4.0.4 on PC
##R version 4.1.3 (2022-03-10) -- "One Push-Up"
setwd("/data/gaos2/2022H0150SAAProject/AnalysisNew/BMRPCA/")
rm(list = ls())
library(Seurat)
load("../seuratObjList_BM_step0.RData")

cellNum<-lapply(X = seuratObjList_BM, FUN = function(x) {
  x <- dim(x@meta.data)[1]
})


seuratObjList_BM <- lapply(X = seuratObjList_BM, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seuratObjList_BM)

seuratObjList_BM <- lapply(X = seuratObjList_BM, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
####Use the HD as reference 3:6
sAA.anchors <- FindIntegrationAnchors(object.list = seuratObjList_BM, reference = 3:6, reduction = "rpca",dims = 1:50)
# this command creates an 'integrated' data assay
sAA.combined <- IntegrateData(anchorset = sAA.anchors, dims = 1:50)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(sAA.combined) <- "integrated"
save(sAA.combined, file = "seuratObjList_BM_step1.RData")

sAA.combined <- ScaleData(sAA.combined, verbose = FALSE)
sAA.combined <- RunPCA(sAA.combined, verbose = FALSE)
sAA.combined <- RunUMAP(sAA.combined, dims = 1:50)
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