library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(tidyverse)

load("sce.RData")
##--------------------------------------##
##    Differential abundance testing    ##
##--------------------------------------##

##--- Step1: Create a Milo object
milo <- Milo(sce)
rm(sce)
##--- Step2: Construct KNN graph
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA") 
proportion<-0.1
##--- Step3: Defining representative neighbourhoods on the KNN graph
milo <- makeNhoods(milo, prop = proportion, k = 30, d = 30, refined = TRUE, reduced_dims = "PCA")   
# check average neighbourhood size: over 5 x N_samples is required 
p.NhSize = plotNhoodSizeHist(milo)
ggsave(file = paste0("Distribution_Neighbourhood.size_ALL.png"), plot = p.NhSize, width = 10, height = 7, dpi = 100)

##--- Step4: Counting cells in neighbourhoods
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="sampleID")

##--- Step6: Computing neighbourhood connectivity---##
milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")
milo <- buildNhoodGraph(milo)

save(milo, file = "milo.RData")
