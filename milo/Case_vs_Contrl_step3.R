rm(list = ls())
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(tidyverse)
load("milo.RData")

head(nhoodCounts(milo))

nhoods<-milo@nhoods  ###This includes the information of cell assignment
dim(nhoods)
nhoodIndex<-milo@nhoodIndex
nhoodGraph<-milo@nhoodGraph
####Add status
milo@colData$status<-"PT"
milo@colData$status[grepl("_HD",milo@colData$sampleID)]<-"HD"
table(milo@colData$status)

#Now we are all set to test for differential abundance in neighbourhoods. We implement this hypothesis testing in a 
#generalized linear model (GLM) framework, specifically using the Negative Binomial GLM implementation in edgeR.
#We first need to think about our experimental design. The design matrix should match each sample to the experimental 
#condition of interest for DA testing. In this case, we want to detect DA between embryonic stages, stored in the stage 
#column of the dataset colData. We also include the sequencing.batch column in the design matrix. This represents a known 
#technical covariate that we want to account for in DA testing.
## Convert batch info from integer to factor
embryo_design <- data.frame(colData(milo))[,c("sampleID", "status")]
embryo_design <- distinct(embryo_design)
rownames(embryo_design) <- embryo_design$sampleID

da_results <- testNhoods(milo, design = ~ status, design.df = embryo_design, reduced.dim="PCA")
head(da_results)


## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo, dimred = "UMAP", colour_by="status", text_by = "celltype",text_size = 3, point_size=0.5) + guides(fill="none")
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.1) 
png("umap_1.png", width = 2000, height = 1000, res = 100)
gridExtra::grid.arrange(umap_pl,nh_graph_pl, ncol=2)
dev.off()

## Warning: ggrepel: 5 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
plotNhoodGroups(milo, da_results, layout="umap") 
set.seed(42)
da_results_1 <- groupNhoods(milo, da_results, max.lfc.delta = 5, overlap=5)
png("umap_2.png", width = 2000, height = 1000, res = 100)
plotNhoodGroups(milo, da_results_1, layout="UMAP")
dev.off()

