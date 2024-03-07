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
nhoodCounts<-milo@nhoodCounts
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

####Let us check the hoodInformation
hoods<-colnames(nhoods)
nhoods_node_names<-c()
index<-0
for(hood in hoods){
  cat(index);cat("\n");index<-index+1;
  nhoods_node<-nhoods[, hood]
  nhoods_node<-nhoods_node[nhoods_node>0]
  nhoods_node_names<-c(nhoods_node_names, paste(names(nhoods_node), collapse = ":"))
}
names(nhoods_node_names)<-hoods
###Add da_results
da_results$Nhoodname<-as.character(unlist(nhoodIndex)[da_results$Nhood])
da_results$nhoods_node_names<-nhoods_node_names[da_results$Nhoodname]
da_results$countForCheck<-rowSums(nhoodCounts)
da_results$countCell<-unlist(lapply(strsplit(da_results$nhoods_node_names, ":"), length))
da_results<-da_results[order(da_results$PValue),]
save(da_results, file = "da_results.RData")
write.csv(da_results, file = "da_results.csv", quote = FALSE)
