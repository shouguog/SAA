###R4.0.4
library(Seurat)
library(Matrix)
rm(list=ls())
###Let us read the sinature genes
cutoff=0.001;
cutoffcluster=0.01;

library(GSEABase)
#From NHGRI UAB1 
#geneSets <- getGmt("../../../BMSeurat8runNoNorm/HCAtypeassignment/mmc4_200genes_BM_oneCD34.gmt")
geneSets <- getGmt("../res.0.5/mmc4_200genes_BM_oneCD34.gmt")
allgenes<-c()
for(ii in 1:16){
  allgenes<-c(allgenes, geneSets[[ii]]@geneIds)
}
#####All genes in study###
load("../../seuratData_BM_step0.RData")
overGenes<-rownames(seuratData_BM[[1]])

zz <- file("signaureGenesOverlapOneCD34.txt", "w")  # open an output file connection
cat("signature\tcluster\toverlapNum\toverlapGene\tsignaturegene\tclustergene\tallsignaturegene\n", file=zz)
for(sig in 1:16){
	sigGene<-as.character(geneSets[[sig]]@geneIds)
	sigGeneTop<-toupper(sigGene)
	for(cluster in c(0:46)){
  		###Get the result of ident_1
  		list1<-read.csv(paste("../BM.0.8_cluster_", cluster, ".csv", sep=""), header=T, row.names=1)
  		list2<-list1[list1[, "avg_log2FC"]>0.2 & list1[, "p_val"]<cutoffcluster, ]
  		clusGene<-toupper(rownames(list2))
  		clusGene<-intersect(clusGene, overGenes)
  		##Let us write the result
		cat(geneSets[[sig]]@setName, file=zz)
		cat("\t", file=zz)
		cat(cluster, file=zz)
		cat("\t", file=zz)
		overlappinggene<-intersect(sigGeneTop, clusGene)
		cat(length(overlappinggene), file=zz)
		cat("\t|", file=zz)
		for(gene in overlappinggene){
			cat(gene, file=zz)
			cat("|", file=zz)
		}
		cat("\t", file=zz)
		cat(length(intersect(sigGeneTop, overGenes)), file=zz)
		cat("\t", file=zz)
		cat(length(clusGene), file=zz)
		cat("\t", file=zz)
		cat(length(overGenes), file=zz)
		cat("\n", file=zz)
		}
}
close(zz)
