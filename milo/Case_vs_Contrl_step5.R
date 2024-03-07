###we want to add cell types and HD PT
rm(list = ls())
load("da_results.RData")
load("../plotResultAliceAnnotation/metaData_AliceAnnotation_CD4CD8.RData")
da_results<-da_results[da_results$FDR<0.3,]
metaData$status<-"PT"
metaData$status[grepl("_HD", metaData$sampleID)]<-"HD"
da_results$statuslist<-""
da_results$statusHDCount[ii]<-0
da_results$statusPTCount[ii]<-0
da_results$celltypelist<-""
for(ii in 1:dim(da_results)[1]){
  cat(ii);cat("\n")
  metaData_DA<-metaData[rownames(metaData) %in% strsplit(da_results$nhoods_node_names[ii], ":")[[1]],]
  da_results$statuslist[ii]<-paste(sort(metaData_DA$status), collapse = ":")
  da_results$celltypelist[ii]<-paste(sort(metaData_DA$celltype), collapse = ":")
  metaData_DA_HD<-metaData_DA[metaData_DA$status=="HD",]
  metaData_DA_PT<-metaData_DA[metaData_DA$status=="PT",]
  da_results$statusHDCount[ii]<-dim(metaData_DA_HD)[1]
  da_results$statusPTCount[ii]<-dim(metaData_DA_PT)[1]
}

save(da_results, file = "da_results_withcelltype_status.RData")
write.csv(da_results, file = "da_results_withcelltype_status.csv", quote = FALSE)
