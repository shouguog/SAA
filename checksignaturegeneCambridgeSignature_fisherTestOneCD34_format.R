###This is used to reformat the Fisher result

rm(list=ls())

zz1 <- file("signaureGenesOverlapFisherFormatOneCD34.txt", "w")  # open an output file connection
cat("cluster\tcelltype\n", file=zz1);
result<-read.table("signaureGenesOverlapFisherOneCD34.txt",header=T) 
for(cluster in c(0:46)){
  ####for each cluster we record top 3
  result_cluster<-result[result[, "cluster"]==cluster,]
  result_cluster<-result_cluster[order(result_cluster[, "pvalue"]),]
  cat(cluster, file=zz1)
  cat("\t", file=zz1)
  defined=0;
  for(ii in 1:3){
    if(result_cluster[ii, "pvalue"]<1e-5){
      defined=1;
      cat(as.character(result_cluster[ii, "signature"]), file=zz1)
      cat("(", file=zz1)
      cat(result_cluster[ii, "pvalue"], file=zz1)
      cat(")|", file=zz1)
    }
  }
  if(defined==0){
    cat("No type Defination", file=zz1)
  }
  cat("\n", file=zz1)
}
close(zz1)


