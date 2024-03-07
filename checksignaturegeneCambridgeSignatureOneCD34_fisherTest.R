###This is used to check the expression of clusters
rm(list=ls())
zz <- read.table("signaureGenesOverlapOneCD34.txt", header=T)
zz1 <- file("signaureGenesOverlapFisherOneCD34.txt", "w")  # open an output file connection
cat("signature\tcluster\toverlapNum\toverlapGene\tsignaturegene\tclustergene\tallsignaturegene\tpvalue\n", file=zz1)
for(ii in 1:dim(zz)[1]){
  for(kk in 1:7){
    cat(as.character(zz[ii,kk]), file=zz1)
    cat("\t", file=zz1)
  }
  Convictions <-  matrix(c(zz[ii,3], zz[ii,6]-zz[ii,3], zz[ii,5]-zz[ii,3], zz[ii,7]-zz[ii,5]-zz[ii,6]+zz[ii,3]),
                         nrow = 2,dimnames =list(c("Signature", "NoSignature"),c("ClusterOverExp", "ClusterNoOverExp")))
  Convictions
  cat(fisher.test(Convictions, alternative = "greater")$p.value, file=zz1)
  cat("\n", file=zz1)
}
close(zz1)