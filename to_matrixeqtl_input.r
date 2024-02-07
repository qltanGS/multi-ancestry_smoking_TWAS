library(dplyr)
library(tidyverse)
library(data.table)
library(optparse)

for (i in 1:26){
    dosage <- fread(paste0("./projects/cross_tissue/rep/data/geno/tmp/geuvadis421_chr", i,".raw")) %>% as.data.frame()
    dosage <- dosage[,-c(1,3,4,5,6)]
    colnames(dosage) <- as.character(sapply(colnames(dosage), function(x) strsplit(x, "_")[[1]][1]))
    rownames(dosage) <- dosage$IID
    dosage <- dosage[,-1]
    dosage <- as.data.frame(t(dosage))
    dosage <- cbind(rownames(dosage),dosage)
    colnames(dosage)[1] <- "snpid"
    write.table(dosage, col.names=T, row.names=F, quote=F, sep='\t', file = paste0("./summary_stat_prediction/geuvadis/eqtl/input_snp_chr",i,".txt"))
    #fwrite(dosage, col.names=T, row.names=F, quote=F, sep='\t',file=paste0("./summary_stat_prediction/geuvadis/eqtl/input_snp_chr",i,".txt"))
    gc()
}

cov_1 <- fread("./projects/cross_tissue/rep/data/geno/geuvadis421.fam", header=F) %>% as.data.frame()
cov_2 <- fread("./projects/cross_tissue/rep/data/geno/tmp/geuvadis421_pca_result.eigenvec", header=F) %>% as.data.frame()
cov <- data.frame(cov_2[,c(3:7)],cov_1$V5)
colnames(cov) <- c("pc1","pc2","pc3","pc4","pc5","sex")
rownames(cov) <- cov_1$V2
cov <- as.data.frame(t(cov))
cov <- cbind(rownames(cov),cov)
colnames(cov)[1] <- "id"
fwrite(cov, col.names=T, row.names=F, quote=F, sep='\t',file="./summary_stat_prediction/geuvadis/eqtl/input_cov.txt")