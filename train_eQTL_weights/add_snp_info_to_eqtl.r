library(dplyr)
library(tidyverse)
library(data.table)

eqtl_raw <- fread("./summary_stat_prediction/geuvadis/eqtl/geuvadis_eQTL_cis.txt") %>% as.data.frame()
colnames(eqtl_raw) <- c('SNP','gene','beta','t_stat','pval','FDR')

freq <- fread("./projects/cross_tissue/rep/data/geno/geuvadis421.frq") %>% as.data.frame()
freq <- freq[match(eqtl_raw$SNP, freq$SNP),]
identical(eqtl_raw$SNP, freq$SNP)

bim <- fread("./projects/cross_tissue/rep/data/geno/geuvadis421.bim", header=F) %>% as.data.frame()
colnames(bim) <- c('chr','rsid','a',"pos",'eff_allele',"ref_allele")
bim <- bim[bim$rsid %in% eqtl_raw$SNP,]
bim <- bim[match(eqtl_raw$SNP, bim$rsid),]
identical(eqtl_raw$SNP, bim$rsid)
se <- sqrt(((eqtl_raw$beta)^2)/qchisq(eqtl_raw$pval,1,lower.tail=F))

new_eqtl <- data.frame(eqtl_raw$gene,eqtl_raw$SNP,bim$chr, bim$pos,bim$eff_allele,bim$ref_allele, freq$MAF,'N_sample',se,eqtl_raw[,c(3:6)])
colnames(new_eqtl) <- c("Gene","rsid","chr","pos",'eff_allele',"ref_allele","MAF",'N_sample','SE','beta','t_stat','pval','FDR')
new_eqtl$N_sample <- 421

fwrite(new_eqtl, row.names=F, col.names=T, quote=F, sep='\t',file = "./summary_stat_prediction/geuvadis/eqtl/geuvadis_eQTL_with_info.txt")
