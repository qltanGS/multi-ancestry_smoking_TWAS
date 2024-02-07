library(optparse)
library(data.table)
library(tidyr)
library(dbplyr)
library(RSQLite)


#file list to import
weights_list <- dir("./summary_stat_prediction/geuvadis/prscs_weights/weights/", pattern = "*_weight.txt$",recursive = T)
covs_list <- dir("./summary_stat_prediction/geuvadis/prscs_weights/cov/", pattern = "*_covariance.txt$",recursive = T)

#rbind weight files
weights <- list()
for (i in 1:length(weights_list)){
  tmp_i <- fread(paste0("./summary_stat_prediction/geuvadis/prscs_weights/weights/",weights_list[i])) %>% as.data.frame()
  tmp_i$n_snp_used <- nrow(tmp_i)
  weights[[i]] <- tmp_i
}
weights <- do.call(rbind, weights)
fwrite(weights, col.names=T, row.names=F, quote=F, sep='\t',file=paste0("./summary_stat_prediction/geuvadis/prscs_weights/db/","geuvadis_weights.txt"))

#for predixcan cov file
covs <- list()
for (i in 1:length(covs_list)){
  tmp_i <- fread(paste0("./summary_stat_prediction/geuvadis/prscs_weights/cov/",covs_list[i])) %>% as.data.frame()
  covs[[i]] <- tmp_i
}
covs <- do.call(rbind, covs)
fwrite(covs, col.names=T, row.names=F, quote=F, sep='\t',
       file = paste0("./summary_stat_prediction/geuvadis/prscs_weights/db/geuvadis_cov.txt"))


genome_info = read.table('./anno/gencode/37/gencode.v32.GRCh37.txt', header = T,stringsAsFactors = F)
genome_info$gene_pos = (genome_info$left+genome_info$right)/2

all_info = merge(weights, genome_info[,c('geneid','genename','gene_pos','left','right')], by.x="GENE", by.y="geneid")  #all 37
all_info$model_id = sapply(all_info$GENE, function(x) as.numeric(sub('ENSG','',x)))

#gene list
gene_list <- all_info[,c('GENE','genename','model_id')]
gene_list <- gene_list[!duplicated(gene_list$genename),]

#FOCUS db file
print("INFO -------- generate weight file")
file_1_weight <- data.frame("id", all_info[,c('RSID','CHR','POS','A1','A2','Weight')],
                            std_error = NA, all_info[,c('model_id')])
colnames(file_1_weight) <- c("id","snp","chrom","pos","effect_allele","alt_allele","effect","std_error","model_id")
file_1_weight$id <- 1:nrow(file_1_weight)

print("INFO -------- generate model file")
file_2_model <- data.frame(id = unique(file_1_weight$model_id),
                           inference = "PRSCS",
                           ref_id = "1",
                           mol_id = unique(file_1_weight$model_id))

print("INFO -------- generate modelattribute file")
file_3_modelattribute <- data.frame(matrix(ncol=4,nrow=2*nrow(gene_list)))
colnames(file_3_modelattribute) = c('id','attr_name','value','model_id')
file_3_modelattribute[,'id'] <- 1:nrow(file_3_modelattribute)
file_3_modelattribute[,'attr_name'] <- ifelse(file_3_modelattribute[,1] %% 2 == 0,"cv.R2.pval","cv.R2")
file_3_modelattribute[,'value'] <- NA
file_3_modelattribute[seq(1,nrow(file_3_modelattribute),2),'model_id'] <- gene_list$model_id
file_3_modelattribute[seq(0,nrow(file_3_modelattribute),2),'model_id'] <- gene_list$model_id

print("INFO -------- generate molecularfeature file")
file_4_molecularfeature <- data.frame(all_info[,c('model_id','GENE')],
                                      ens_tx_id = "9606",
                                      all_info[,c('genename')],
                                      type = "NA",
                                      all_info[,c('CHR','left','right')])
colnames(file_4_molecularfeature) <- c('id',"ens_gene_id","ens_tx_id","mol_name","type","chrom","tx_start","tx_stop")
file_4_molecularfeature <- file_4_molecularfeature[!duplicated(file_4_molecularfeature$ens_gene_id),]

print("INFO -------- generate refpanel file")
file_5_refpanel <- data.frame(matrix(ncol=4,nrow=1))
file_5_refpanel[,1] <- "id"
file_5_refpanel[,2] <- "MetaBrain"
file_5_refpanel[,3] <- "Brain"
file_5_refpanel[,4] <- "rnaseq"
colnames(file_5_refpanel) <- c("id","ref_name","tissue","assay")
file_5_refpanel$id <- rownames(file_5_refpanel)

weight <- file_1_weight
model <- file_2_model
modelattribute <- file_3_modelattribute
molecularfeature <- file_4_molecularfeature 
refpanel <- file_5_refpanel

db<-dbConnect(RSQLite::SQLite(), paste0(opt$main_dir,"alpha_",alpha,"_",pvalue_cutoff,"_MA_FOCUS.db"))
dbWriteTable(db, "weight", weight)
dbWriteTable(db, "molecularfeature", molecularfeature)
dbWriteTable(db, "refpanel", refpanel)
dbWriteTable(db, "modelattribute", modelattribute)
dbWriteTable(db, "model", model)
dbDisconnect(db)


#for predixcan db file
predic_weight <- all_info[,c('RSID','GENE','Weight','A2','A1')]
colnames(predic_weight) <- c('rsid','gene','weight',"ref_allele",'eff_allele')

gene_info <- all_info[!duplicated(all_info$GENE),]

predic_extra = data.frame(gene = gene_info[,'GENE'],
                          genename = gene_info[,'genename'],
                          pred.perf.R2 = NA,
                          n.snps.in.model = gene_info[,'n_snp_used'],
                          pred.perf.pval = NA,
                          pred.perf.qval = NA)

pred_db<-dbConnect(RSQLite::SQLite(), paste0("./summary_stat_prediction/geuvadis/prscs_weights/db/","geuvadis_prscs.db"))
weights <- predic_weight
extra <- predic_extra
dbWriteTable(pred_db, "weights", weights)
dbWriteTable(pred_db, "extra", extra)
dbDisconnect(pred_db)
