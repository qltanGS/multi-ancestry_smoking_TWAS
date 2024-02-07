#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 X11/.20190311 R/3.6.0
args<-as.numeric(commandArgs(TRUE))

library(RSQLite)
library(data.table)
library(foreach) 
library(doParallel)
library(tidyr)
library(dbplyr)


#获取基因列表
geuvadis_plink <- "./projects/cross_tissue/rep/data/geno/geuvadis421"

#weights_list <- dir("./summary_stat_prediction/geuvadis/prscs_weights/weights/", pattern = "*_weight.txt$",recursive = T)
#gene_list <- gsub("_weight.txt","",weights_list)
con <-dbConnect(RSQLite::SQLite(),dbname='./summary_stat_prediction/geuvadis/individule/db/GEUVADIS_LCLs.db')
weights=dbReadTable(con,'weights')
gene_list <- weights$gene[!duplicated(weights$gene)]

plink_score <- function(gene_id){
    #tmp_weight <- fread(paste0("./summary_stat_prediction/geuvadis/prscs_weights/weights/",gene_id,"_weight.txt")) %>% as.data.frame()
    tmp_weight <- weights[weights$gene == gene_id,]
    write.table(tmp_weight, col.names=T, row.names=F, quote=F, sep='\t', file=paste0("./summary_stat_prediction/geuvadis/individule/pred_expr/weights/", gene_id, "_weight.txt"))
    write.table(tmp_weight$rsid, row.names=F, col.names=F, quote=F, file=paste0("./summary_stat_prediction/geuvadis/individule/pred_expr/tmp/",gene_id,"_rs_id.txt"))

    plink_cmd <- paste0("plink --bfile ./projects/cross_tissue/rep/data/geno/geuvadis421 --score ./summary_stat_prediction/geuvadis/individule/pred_expr/weights/",
                        gene_id,"_weight.txt"," 1 5 3 --extract ./summary_stat_prediction/geuvadis/individule/pred_expr/tmp/",gene_id,"_rs_id.txt --out ./summary_stat_prediction/geuvadis/individule/pred_expr/pred/",
                        gene_id)
    print(plink_cmd)
    system(plink_cmd,ignore.stdout=T,ignore.stderr=T,wait=T)
}
parallel = TRUE
n = args
i_start <- (n-1)*500+1
if (n < 7){
  i_end <- i_start+499
}else{
  i_end <- 3218
}
print(i_start)

#run
if(length(unique(gene_list))>2){
  if(parallel){
    registerDoParallel(cores = detectCores())
    #registerDoParallel(cores = max(detectCores()-1, 1)) for vgi01
    foreach(i = 1:length(gene_list)) %dopar% {
      plink_score(gene_list[i])
      
      stopImplicitCluster()
    }
  }else{
    for (i in 1:length(gene_list)){
      print(paste0(i,'/',gene_list[i]))
      plink_score(gene_list[i])
    }
  }
}else{
  print('no more genes')
}

