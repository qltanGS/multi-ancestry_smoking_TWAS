library(optparse)
library(data.table)
library(tidyr)
library(dbplyr)
library(RSQLite)

#./summary_stat_prediction/geuvadis/pred_exp/
score_list <- dir("./summary_stat_prediction/geuvadis/individule/pred_expr/pred/", pattern = "*.profile$",recursive = T)

pred_exp <- list()
for (i in 1:length(score_list)){
    tmp_score <- fread(paste0("./summary_stat_prediction/geuvadis/individule/pred_expr/pred/",score_list[i])) %>% as.data.frame()
    tmp_score <- tmp_score[,c(2,6)]
    colnames(tmp_score)[2] <- gsub(".profile","",score_list[i])
    rownames(tmp_score) <- tmp_score[,1]
    pred_exp[[i]] <- tmp_score[,2]
}
pred_exp <- do.call(cbind, pred_exp)

tmp <- fread("./summary_stat_prediction/geuvadis/individule/pred_expr/pred/ENSG00000138073.profile") %>% as.data.frame()

rownames(pred_exp) <- tmp$IID
colnames(pred_exp) <- gsub(".profile","",score_list)
pred_exp <- as.data.frame(t(pred_exp))
pred_exp <- cbind(rownames(pred_exp),pred_exp)
colnames(pred_exp)[1] <- "gene_id"
#./summary_stat_prediction/geuvadis/individule/pred_expr/individule_pred_exp.txt
#./summary_stat_prediction/geuvadis/prscs_weights/prscs_pred_exp.txt
fwrite(pred_exp, row.names=F, col.names=T, quote=F, sep='\t', file="./summary_stat_prediction/geuvadis/individule/pred_expr/individule_pred_exp.txt")


