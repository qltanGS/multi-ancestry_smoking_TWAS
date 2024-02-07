library(data.table)
library(tidyr)
library(dbplyr)
library(ggplot2)
library(patchwork)

prscs <- fread("./summary_stat_prediction/geuvadis/prscs_weights/prscs_pred_exp.txt") %>% as.data.frame()
indiv <- fread("./summary_stat_prediction/geuvadis/individule/pred_expr/individule_pred_exp.txt") %>% as.data.frame()
expr <- fread("./projects/cross_tissue/rep/data/exp/exp_wb.txt") %>% as.data.frame()

rownames(prscs) <- prscs$gene_id
rownames(indiv) <- indiv$gene_id
rownames(expr) <- expr$geneid

prscs <- prscs[,-1]
indiv <- indiv[,-1]
expr <- expr[,-1]

gene_id <- Reduce(intersect, list(rownames(prscs),
                                    rownames(indiv),
                                    rownames(expr)))
sample_id   <- Reduce(intersect, list(colnames(prscs),
                                    colnames(indiv),
                                    colnames(expr)))

prscs <- prscs[gene_id, sample_id]
indiv <- indiv[gene_id, sample_id]
expr <- expr[gene_id, sample_id]

summ_pred <- as.data.frame(t(prscs))
indi_pred <- as.data.frame(t(indiv))
observed  <- as.data.frame(t(expr))

#calculate r^2 of summ_pred VS. indi_pred
summ_vs_indi = data.frame(matrix(nrow = length(gene_id), ncol = 2))
for (i in 1:length(gene_id)){
    cal1 <- summary(lm(summ_pred[,i]~indi_pred[,i]))
    summ_vs_indi[i,1] <- cal1$r.squared

    cal2 <- cal1$fstatistic
    summ_vs_indi[i,2] <- pf(cal2["value"], cal2["numdf"], cal2["dendf"], lower.tail = FALSE)
}
colnames(summ_vs_indi) <- c("r_squared","p_val")

#calculate r^2 of summ_pred VS. Observed
summ_vs_obser = data.frame(matrix(nrow = length(gene_id), ncol = 2))
for (i in 1:length(gene_id)){
    cal1 <- summary(lm(observed[,i]~summ_pred[,i]))
    summ_vs_obser[i,1] <- cal1$r.squared

    cal2 <- cal1$fstatistic
    summ_vs_obser[i,2] <- pf(cal2["value"], cal2["numdf"], cal2["dendf"], lower.tail = FALSE)
}
colnames(summ_vs_obser) <- c("r_squared","p_val")

#calculate r^2 of indi_pred VS. Observed
indi_vs_obser = data.frame(matrix(nrow = length(gene_id), ncol = 2))
for (i in 1:length(gene_id)){
    cal1 <- summary(lm(observed[,i]~indi_pred[,i]))
    indi_vs_obser[i,1] <- cal1$r.squared

    cal2 <- cal1$fstatistic
    indi_vs_obser[i,2] <- pf(cal2["value"], cal2["numdf"], cal2["dendf"], lower.tail = FALSE)
}
colnames(indi_vs_obser) <- c("r_squared","p_val")

#get datafram of r^2_summ_vs_obser & r^2_indi_vs_obser
comp_df <- data.frame(summ_vs_obser$r_squared,indi_vs_obser$r_squared )
colnames(comp_df) <- c("Ob_summ","Ob_indi")


p1 <- ggplot(comp_df, aes(x=Ob_summ, y=Ob_indi)) +
      geom_point() +
      geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
      ggtitle("R^2_summ_Obser vs. R^2_indi_Obser") +
      theme_classic()
#ggsave(p1, filename = "./summary_stat_prediction/geuvadis/compared_r_squared.png")


#Correlation coefficient between individule and summary
cor_indi_summ <- list()
for (i in 1:length(gene_id)){
    tmp_cor <- cor.test(summ_pred[,i],indi_pred[,i])
    cor_indi_summ[[i]] <- tmp_cor$estimate
}
cor_indi_summ <- as.numeric(unlist(cor_indi_summ))
plot_df <- data.frame(colnames(summ_pred),cor_indi_summ)
p2 <- ggplot(plot_df, aes(x = cor_indi_summ)) +
        geom_histogram()+
        ggtitle("r of summary_pred VS. individule_pred")
#ggsave(p2, file="./summary_stat_prediction/geuvadis/summ_indiv_cor.png")

pp <- p1|p2
ggsave(pp, filename = "./summary_stat_prediction/geuvadis/compare_geuvadis.png", width = 14)
save(summ_pred, indi_pred, observed, sample_id, gene_id, summ_vs_indi, summ_vs_obser,indi_vs_obser,plot_df,cor_indi_summ,  file = "./summary_stat_prediction/geuvadis/compare_work.Rdata")
