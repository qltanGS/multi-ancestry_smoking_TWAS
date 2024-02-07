library(MatrixEQTL)
library(data.table)
library(tidyr)
library(dbplyr)
library(optparse)

option_list = list(
  #general
  make_option("--chr", action="store", default=NA, type='character')
)
opt = parse_args(OptionParser(option_list=option_list))


errorCovariance = numeric()
useModel = modelLINEAR

snp_info  <- fread("./summary_stat_prediction/PRScs/liver_data/prs_full_snp/eqtl/snp_pos_info.txt") %>% as.data.frame()
gene_info <- fread("./summary_stat_prediction/PRScs/liver_data/prs_full_snp/eqtl/gene_pos_info.txt") %>% as.data.frame()

snp = SlicedData$new()
snp$fileDelimiter = "\t"     # the TAB character
snp$fileOmitCharacters = "NA" # denote missing values;
snp$fileSkipRows = 1          # one row of column labels
snp$fileSkipColumns = 1       # one column of row labels
snp$fileSliceSize = 2000     # read file in pieces of 2,000 rows
snp$LoadFile( paste0("./summary_stat_prediction/geuvadis/eqtl/matrixeqtl_input/input_snp_chr",opt$chr,".txt") )

exp = SlicedData$new()
exp$fileDelimiter = "\t"      # the TAB character
exp$fileOmitCharacters = "NA" # denote missing values;
exp$fileSkipRows = 1          # one row of column labels
exp$fileSkipColumns = 1      # one column of row labels
exp$fileSliceSize = 2000      # read file in pieces of 2,000 rows
exp$LoadFile("./summary_stat_prediction/geuvadis/eqtl/matrixeqtl_input/input_exp.txt")

cov = SlicedData$new()
cov$LoadFile("./summary_stat_prediction/geuvadis/eqtl/matrixeqtl_input/input_cov.txt")

pvOutputThreshold_cis = 0.05;
pvOutputThreshold_tra = 1e-10000;
cisDist = 1e6;


wb_eqtl = Matrix_eQTL_main(
  snps = snp, 
  gene = exp, 
  cvrt = cov,
  output_file_name = './summary_stat_prediction/geuvadis/eqtl/geuvadis_eQTL_trans.txt',
  pvOutputThreshold     = 1e-1000,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = paste0("./summary_stat_prediction/geuvadis/eqtl/geuvadis_eQTL_cis_chr",opt$chr,".txt"),
  pvOutputThreshold.cis = 0.05,
  snpspos = snp_info, 
  genepos = gene_info,
  cisDist = cisDist,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);


