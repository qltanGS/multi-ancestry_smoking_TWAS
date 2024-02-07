#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=plink_geuvadis
#SBATCH --array=1-7


ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 X11/.20190311 R/3.6.0


Rscript ./summary_stat_prediction/geuvadis/indi_plink_score.r ${SLURM_ARRAY_TASK_ID} 

