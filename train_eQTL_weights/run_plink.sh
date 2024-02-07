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
#SBATCH --array=1-21



ml GCC OpenMPI R Anaconda3/2022.05


Rscript ./summary_stat_prediction/geuvadis/plink_score.r ${SLURM_ARRAY_TASK_ID} 

