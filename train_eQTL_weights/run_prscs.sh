#!/bin/bash
#SBATCH --mail-user=  #############################
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=100G
#SBATCH --job-name=prscs_geuvadis
#SBATCH --array=1-29




ml GCC OpenMPI R Anaconda3/2022.05


Rscript ./summary_stat_prediction/geuvadis/prscs_eqtl.r ${SLURM_ARRAY_TASK_ID} 

