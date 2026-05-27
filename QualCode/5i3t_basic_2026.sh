#!/bin/bash
#SBATCH --job-name=5i3t_basic
#SBATCH --account=yx18-ic
#SBATCH --partition=IllinoisComputes
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-user=hk33@illinois.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load R/4.5.1
export R_LIBS_USER="/u/hk33/R/x86_64-redhat-linux-gnu-library/4.5"

Rscript 5i3t_basic_2026.R
