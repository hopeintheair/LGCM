#!/bin/bash
#SBATCH --job-name=lgcm_sim
#SBATCH --account=yx18-ic
#SBATCH --partition=IllinoisComputes
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-user=hk33@illinois.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-180                                      
#SBATCH --output=/projects/illinois/educ/ep/yx18/Hahyeong/Dissertation/2026/logs_newer/job%a.out
#SBATCH --error=/projects/illinois/educ/ep/yx18/Hahyeong/Dissertation/2026/logs_newer/job%a.err

PROJECT_DIR="/projects/illinois/educ/ep/yx18/Hahyeong/Dissertation/2026"

export R_LIBS_USER="/u/hk33/R/x86_64-redhat-linux-gnu-library/4.5"
module load R/4.5.1

mkdir -p "${PROJECT_DIR}/logs_newer" "${PROJECT_DIR}/results_newer"

cd "${PROJECT_DIR}" || exit 1

Rscript "${PROJECT_DIR}/run_newer.R" "${SLURM_ARRAY_TASK_ID}"