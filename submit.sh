#!/bin/bash
# submit.sh — SLURM array job
#SBATCH --job-name=lgcm_sim
#SBATCH --array=1-48                # 48 conditions (= nrow of conditions.csv)
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00           
#SBATCH --output=logs/job%a.out
#SBATCH --error=logs/job%a.err

# ── Edit these for your cluster ──────────────────────────────────────────────
export PROJECT_DIR=$HOME/LGCM       # path to this folder on the cluster
module load R/4.3.1                 # check available: module avail R

# ── Run ──────────────────────────────────────────────────────────────────────
mkdir -p $PROJECT_DIR/logs $PROJECT_DIR/results

Rscript $PROJECT_DIR/run_one.R $SLURM_ARRAY_TASK_ID
