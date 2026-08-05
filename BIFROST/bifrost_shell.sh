#!/bin/bash
#
#SBATCH -p long
#SBATCH --mem=1G
#SBATCH --job-name=bifrost
#SBATCH --cpus-per-task=8
#SBATCH -o /home/rfelice/me/bifrost/bifrost1_%A.out
#SBATCH -e /home/rfelice/me/bifrost/bifrost1_%A.err
#SBATCH --mail-user=r.felice@nhm.ac.uk
#SBATCH --mail-type=ALL

source /mnt/shared/scratch/rfelice/apps/conda/etc/profile.d/conda.sh

conda activate r_r

Rscript --no-save /home/rfelice/me/bifrost/scripts/Bifrost_analysis_cluster.R
