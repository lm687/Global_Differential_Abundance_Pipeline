#!/bin/sh
#SBATCH --job-name=analyse_betas
##SBATCH --out=sbatches/analyse_betas$RANDOM.out
#SBATCH --cpus-per-task=1
#SBATCH --time=9:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/extract_betas.R --model_and_feature LNM_signatures --grep_posteriors signatures_15000_LNMROO.RData


conda deactivate



