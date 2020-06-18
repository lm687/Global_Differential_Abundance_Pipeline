#!/bin/sh
#SBATCH --job-name=Stan_RBreast-LobularCA_signatures
#SBATCH --out=sbatches/PPC_Breast-LobularCA_signatures.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:00
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/posterior_predictive_checks.R --files_posterior 'Breast-LobularCA_signatures_20000_LNMROO.RData Breast-LobularCA_signatures_ROO.RData Breast-LobularCA_signatures_20000_MROO.RData Breast-LobularCA_signatures_20000_DMROO.RData'

conda deactivate
