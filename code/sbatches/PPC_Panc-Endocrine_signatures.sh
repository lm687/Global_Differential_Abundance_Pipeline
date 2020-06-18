#!/bin/sh
#SBATCH --job-name=Stan_RPanc-Endocrine_signatures
#SBATCH --out=sbatches/PPC_Panc-Endocrine_signatures.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:00
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/posterior_predictive_checks.R --files_posterior 'Panc-Endocrine_signatures_20000_DMROO.RData Panc-Endocrine_signatures_20000_LNMROO.RData Panc-Endocrine_signatures_ROO.RData'

conda deactivate