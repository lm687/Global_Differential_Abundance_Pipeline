#!/bin/sh
#SBATCH --job-name=Stan_RKidney-ChRCC_signatures
#SBATCH --out=sbatches/PPC_Kidney-ChRCC_signatures.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:00
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/posterior_predictive_checks.R --files_posterior 'Kidney-ChRCC_signatures_ROO.RData Kidney-ChRCC_signatures_20000_LNMROO.RData Kidney-ChRCC_signatures_20000_DMROO.RData Kidney-ChRCC_signatures_20000_MROO.RData'

conda deactivate
