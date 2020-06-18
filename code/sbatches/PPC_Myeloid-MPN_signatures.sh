#!/bin/sh
#SBATCH --job-name=Stan_RMyeloid-MPN_signatures
#SBATCH --out=sbatches/PPC_Myeloid-MPN_signatures.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:00
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/posterior_predictive_checks.R --files_posterior 'Myeloid-MPN_signatures_20000_DMROO.RData Myeloid-MPN_signatures_20000_LNMROO.RData Myeloid-MPN_signatures_ROO.RData Myeloid-MPN_signatures_20000_MROO.RData'

conda deactivate