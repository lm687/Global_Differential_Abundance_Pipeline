#!/bin/sh
#SBATCH --job-name=Stan_RBiliary-AdenoCA_signatures
#SBATCH --out=sbatches/sbatches/PPC_Biliary-AdenoCA_signatures.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:00
#SBATCH --no-requeue
#SBATCH -p general

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/posterior_predictive_checks.R --files_posterior 'Biliary-AdenoCA_signatures_20000_MROO.RData'

conda deactivate
