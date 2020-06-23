#!/bin/sh
#SBATCH --job-name=Stan_ROvary-AdenoCA_signatures
#SBATCH --out=sbatches/PPC_Ovary-AdenoCA_signatures.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:00
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/posterior_predictive_checks.R --files_posteriors '../data/inference/Ovary-AdenoCA_signatures_20000_DMROO.RData ../data/inference/Ovary-AdenoCA_signatures_20000_LNMROO.RData ../data/inference/Ovary-AdenoCA_signatures_20000_MROO.RData'

conda deactivate
