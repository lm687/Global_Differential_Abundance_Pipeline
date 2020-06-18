#!/bin/sh
#SBATCH --job-name=Stan_REso-AdenoCA_nucleotidesubstitution1
#SBATCH --out=sbatches/PPC_Eso-AdenoCA_nucleotidesubstitution1.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:00
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/posterior_predictive_checks.R --files_posterior 'Eso-AdenoCA_nucleotidesubstitution1_20000_DMROO.RData'

conda deactivate