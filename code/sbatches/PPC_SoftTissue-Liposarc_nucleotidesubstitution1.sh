#!/bin/sh
#SBATCH --job-name=Stan_RSoftTissue-Liposarc_nucleotidesubstitution1
#SBATCH --out=sbatches/PPC_SoftTissue-Liposarc_nucleotidesubstitution1.out
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:00
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/posterior_predictive_checks.R --files_posteriors '../data/inference/SoftTissue-Liposarc_nucleotidesubstitution1_20000_DMROO.RData'

conda deactivate
