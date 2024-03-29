#!/bin/sh
#SBATCH --job-name=Stan_Rsbatches/convergence_SoftTissue-Liposarc_signatures_20000_LNMROO.sh
#SBATCH --out=sbatches/convergence_SoftTissue-Liposarc_signatures_20000_LNMROO.out
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

 Rscript --vanilla 3_analysis/assess_convergence.R  --file_posterior ../data/inference/SoftTissue-Liposarc_signatures_20000_LNMROO.RData

conda deactivate
