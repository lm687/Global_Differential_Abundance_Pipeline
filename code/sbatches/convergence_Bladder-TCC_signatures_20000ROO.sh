#!/bin/sh
#SBATCH --job-name=Stan_Rsbatches/convergence_Bladder-TCC_signatures_20000ROO.sh
#SBATCH --out=sbatches/convergence_Bladder-TCC_signatures_20000ROO.out
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

 ~/.conda/envs/rstan_env/bin/Rscript --vanilla 3_analysis/assess_convergence.R  --file_posterior ../data/inference/Bladder-TCC_signatures_20000ROO.RData

conda deactivate
