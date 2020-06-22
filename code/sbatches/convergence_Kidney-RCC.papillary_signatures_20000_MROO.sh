#!/bin/sh
#SBATCH --job-name=Stan_Rsbatches/convergence_Kidney-RCC.papillary_signatures_20000_MROO.sh
#SBATCH --out=sbatches/convergence_Kidney-RCC.papillary_signatures_20000_MROO.out
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

<<<<<<< HEAD
 ~/.conda/envs/rstan_env/bin/Rscript --vanilla 3_analysis/assess_convergence.R  --file_posterior ../data/inference/Kidney-RCC.papillary_signatures_20000_MROO.RData
=======
 Rscript --vanilla 3_analysis/assess_convergence.R  --file_posterior ../data/inference/Kidney-RCC.papillary_signatures_20000_MROO.RData
>>>>>>> 6b3ba668aa912fca359d134121e426c2ee2d138a

conda deactivate
