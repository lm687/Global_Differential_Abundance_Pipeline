#!/bin/sh
#SBATCH --job-name=sim_inference
#SBATCH --out=sbatches/sim_inference.out
#SBATCH --cpus-per-task=1
#SBATCH --time=9:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env_analysis

Rscript --vanilla 3_analysis/simulation_results.R

conda deactivate



