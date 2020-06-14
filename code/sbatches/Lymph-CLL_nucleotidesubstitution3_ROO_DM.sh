#!/bin/sh
#SBATCH --job-name=Stan_R../data/roo/Lymph-CLL_nucleotidesubstitution3_ROO.RDS
#SBATCH --out=sbatches/Lymph-CLL_nucleotidesubstitution3_ROO_DM.out
#SBATCH --cpus-per-task=1
#SBATCH --time=9:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env

 ~/.conda/envs/rstan_env/bin/Rscript --vanilla 2_inference/fit_PCAWG.R --cancertype Lymph-CLL --typedata nucleotidesubstitution3 --infile ../data/roo/Lymph-CLL_nucleotidesubstitution3_ROO.RDS --output ../data/inference/Lymph-CLL_nucleotidesubstitution3_20000_DMROO.RData --iterations 20000 --model DM

conda deactivate
