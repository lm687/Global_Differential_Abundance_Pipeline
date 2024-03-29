#!/bin/sh
#SBATCH --job-name=Stan_R../data/roo/Liver-HCC_nucleotidesubstitution3_ROO.RDS
#SBATCH --out=sbatches/Liver-HCC_nucleotidesubstitution3_ROO_LNM.out
#SBATCH --cpus-per-task=1
#SBATCH --time=9:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env

 ~/.conda/envs/rstan_env/bin/Rscript --vanilla 2_inference/fit_PCAWG.R --cancertype Liver-HCC --typedata nucleotidesubstitution3 --infile ../data/roo/Liver-HCC_nucleotidesubstitution3_ROO.RDS --output ../data/inference/Liver-HCC_nucleotidesubstitution3_20000_LNMROO.RData --iterations 20000 --model LNM

conda deactivate
