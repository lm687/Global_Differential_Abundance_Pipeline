#!/bin/sh
#SBATCH --job-name=Stan_R../data/roo/Kidney-RCC.papillary_nucleotidesubstitution1_ROO.RDS
#SBATCH --out=sbatches/Kidney-RCC.papillary_nucleotidesubstitution1_ROO.out
#SBATCH --cpus-per-task=1
#SBATCH --time=9:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env

 ~/.conda/envs/rstan_env/bin/Rscript --vanilla 2_inference/fit_PCAWG.R --cancertype ../data/roo/Kidney-RCC.papillary --typedata nucleotidesubstitution1 --infile ../data/roo/Kidney-RCC.papillary_nucleotidesubstitution1_ROO.RDS --output ../data/inference/Kidney-RCC.papillary_nucleotidesubstitution1_20000ROO.RData --iterations 20000 --model M

conda deactivate
