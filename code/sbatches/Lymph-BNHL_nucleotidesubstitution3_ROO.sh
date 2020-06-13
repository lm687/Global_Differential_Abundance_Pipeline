#!/bin/sh
#SBATCH --job-name=Stan_R../data/roo/Lymph-BNHL_nucleotidesubstitution3_ROO.RDS
#SBATCH --out=sbatches/Lymph-BNHL_nucleotidesubstitution3_ROO.out
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -p general

source activate rstan_env

Rscript --vanilla  Rscript --vanilla 2_inference/fit_PCAWG.R --cancertype ../data/roo/Lymph-BNHL --typedata nucleotidesubstitution3 --infile ../data/roo/Lymph-BNHL_nucleotidesubstitution3_ROO.RDS -- output ../data/inference/Lymph-BNHL_nucleotidesubstitution3_ROO.RData --iterations 10000 --model M

conda deactivate
