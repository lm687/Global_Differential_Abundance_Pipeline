#!/bin/sh
#SBATCH --job-name=Stan_R../data/roo/SoftTissue-Liposarc_signatures_ROO.RDS
#SBATCH --out=sbatches/SoftTissue-Liposarc_signatures_ROO.out
#SBATCH --cpus-per-task=1
#SBATCH --time=9:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate rstan_env

 ~/.conda/envs/rstan_env/bin/Rscript --vanilla 2_inference/fit_PCAWG.R --cancertype ../data/roo/SoftTissue-Liposarc --typedata signatures --infile ../data/roo/SoftTissue-Liposarc_signatures_ROO.RDS --output ../data/inference/SoftTissue-Liposarc_signatures_ROO.RData --iterations 10000 --model M

conda deactivate
