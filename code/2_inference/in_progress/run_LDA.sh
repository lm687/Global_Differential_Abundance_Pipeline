#!/bin/sh
#SBATCH --job-name="LDA_simple"
#SBATCH --out="LDA_simple.out"
##SBATCH --cpus-per-task=4
##SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH -p skylake
#SBATCH -A MARKOWETZ-SL3-CPU

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate snakemake-globalDA
Rscript --vanilla LDA.R 
conda deactivate

