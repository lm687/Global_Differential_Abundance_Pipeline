#!/bin/sh
#
#SBATCH --job-name=stanfullREDMDL
#SBATCH --output=stanfullREDMDL.out
#SBATCH --ntasks=4

Rscript running_stan2.R

