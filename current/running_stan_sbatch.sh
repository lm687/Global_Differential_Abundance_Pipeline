#!/bin/sh
#
#SBATCH --job-name=stanfullREDMSL
#SBATCH --output=stanfullREDMSL.out
#SBATCH --ntasks=1

Rscript running_stan.R


