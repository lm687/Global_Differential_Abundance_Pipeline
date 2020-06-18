import os

folder_inference = '../data/inference/'
subfiles = set([ '_'.join(x.split('_')[0:2]) for x in os.listdir(folder_inference) if 'RData' in x ])

for s in subfiles:
    all_files = ' '.join([folder_inference+j for j in os.listdir(folder_inference) if str(s) in j])

    script = "Rscript --vanilla 3_analysis/posterior_predictive_checks.R --files_posteriors '{}'".format(all_files)

    fileout="sbatches/PPC_"+str(s)+".sh"
    a = open(fileout, "w")
    a.write('''#!/bin/sh
#SBATCH --job-name=Stan_R{}
#SBATCH --out={}
#SBATCH --cpus-per-task=1
#SBATCH --time={}
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate {}

{}

conda deactivate
'''.format(str(s), fileout.replace(".sh", ".out"), "0:40:00", "rstan_env_analysis", script))
    a.close()

    print("sbatch "+fileout)    

