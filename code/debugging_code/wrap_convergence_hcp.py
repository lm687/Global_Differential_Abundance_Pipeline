


import os
  
#print(os.listdir("../data/roo/"))

nits = 20000
for i in ['../data/inference/'+x for x in os.listdir("../data/inference/") if 'ROO' in x]:

    script = "~/.conda/envs/rstan_env/bin/Rscript --vanilla 3_analysis/assess_convergence.R  --file_posterior {}".format(i)

    fileout="sbatches/convergence_"+os.path.basename(i).replace(".RData", ".sh")
    a = open(fileout, "w")
    a.write('''#!/bin/sh
#SBATCH --job-name=Stan_R{}
#SBATCH --out={}
#SBATCH --cpus-per-task=1
#SBATCH --time={}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -A MARKOWETZ-SL3-CPU
#SBATCH -p skylake

module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
source activate {}

 {}

conda deactivate
'''.format(fileout, fileout.replace(".sh", ".out"), "2:00:00", "rstan_env_analysis", script))
    a.close()
    print("sbatch "+fileout)



