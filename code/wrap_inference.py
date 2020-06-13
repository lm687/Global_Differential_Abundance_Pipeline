


import os
  
#print(os.listdir("../data/roo/"))

for i in ['../data/roo/'+x for x in os.listdir("../data/roo/") if 'ROO' in x]:
    script = "Rscript --vanilla 2_inference/fit_PCAWG.R --cancertype {} --typedata {} --infile {} -- output {} --iterations 10000 --model M".format(i.split('_')[0], (i.split('_')[1]), i, i.replace("roo/", "inference/").replace(".RDS", ".RData"))

    fileout="sbatches/"+os.path.basename(i).replace(".RDS", ".sh")
    a = open(fileout, "w")
    a.write('''#!/bin/sh
#SBATCH --job-name=Stan_R{}
#SBATCH --out=sbatches/{}
#SBATCH --cpus-per-task=1
#SBATCH --time={}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lena.morrill@cruk.cam.ac.uk
#SBATCH --no-requeue
#SBATCH -p general

source activate {}

Rscript --vanilla  {}

conda deactivate
'''.format(i, os.path.basename(i).replace(".RDS", ".out"), "5:00:00", "rstan_env", script))
    a.close()

    print("sbatch "+fileout)


#Rscript --vanilla 2_inference/fit_PCAWG.R --cancertype Prost-AdenoCA --typedata signatures         --infile ../data/roo/Prost-AdenoCA_signatures_ROO.RDS --output ../data/inference/M_6000_Prost-AdenoCA_signatures.Rdata --niterations 6000 --model M



