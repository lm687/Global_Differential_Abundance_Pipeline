source activate snakemake-globalDA
snakemake --cluster "sbatch -p general -t 12:00:00 --cores 1" --jobs 40 --printshellcmds
conda deactivate

