#!/bin/sh
# properties = {"type": "single", "rule": "simulate_for_comparison_multiple", "local": false, "input": ["../data/assessing_models_simulation/GenerationInoRE"], "output": ["../data/assessing_models_simulation/datasets/multiple_GenerationInoRE_50_100_80_6_0.01_NA_NA_NA_dataset0.RDS"], "wildcards": {"datasetgeneration": "GenerationInoRE", "n": "50", "nlambda": "100", "lmbda": "80", "d": "6", "beta_intensity": "0.01", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "0"}, "params": {"datasetgeneration": "GenerationInoRE", "n": "50", "nlambda": "100", "lmbda": "80", "d": "6", "beta_intensity": "0.01", "itnum": "0", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA"}, "log": [], "threads": 1, "resources": {}, "jobid": 1355, "cluster": {}}
 cd /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code && \
PATH='/scratcha/fmlab/morril01/software/miniconda3/envs/snakemake-globalDA/bin':$PATH /scratcha/fmlab/morril01/software/miniconda3/envs/snakemake-globalDA/bin/python3.9 \
-m snakemake ../data/assessing_models_simulation/datasets/multiple_GenerationInoRE_50_100_80_6_0.01_NA_NA_NA_dataset0.RDS --snakefile /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.m52_au9w ../data/assessing_models_simulation/GenerationInoRE --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules simulate_for_comparison_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.m52_au9w/1355.jobfinished || (touch /mnt/scratcha/fmlab/morril01/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.m52_au9w/1355.jobfailed; exit 1)

