#!/bin/sh
# properties = {"type": "single", "rule": "simulate_for_comparison_multiple", "local": false, "input": ["../data/assessing_models_simulation/GenerationJnorm"], "output": ["../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_20_100_10_3_0.6_NA_NA_NA_dataset0.RDS"], "wildcards": {"datasetgeneration": "GenerationJnorm", "n": "20", "nlambda": "100", "lmbda": "10", "d": "3", "beta_intensity": "0.6", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "0"}, "params": {"datasetgeneration": "GenerationJnorm", "n": "20", "nlambda": "100", "lmbda": "10", "d": "3", "beta_intensity": "0.6", "itnum": "0", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA"}, "log": [], "threads": 1, "resources": {}, "jobid": 5901, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_20_100_10_3_0.6_NA_NA_NA_dataset0.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw ../data/assessing_models_simulation/GenerationJnorm --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules simulate_for_comparison_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw/5901.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw/5901.jobfailed; exit 1)

