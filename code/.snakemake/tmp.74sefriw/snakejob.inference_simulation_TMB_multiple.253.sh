#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_50_100_10_3_0_NA_NA_NA_dataset0.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm_50_100_10_3_0_fullREM_NA_NA_NA_dataset0.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationJnorm", "n": "50", "nlambda": "100", "lmbda": "10", "d": "3", "beta_intensity": "0", "model": "fullREM", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "0"}, "params": {"model": "fullREM", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationJnorm_50_100_10_3_0_ModelfullREM_NA_NA_NA_dataset0.log"], "threads": 1, "resources": {}, "jobid": 253, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm_50_100_10_3_0_fullREM_NA_NA_NA_dataset0.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw ../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_50_100_10_3_0_NA_NA_NA_dataset0.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw/253.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw/253.jobfailed; exit 1)

