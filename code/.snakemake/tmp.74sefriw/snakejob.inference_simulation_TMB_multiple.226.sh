#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_20_100_80_5_0.8_NA_NA_NA_dataset0.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm_20_100_80_5_0.8_fullREM_NA_NA_NA_dataset0.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationJnorm", "n": "20", "nlambda": "100", "lmbda": "80", "d": "5", "beta_intensity": "0.8", "model": "fullREM", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "0"}, "params": {"model": "fullREM", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationJnorm_20_100_80_5_0.8_ModelfullREM_NA_NA_NA_dataset0.log"], "threads": 1, "resources": {}, "jobid": 226, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm_20_100_80_5_0.8_fullREM_NA_NA_NA_dataset0.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw ../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_20_100_80_5_0.8_NA_NA_NA_dataset0.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw/226.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw/226.jobfailed; exit 1)

