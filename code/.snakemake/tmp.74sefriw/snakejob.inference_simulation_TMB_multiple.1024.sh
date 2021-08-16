#!/bin/sh
# properties = {"type": "single", "rule": "inference_simulation_TMB_multiple", "local": false, "input": ["../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_10_100_10_7_0.02_NA_NA_NA_dataset0.RDS"], "output": ["../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm_10_100_10_7_0.02_fullREDMsinglelambda_NA_NA_NA_dataset0.RDS"], "wildcards": {"optimiser": "nlminb", "datasetgeneration": "GenerationJnorm", "n": "10", "nlambda": "100", "lmbda": "10", "d": "7", "beta_intensity": "0.02", "model": "fullREDMsinglelambda", "fixed_beta_intercept": "NA", "fixed_beta_slope": "NA", "sdRE_input": "NA", "itnum": "0"}, "params": {"model": "fullREDMsinglelambda", "optimiser": "nlminb"}, "log": ["logs/inference/nlminb/simulation_GenerationJnorm_10_100_10_7_0.02_ModelfullREDMsinglelambda_NA_NA_NA_dataset0.log"], "threads": 1, "resources": {}, "jobid": 1024, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/assessing_models_simulation/inference_results/TMB/nlminb/multiple_GenerationJnorm_10_100_10_7_0.02_fullREDMsinglelambda_NA_NA_NA_dataset0.RDS --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw ../data/assessing_models_simulation/datasets/multiple_GenerationJnorm_10_100_10_7_0.02_NA_NA_NA_dataset0.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference_simulation_TMB_multiple --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw/1024.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.74sefriw/1024.jobfailed; exit 1)

