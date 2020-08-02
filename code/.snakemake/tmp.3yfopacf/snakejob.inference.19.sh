#!/bin/sh
# properties = {"type": "single", "rule": "inference", "local": false, "input": ["../data/roo/Kidney-RCC.clearcell_signatures_ROO.RDS"], "output": ["../data/inference/Kidney-RCC.clearcell_signatures_15000_DMROO.RData"], "wildcards": {"cancer_type": "Kidney-RCC.clearcell", "feature_type": "signatures", "nits": "15000", "model": "DM"}, "params": {"cancer_type": "Kidney-RCC.clearcell", "feature_type": "signatures", "model": ["DM"], "nits": ["15000"]}, "log": ["logs/inference/Kidney-RCC.clearcell_signatures_15000_DM.log"], "threads": 1, "resources": {}, "jobid": 19, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/inference/Kidney-RCC.clearcell_signatures_15000_DMROO.RData --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.3yfopacf ../data/roo/Kidney-RCC.clearcell_signatures_ROO.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.3yfopacf/19.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.3yfopacf/19.jobfailed; exit 1)

