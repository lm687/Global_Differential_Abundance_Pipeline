#!/bin/sh
# properties = {"type": "single", "rule": "inference", "local": false, "input": ["../data/roo/Liver-HCC_signatures_ROO.RDS"], "output": ["../data/inference/Liver-HCC_signatures_15000_DMROO.RData"], "wildcards": {"cancer_type": "Liver-HCC", "feature_type": "signatures", "nits": "15000", "model": "DM"}, "params": {"cancer_type": "Liver-HCC", "feature_type": "signatures", "model": ["DM"], "nits": ["15000"]}, "log": ["logs/inference/Liver-HCC_signatures_15000_DM.log"], "threads": 1, "resources": {}, "jobid": 21, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/inference/Liver-HCC_signatures_15000_DMROO.RData --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.3yfopacf ../data/roo/Liver-HCC_signatures_ROO.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.3yfopacf/21.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.3yfopacf/21.jobfailed; exit 1)

