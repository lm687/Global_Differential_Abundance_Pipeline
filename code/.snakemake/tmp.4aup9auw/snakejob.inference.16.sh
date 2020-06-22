#!/bin/sh
# properties = {"type": "single", "rule": "inference", "local": false, "input": ["../data/roo/Eso-AdenoCA_signatures_ROO.RDS"], "output": ["../data/inference/Eso-AdenoCA_signatures_20000_MROO.RData"], "wildcards": {"cancer_type": "Eso-AdenoCA", "feature_type": "signatures", "nits": "20000", "model": "M"}, "params": {"cancer_type": "Eso-AdenoCA", "feature_type": "signatures", "model": ["M"], "nits": ["20000"]}, "log": ["logs/inference/Eso-AdenoCA_signatures_20000_M.log"], "threads": 1, "resources": {}, "jobid": 16, "cluster": {}}
 cd /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code && \
/home/lm687/.conda/envs/snakemake-globalDA/bin/python3.8 \
-m snakemake ../data/inference/Eso-AdenoCA_signatures_20000_MROO.RData --snakefile /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.4aup9auw ../data/roo/Eso-AdenoCA_signatures_ROO.RDS --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules inference --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.4aup9auw/16.jobfinished || (touch /rds/user/lm687/hpc-work/Global_Differential_Abundance_Pipeline/code/.snakemake/tmp.4aup9auw/16.jobfailed; exit 1)

