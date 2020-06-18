## Ignored folders
The data folder data/ is restricted and therefore not available

## Modifying vcfs

0. Modify VCF files to get only the important information (position, mutation type), and add the flanking bases (embed to count space)
```
sh code/0_modify_vcf/get_flanking_and_mut.sh
```

Creates files with suffix `_merged`.


1. Create ROO objects

### Signature definitions
https://www.synapse.org/#!Synapse:syn11738319






## Conda environment

```
source activate snakemake-globalDA
...
conda deactivate
```


### Config file
first of all
`config_PCAWG_working.yaml`

snakemake -p ../data/restricted/pcawg/pcawg_restricted_snv/fec30898-f86b-4207-aa78-de77142c8f50.consensus.20160830.somatic.snv_mnv.vcf.gz.tbi


## Samples in PCAWG for which there are no VCFs
These files appear in the metadata and may have them in the mutccf file, but I don't have their VCF, which is the only file that contains what mutation it is (in mutccf you can have the position and CCF, but not mutation type).

../results/reports/faulty_samples

