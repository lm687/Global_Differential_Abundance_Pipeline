vcf_filename=$1
vcf_filename_tbi=$2 ## not actually used, but since the vcf needs to be indexed first, we use this as input for snakemake to infer the dependency
ccf_file=$3
out_file_final=$4

raw_vcf_filename=`echo $vcf_filename | sed 's/.consensus.20160830.somatic.snv_mnv.vcf.gz//' | sed 's/pcawg_snv/consensus_subclonal_reconstruction_mutccf_20170325/'`

python3 0_embed_in_count_space/getflanking.py --file $vcf_filename > ${out_file_final}tmpflanking

vcf-query -f '%POS,%QUAL,%INFO/VAF\n' $vcf_filename |  sed 's/.*,//' > ${out_file_final}tmpVAF2

paste ${out_file_final}tmpflanking ${out_file_final}tmpVAF2 | column -s $'\t' -t > ${out_file_final}tmpout

python3 0_embed_in_count_space/match_ccf.py --file_VAF ${out_file_final}tmpout --file_ccf ${ccf_file} --out_file $out_file_final

rm ${out_file_final}tmpflanking
rm ${out_file_final}tmpVAF2
rm  ${out_file_final}tmpout

