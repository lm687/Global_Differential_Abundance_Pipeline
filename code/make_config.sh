#!/bin/bash

name_config="config_PCAWG.yaml"
name_failed_samples="../text/readme/1_faulty_samples"
echo "Creating config file $name_config ... "

echo "Samples for which there is no VCF file\nThese files appear in the metadata and may have them in the mutccf file, but I don't have their VCF, which is the only file that contains what mutation it is (in mutccf you can have the position and CCF, but not mutation type).\n\n `````` " > $name_failed_samples


## Samples
echo "samples:" > $name_config
flder="../data/restricted/pcawg/pcawg_restricted_snv/"

for i in ${flder}*gz; do
	bs=`(basename $i)`
	echo "    " ${bs::36} ": " $bs >>  $name_config
done


## Sample groups
echo "\nsample_groups:" >>  $name_config
awk '{print "- "$2}'  ../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt  | tail -n +2 | sort -u >>  $name_config


## Types of features
echo "\nfeature_types:
-    signatures
-    nucleotidesubstitution1
-    nucleotidesubstitution3" >> $name_config

## To which cancer type samples belong
printf "\n%s" "grouped_samples:" >>  $name_config
for i in `awk '{print $2}'  ../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt  | tail -n +2 | sort -u`; do  printf "\n    %s : " "$i" >> $name_config; 
flenames=`grep $i ../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt | awk '{printf "../data/restricted/pcawg/pcawg_restricted_snv_counts/"$1 " "}'`
for flename in $flenames; do
	nme="../data/restricted/pcawg/pcawg_restricted_snv/"$(basename $flename)".consensus.20160830.somatic.snv_mnv.vcf.gz"
echo $nme
## check if file exists. if it doesn't append to some extra file containing failed samples
if [ -f $nme ]; then
	printf "%s " $flename >> $name_config
else
	echo $flename >> $name_failed_samples
fi
done
done


echo "... Config file $name_config created."
echo "``````" >> > $name_failed_samples
