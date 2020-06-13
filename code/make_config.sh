#!/bin/bash

name_config="config_PCAWG.yaml"
echo "Creating config file $name_config ... "

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
for i in `awk '{print $2}'  ../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt  | tail -n +2 | sort -u`; do  printf "\n    %s : " "$i" >> $name_config;  grep $i ../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt | awk '{printf "../data/restricted/pcawg/pcawg_restricted_snv_counts/"$1 " "}' >> $name_config; done


echo "... Config file $name_config created."
