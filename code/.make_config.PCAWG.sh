echo "samples:" >> config_PCAWG.yaml

flder="../data/restricted/pcawg/pcawg_restricted_snv/"

for i in ${flder}*gz; do
	bs=`(basename $i)`
	echo "    " ${bs::36} ": " $bs >> config_PCAWG.yaml
done





