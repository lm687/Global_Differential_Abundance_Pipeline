echo "sample_groups:" >> config_PCAWG.yaml
awk '{print "    "$2}'  ../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt  | tail -n +2 | sort -u >> config_PCAWG.yaml






