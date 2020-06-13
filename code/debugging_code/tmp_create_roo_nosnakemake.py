import pandas as pd

pd = pd.read_csv("tmp_config_PCAWG_part.yaml")
for f in ['nucleotidesubstitution1', 'nucleotidesubstitution3', 'signatures']:
	for i in pd['grouped_samples:']:
		splt_str = i.split(":")
		# print(splt_str)
		print("Rscript 1_create_ROO/create_ROO_split.R --input_files '{inpt}' --cancer_type {ct} --feature_type {feat} --output ../data/roo/{ct}_{feat}_ROO.RDS".format(ct = splt_str[0].strip(' '), inpt = splt_str[1].strip(' '), feat = f))
