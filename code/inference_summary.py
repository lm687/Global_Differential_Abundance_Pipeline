## Read the output from slurm to see if there were divergent transitions after warmup, high R hat, low ESS, etc.

## Run: python inference_summary.py  | sort > ../data/diagnostics/inference_diagnostics.csv

import os
import re
import numpy as np

fles = [i for i in os.listdir() if "slurm-" in i]

def to_float_time(i):
	return(float(re.sub(r' seconds.*', '', re.sub(r'^.*?Elapsed Time: ', '', i), flags=re.DOTALL)))
def to_float_time2(i):
	return((re.sub(r'Chain [0-9]+:', '', re.sub(r' seconds.*', '',i, flags=re.DOTALL)).strip()))
	# (re.sub(r' seconds.*', '', re.sub(r'^.*?: ', '', i), flags=re.DOTALL))
	# (re.sub(r' seconds.*', '',i, flags=re.DOTALL))

# print("00 CT \t features \t divergent transitions \t Rhat high \t ESS")
print("*slurm,CT,features,nits,model,divergent transitions,Rhat high,ESS,Cancelled (time limit), warmup time (s), sampling time (s), total time (s)")
for fle in fles:
	# print(fle)
	list_outfile = None; #line_elapsed_time = None;
	time_limit = '.'; div_trans = '.'; Tail_ESS = '-'; Rhat = '-'
	a = open(fle, "r")
	# contents = a.read()

	# print('Elapsed' in contents)

	## search lines
	warm_up_time = []; sampling_time = []; total_time = []
	for line in a:
		if re.search("output\:", line):
			list_outfile = os.path.basename(line).rstrip("ROO.RData\n").split('_')
		if re.search("divergent transitions", line):
			div_trans = True
		if re.search("Tail-ESS", line):
			Tail_ESS = True
		if re.search("Rhat", line):
			Rhat = True
		if re.search("DUE TO TIME LIMIT", line):
			time_limit = True
		# if re.search("Elapsed Time", line):
		# 	line_elapsed_time.append(line)
		if re.search("seconds \(Warm-up\)", line):
			warm_up_time.append(line)
		if re.search("seconds \(Sampling\)", line):
			sampling_time.append(line)
		if re.search("seconds \(Total\)", line):
			total_time.append(line)


	# print([ to_float_time(i) for i in warm_up_time])
	
	warm_up_time = np.array([ to_float_time(i) for i in warm_up_time]).astype(np.float)
	# print(np.mean(warm_up_time))
	sampling_time = np.array([ to_float_time2(i) for i in sampling_time]).astype(np.float)
	# print(np.mean(sampling_time))
	total_time = np.array([ to_float_time2(i) for i in total_time]).astype(np.float)
	# print(np.mean(total_time))

	# print(total_time is None)
	# print(len(total_time))
	# if not np.isnan(total_time):
	if len(total_time > 0):
		# times_avg = '\t'.join([str(np.mean(warm_up_time)), str(np.mean(sampling_time)), str(np.mean(total_time))])
		# print('\t'.join([str(np.mean(warm_up_time)), str(np.mean(sampling_time))]))
		times_avg = ','.join([str(np.round(np.mean(warm_up_time))),\
		 str(np.round(np.mean(sampling_time))), str(np.round(np.mean(total_time)))])
	else:
		times_avg = ','.join(['-', '-', '-'])

	# print(times_avg)
	# print(a[line_elapsed_time]) ## warmup
	# print(a[line_elapsed_time+1]) ## sampling
	# print(a[line_elapsed_time+2]) ## total
	# print([x for i, x in enumerate(a) if i in line_elapsed_time])
	# print("%s \t %s \t %s \t %s" % ('\t'.join(list_outfile), div_trans, Rhat, Tail_ESS))
	print("%s,%s,%s,%s,%s,%s,%s" % (fle, ','.join(list_outfile), div_trans, Rhat, Tail_ESS, time_limit, times_avg))



